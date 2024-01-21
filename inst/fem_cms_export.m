## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} fem_cms_export(@var{filename}, @var{mesh}, @var{dof_map}, @var{mat_ass}, @var{cms_opt})
## Create modal element files for MBDyn.
##
## @var{filename} @dots{} Two output files will be created: "<@var{filename}>.elm" and "<@var{filename}>.fem"
## "<@var{filename}>.elm" must be included in the main input file of MBDyn via include: "<@var{filename}>.elm";
##
## @var{mesh} @dots{} Finite Element mesh data structure.
##
## @var{dof_map} @dots{} Degree of freedom mapping
##
## @var{mat_ass} @dots{} Reduced order Finite Element matrices returned from fem_cms_create.
##
## @var{cms_opt} @dots{} Substructure options returned from fem_cms_create.
## @seealso{fem_cms_create}
## @end deftypefn

function fem_cms_export(filename, mesh, dof_map, mat_ass, cms_opt)
  if (nargin ~= 5 || nargout > 0)
    print_usage();
  endif

  if (~isfield(cms_opt, "verbose"))
    cms_opt.verbose = false;
  endif

  if (~isfield(cms_opt, "invariants"))
    cms_opt.invariants = true;
  endif

  if (~isfield(cms_opt, "create_binary"))
    cms_opt.create_binary = false;
  endif

  if (~isfield(cms_opt, "use_binary"))
    cms_opt.use_binary = false;
  endif

  if (~isfield(cms_opt, "update_binary"))
    cms_opt.update_binary = false;
  endif

  if (ispc())
    filename(find(filename == '\')) = '/'; ## Required for MBDyn's parser
  endif

  fd = -1;

  unwind_protect
    [fd, msg] = fopen([filename, ".elm"], "w");

    if (fd == -1)
      error("failed to open file \"%s\": %s", [filename, ".elm"], msg);
    endif

    warning("error", "Octave:singular-matrix", "local");

    fprintf(fd, "## cond(Mred)=%.1e\n", cond(mat_ass.Mred));
    fprintf(fd, "## cond(Kred)=%.1e\n\n", cond(mat_ass.Kred));

    fprintf(fd, "joint: %s, modal, %s,\n", cms_opt.element.name, cms_opt.nodes.modal.name);

    if (isfield(cms_opt, "selected_modes"))
      fprintf(fd, "\t%d,\n", numel(cms_opt.selected_modes));
      fprintf(fd, "\tlist");
      fprintf(fd, ", %d", cms_opt.selected_modes);
      fprintf(fd, ",\n");
    else
      fprintf(fd, "\t%d,\n", columns(mat_ass.Tred));
    endif

    fprintf(fd, "\tfrom file,\n");
    fprintf(fd, "\tdamping from file,\n");
    fprintf(fd, "\t\"%s\",\n", [filename, ".fem"]);

    if (cms_opt.create_binary)
      fputs(fd, "\tcreate binary,\n");
    endif

    if (cms_opt.use_binary)
      fputs(fd, "\tuse binary,\n");
    endif

    if (cms_opt.update_binary)
      fputs(fd, "\tupdate binary,\n");
    endif

    fprintf(fd, "\torigin node, %d,\n", cms_opt.nodes.modal.number);

    inum_nodes_itf = int32(0);

    for i=1:numel(cms_opt.nodes.interfaces)
      if (numel(cms_opt.nodes.interfaces(i).name))
        ++inum_nodes_itf;
      endif
    endfor

    fprintf(fd, "\t%d", inum_nodes_itf);

    for i=1:numel(cms_opt.nodes.interfaces)
      if (numel(cms_opt.nodes.interfaces(i).name))
        fprintf(fd, ",\n\t\t%d, %s, null", cms_opt.nodes.interfaces(i).number, cms_opt.nodes.interfaces(i).name);
      endif
    endfor

    if (~cms_opt.invariants)
      fputs(fd, ",\n\tuse invariant 9"); ## MBDyn will not compute Inv9 by default
    endif

    fprintf(fd, ";\n");
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect

  if (cms_opt.verbose)
    fprintf(stderr, "creating file \"%s\" ...\n", [filename, ".fem"]);
    tic();
  endif

  if (cms_opt.invariants)
    idx_node_output = [cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
  else
    idx_node_output = 1:rows(mesh.nodes);
  endif

  Phi = zeros(6 * numel(idx_node_output), columns(mat_ass.Tred));
  ridx_node = zeros(max(dof_map.idx_node), 1, "int32");
  ridx_node(dof_map.idx_node) = 1:numel(dof_map.idx_node);

  for i=1:columns(dof_map.ndof)
    idx_dof_output = dof_map.ndof(idx_node_output, i);
    idx_act_dof = find(idx_dof_output > 0);
    idx_glob_dof = idx_dof_output(idx_act_dof);
    Phi((idx_act_dof - 1) * 6 + i, :) = mat_ass.Tred(ridx_node(idx_glob_dof), :);
  endfor

  if (isfield(mat_ass, "KTAU0red"))
    KTAU0red = mat_ass.KTAU0red;
    index_KTAU0red = cms_opt.index_KTAU0red;
  else
    KTAU0red = [];
    index_KTAU0red = [];
  endif

  for i=1:numel(index_KTAU0red)
    if (~(index_KTAU0red(i) >= 1 && index_KTAU0red(i) <= 12 + numel(cms_opt.nodes.interfaces) * 6))
      error("invalid index for cms_opt.index_KTAU0red(%d)=%d", i, index_KTAU0red(i));
    endif
  endfor

  if (cms_opt.invariants)
    ## position of the modal node
    X0 = mesh.nodes(cms_opt.nodes.modal.number, 1:3).';

    ## centre of gravity with respect to the global FEM reference frame
    Xgc = mat_ass.S / mat_ass.dm;

    ## moment of inertia with respect to the center of gravity
    Jgc = mat_ass.J + (skew(Xgc) * skew(Xgc)) * mat_ass.dm;

    mbdyn_pre_write_fem_data([filename, ".fem"], ...
                             mat_ass.Mred, ...
                             mat_ass.Dred, ...
                             mat_ass.Kred, ...
                             Phi, ...
                             mesh.nodes(idx_node_output, 1:3).', ...
                             zeros(columns(mat_ass.Mred), 1), ...
                             zeros(columns(mat_ass.Mred), 1), ...
                             [], ...
                             mat_ass.dm, ...
                             Xgc, ...
                             Jgc, ...
                             idx_node_output, ...
                             mat_ass.Inv3, ...
                             mat_ass.Inv4, ...
                             mat_ass.Inv5, ...
                             mat_ass.Inv8, ...
                             mat_ass.Inv9, ...
                             index_KTAU0red, ...
                             KTAU0red);
  else
    mbdyn_pre_write_fem_data([filename, ".fem"], ...
                             mat_ass.Mred, ...
                             mat_ass.Dred, ...
                             mat_ass.Kred, ...
                             Phi, ...
                             mesh.nodes(:, 1:3).', ...
                             zeros(columns(mat_ass.Mred), 1), ...
                             zeros(columns(mat_ass.Mred), 1), ...
                             mat_ass.diagM, ...
                             [], ...
                             [], ...
                             [], ...
                             [], ...
                             [], ...
                             [], ...
                             [], ...
                             [], ...
                             [], ...
                             index_KTAU0red, ...
                             KTAU0red);
  endif

  if (cms_opt.verbose)
    fprintf(stderr, "file \"%s\" created ...\n", [filename, ".fem"]);
    toc();
  endif
endfunction

