## Copyright (C) 2019(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{sol_dyn} = fem_post_cms_sol_import(@var{output_file}, @var{cms_data})
## Import modal element solutions from MBDyn
##
## @var{output_file} @dots{} MBDyn output filename
##
## @var{cms_data} @dots{} Finite element data including mesh, options, and mode shapes used to create modal element data for MBDyn
##
## @var{sol_dyn} @dots{} Modal solution imported from MBDyn
##
## @seealso{fem_cms_export, fem_post_cms_expand, fem_post_cms_sol_merge, fem_post_mesh_merge, fem_post_mesh_export}
## @end deftypefn

function sol_dyn = fem_post_cms_sol_import(output_file, cms_data)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif

  log_dat = mbdyn_post_load_log(output_file);

  [t, TStep, NIter, ResErr, SolErr, SolConv, OutputFlag] = mbdyn_post_load_output_out(output_file, 1024, false);

  t = t(find(OutputFlag));

  [elem_id, q] = mbdyn_post_load_output_mod(output_file, numel(t));

  istr_node_idx_o = zeros(1, numel(cms_data), "int32");
  istr_node_idx_l = zeros(1, numel(cms_data), "int32");
  imodal_node_label = zeros(1, numel(cms_data), "int32");
  ielem_label = zeros(1, numel(cms_data), "int32");
  ielem_idx = zeros(1, numel(cms_data), "int32");

  sol_dyn.t = t;

  empty_cell = cell(1, numel(cms_data));

  sol_dyn.bodies = struct("X", empty_cell, "R", empty_cell, "q", empty_cell);

  for i=1:numel(cms_data)
    imodal_node_label(i) = getfield(log_dat.vars, cms_data(i).cms_opt.nodes.modal.name);
    ielem_label(i) = getfield(log_dat.vars, cms_data(i).cms_opt.element.name);
    ielem_idx_tmp = find(elem_id == ielem_label(i));

    if (numel(ielem_idx_tmp) ~= 1)
      error("modal element %d(%s) not present in file \"%s\"", ...
            ielem_label(i), ...
            cms_data(i).cms_opt.element.name, ...
            output_file);
    endif

    ielem_idx(i) = ielem_idx_tmp;

    istr_node_idx_l_tmp = find(imodal_node_label(i) == [log_dat.nodes.label]);

    if (numel(istr_node_idx_l_tmp) ~= 1)
            error("modal node %d(%s) not present in log file \"%s\"", ...
            imodal_node_label(i), ...
            cms_data(i).cms_opt.nodes.modal.name, ...
            output_file);
    endif

    istr_node_idx_l(i) = istr_node_idx_l_tmp;
  endfor

  [str_node_id, trajectory] = mbdyn_post_load_output_mov(output_file, imodal_node_label, numel(t));

  for i=1:numel(cms_data)
    istr_node_idx_o_tmp = find(imodal_node_label(i) == str_node_id);

    if (numel(istr_node_idx_o_tmp) ~= 1)
      error("modal node %d(%s) not present in file \"%s\"", ...
            imodal_node_label(i), ...
            cms_data(i).cms_opt.nodes.modal.name, ...
            output_file);
    endif

    istr_node_idx_o(i) = istr_node_idx_o_tmp;
  endfor

  for i=1:numel(cms_data)
    sol_dyn.bodies(i).nodes.modal.label = imodal_node_label(i);
    sol_dyn.bodies(i).nodes.modal.name = cms_data(i).cms_opt.nodes.modal.name;
    sol_dyn.bodies(i).element.label = ielem_label(i);
    sol_dyn.bodies(i).element.name = cms_data(i).cms_opt.element.name;

    sol_dyn.bodies(i).X0 = log_dat.nodes(istr_node_idx_l(i)).X0;
    sol_dyn.bodies(i).R0 = log_dat.nodes(istr_node_idx_l(i)).R0;
    sol_dyn.bodies(i).X = trajectory{istr_node_idx_o(i)}(:, 1:3).';

    switch (log_dat.nodes(istr_node_idx_l(i)).orientation_description)
      case "euler123"
        rotfunc = @euler123_to_rotation_matrix;
      case "euler313"
        rotfunc = @euler313_to_rotation_matrix;
      case "euler321"
        rotfunc = @euler321_to_rotation_matrix;
      case "phi"
        rotfunc = @rotation_vector_to_rotation_matrix;
      otherwise
        error("orientation description \"%s\" not supported", ...
              log_dat.nodes(istr_node_idx_l(i)).orientation_description);
    endswitch

    sol_dyn.bodies(i).R = feval(rotfunc, trajectory{istr_node_idx_o(i)}(:, 4:6).');
    sol_dyn.bodies(i).q = q{ielem_idx(i)};
  endfor
endfunction
