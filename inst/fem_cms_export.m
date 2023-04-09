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

%!test
%! ## TEST 1
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;

%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 1e-7 / (SI_unit_second);

%! param.h = 2.5e-3 / SI_unit_meter; ## mesh size
%! param.l = 350e-3 / SI_unit_meter; ## bearing distance
%! param.d = 10e-3 / SI_unit_meter; ## shaft diameter
%! param.D = 150e-3 / SI_unit_meter; ##disk diameter
%! param.w = 15e-3 / SI_unit_meter; ## disk width
%! param.o = 75e-3 / SI_unit_meter; ## disk offset

%! param.ecg = 1e-3 * param.D;
%! param.dm = 1e-6 / SI_unit_kilogram;
%! param.E = 210000e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! m1 = param.D^2 * pi / 4 * param.w * param.rho;
%! param.dr = param.ecg * (m1 + param.dm) / param.dm;
%! param.omega0 = 1000 * pi / 30 / (1 / SI_unit_second);
%! param.omega1 = 10000 * pi / 30 / (1 / SI_unit_second);
%! param.n = 10;

%! cms_opt.algorithm = "shift-invert";
%! cms_opt.refine_max_iter = int32(3);
%! cms_opt.element.name = "elem_id_rotor";
%! cms_opt.nodes.modal.name = "node_id_rotor";
%! cms_opt.nodes.interfaces(1).name = "node_id_bearing1";
%! cms_opt.nodes.interfaces(2).name = "node_id_bearing2";
%! cms_opt.number_of_threads = int32(4);
%! cms_opt.verbose = false;
%! options.number_of_modes = int32(10);
%! options.scale_def = 10e-3;
%! options.geo_tol = sqrt(eps);
%! options.code.use_package = false;
%! options.mbdyn_command = "mbdyn";
%! options.f_run_mbdyn = [true, true];
%! options.verbose = false;

%! function R2j = norm_ref_frame(R2)
%!   e1 = R2(:, 1);
%!   e3 = [0; 0; 1];
%!   e2 = cross(e3, e1);
%!   e1 = cross(e2, e3);
%!   R2j = [e1 / norm(e1), e2 / norm(e2), e3 / norm(e3)];
%! endfunction
%!
%! filename = "";
%!
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   mbdyn_filename_suffix = {"cms", "beam"};
%!   for i=1:numel(mbdyn_filename_suffix)
%!     mbdyn_filenames{i} = [filename, mbdyn_filename_suffix{i}, ".mbdyn"];
%!   endfor
%!
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen(mbdyn_filenames{1}, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_filenames{1});
%!     endif

%!     fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!     fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!     fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!     fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!     fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!     fputs(fd, "set: integer node_id_rotor = 2002;\n");
%!     fputs(fd, "set: integer node_id_bearing1 = 2003;\n");
%!     fputs(fd, "set: integer node_id_bearing2 = 2004;\n");
%!     fputs(fd, "set: integer body_id_unbalance = 3000;\n");
%!     fputs(fd, "set: integer elem_id_rotor = 3005;\n");
%!     fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!     fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!     fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!     fputs(fd, "set: integer elem_id_inertia = 3010;\n");
%!     fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!     fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!     fputs(fd, "set: real initial_time = 0.;\n");
%!     fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: initial_time;\n");
%!     fputs(fd, "        final time: final_time;\n");
%!     fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!     fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!     fputs(fd, "        method: ms, 0.6;\n");
%!     fputs(fd, "        tolerance: 1e-8;\n");
%!     fputs(fd, "        max iterations: 100;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: auto;\n");
%!     fputs(fd, "        output: iterations;\n");
%!     fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always, max iterations, 100;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "        structural nodes:\n");
%!     fputs(fd, "                +1		# modal\n");
%!     fputs(fd, "                +1		# interface 1\n");
%!     fputs(fd, "                +1              # interface 2\n");
%!     fputs(fd, "        ;\n");
%!     fputs(fd, "        joints:\n");
%!     fputs(fd, "                +1		# modal\n");
%!     fputs(fd, "                +1		# bearing1\n");
%!     fputs(fd, "                +1		# bearing2\n");
%!     fputs(fd, "                +1              # drive\n");
%!     fputs(fd, "        ;\n");
%!     fputs(fd, "        rigid bodies: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_shaft,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!     fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null;\n");
%!     fputs(fd, "reference: ref_id_bearing1,\n");
%!     fputs(fd, "        reference, ref_id_shaft,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   -0.5 * l,\n");
%!     fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null;\n");
%!     fputs(fd, "reference: ref_id_bearing2,\n");
%!     fputs(fd, "        reference, ref_id_shaft,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.5 * l,\n");
%!     fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, modal,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_bearing1, static,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, null;\n");
%!     fputs(fd, "        structural: node_id_bearing2, static,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\";\n");
%!     fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!     fputs(fd, "       body: body_id_unbalance,\n");
%!     fputs(fd, "             node_id_rotor,\n");
%!     fputs(fd, "                dm,\n");
%!     fputs(fd, "                reference, node, dr, 0., 0.,\n");
%!     fputs(fd, "                diag, 0., 0., 0.;\n");
%!     fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!     fputs(fd, "                node_id_bearing1,\n");
%!     fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "               position constraint,\n");
%!     fputs(fd, "                        active, active, active,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "               orientation constraint,\n");
%!     fputs(fd, "                        inactive, inactive, inactive,\n");
%!     fputs(fd, "                        null;\n");
%!     fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!     fputs(fd, "               node_id_bearing2,\n");
%!     fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "               position constraint,\n");
%!     fputs(fd, "                        active, active, inactive,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "               orientation constraint,\n");
%!     fputs(fd, "                        inactive, inactive, inactive,\n");
%!     fputs(fd, "                        null;\n");
%!     fputs(fd, "	joint: joint_id_drive, angular velocity,\n");
%!     fputs(fd, "		# node label\n");
%!     fputs(fd, "		node_id_rotor, \n");
%!     fputs(fd, "		# direction\n");
%!     fputs(fd, "		0.,0.,1.,\n");
%!     fputs(fd, "		# angular velocity\n");
%!     fputs(fd, "		reference, drive_id_rotor_speed;\n");
%!     fputs(fd, "        include: \"${MBDYN_ROTOR_DYN_CMS_ELEM_FILE}\";\n");
%!     fputs(fd, "        inertia: elem_id_inertia,\n");
%!     fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!     fputs(fd, "                 output, both;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect

%!   fd = -1;
%!   unwind_protect
%!     fd = fopen(mbdyn_filenames{2}, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_filenames{2});
%!     endif
%!     fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!     fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!     fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!     fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!     fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!     fputs(fd, "set: integer ref_id_shaft_section = 1006;\n");
%!     fputs(fd, "set: integer ref_id_bearing1_center = 1007;\n");
%!     fputs(fd, "set: integer ref_id_bearing2_center = 1008;\n");
%!     fputs(fd, "set: integer node_id_rotor = 2001;\n");
%!     fputs(fd, "set: integer node_id_bearing1 = 2002;\n");
%!     fputs(fd, "set: integer node_id_bearing1_center = 2003;\n");
%!     fputs(fd, "set: integer node_id_bearing2 = 2004;\n");
%!     fputs(fd, "set: integer node_id_bearing2_center = 2005;\n");
%!     fputs(fd, "set: integer body_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer body_id_bearing1 = 3002;\n");
%!     fputs(fd, "set: integer body_id_bearing1_center = 3003;\n");
%!     fputs(fd, "set: integer body_id_bearing2 = 3004;\n");
%!     fputs(fd, "set: integer body_id_bearing2_center = 3005;\n");
%!     fputs(fd, "set: integer beam_id_shaft1 = 4001;\n");
%!     fputs(fd, "set: integer beam_id_shaft2 = 4002;\n");
%!     fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!     fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!     fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!     fputs(fd, "set: integer elem_id_inertia = 4001;\n");
%!     fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!     fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!     fputs(fd, "set: real l1 = 0.5 * l + o - 0.5 * w;\n");
%!     fputs(fd, "set: real l2 = 0.5 * l - o - 0.5 * w;\n");
%!     fputs(fd, "set: real rotor_m = rho * D^2 * pi / 4. * w;\n");
%!     fputs(fd, "set: real rotor_Jx = rotor_m * ((0.5 * D)^2 + w^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real rotor_Jy = rotor_Jx;\n");
%!     fputs(fd, "set: real rotor_Jz = rotor_m * (0.5 * D)^2 / 2.;\n");
%!     fputs(fd, "set: real shaft_rotor1_m = rho * d^2 * pi / 4. * (l1 / 4.);\n");
%!     fputs(fd, "set: real shaft_rotor1_Jx = shaft_rotor1_m * ((0.5 * d)^2 + (l1 / 4.)^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real shaft_rotor1_Jy = shaft_rotor1_Jx;\n");
%!     fputs(fd, "set: real shaft_rotor1_Jz = shaft_rotor1_m * (0.5 * d)^2 / 2.;\n");
%!     fputs(fd, "set: real shaft_rotor2_m = rho * d^2 * pi / 4. * (l2 / 4.);\n");
%!     fputs(fd, "set: real shaft_rotor2_Jx = shaft_rotor2_m * ((0.5 * d)^2 + (l2 / 4.)^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real shaft_rotor2_Jy = shaft_rotor2_Jx;\n");
%!     fputs(fd, "set: real shaft_rotor2_Jz = shaft_rotor2_m * (0.5 * d)^2 / 2.;\n");
%!     fputs(fd, "set: real shaft_bearing1_m = rho * d^2 * pi / 4. * (l1 / 4.);\n");
%!     fputs(fd, "set: real shaft_bearing1_Jx = shaft_bearing1_m * ((0.5 * d)^2 + (l1 / 4.)^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real shaft_bearing1_Jy = shaft_bearing1_Jx;\n");
%!     fputs(fd, "set: real shaft_bearing1_Jz = shaft_bearing1_m * (0.5 * d)^2 / 2.;\n");
%!     fputs(fd, "set: real shaft_bearing1_center_m = rho * d^2 * pi / 4. * (l1 / 2.);\n");
%!     fputs(fd, "set: real shaft_bearing1_center_Jx = shaft_bearing1_center_m * ((0.5 * d)^2 + (l1 / 2.)^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real shaft_bearing1_center_Jy = shaft_bearing1_center_Jx;\n");
%!     fputs(fd, "set: real shaft_bearing1_center_Jz = shaft_bearing1_center_m * ((0.5 * d)^2) / 2.;\n");
%!     fputs(fd, "set: real shaft_bearing2_m = rho * d^2 * pi / 4. * (l2 / 4.);\n");
%!     fputs(fd, "set: real shaft_bearing2_Jx = shaft_bearing2_m * ((0.5 * d)^2 + (l2 / 4.)^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real shaft_bearing2_Jy = shaft_bearing2_Jx;\n");
%!     fputs(fd, "set: real shaft_bearing2_Jz = shaft_bearing2_m * (0.5 * d)^2 / 2.;\n");
%!     fputs(fd, "set: real shaft_bearing2_center_m = rho * d^2 * pi / 4. * (l2 / 2.);\n");
%!     fputs(fd, "set: real shaft_bearing2_center_Jx = shaft_bearing2_center_m * ((0.5 * d)^2 + (l2 / 2.)^2 / 3.) / 4.;\n");
%!     fputs(fd, "set: real shaft_bearing2_center_Jy = shaft_bearing2_center_Jx;\n");
%!     fputs(fd, "set: real shaft_bearing2_center_Jz = shaft_bearing2_center_m * ((0.5 * d)^2) / 2.;\n");
%!     fputs(fd, "set: real shaft_A = d^2 * pi / 4;\n");
%!     fputs(fd, "set: real shaft_As = 9. / 10. * shaft_A;\n");
%!     fputs(fd, "set: real shaft_Iy = d^4 * pi / 64.;\n");
%!     fputs(fd, "set: real shaft_Iz = shaft_Iy;\n");
%!     fputs(fd, "set: real shaft_Ip = shaft_Iy + shaft_Iz;\n");
%!     fputs(fd, "set: real shaft_It = shaft_Ip;\n");
%!     fputs(fd, "set: real shaft_E = E;\n");
%!     fputs(fd, "set: real shaft_nu = nu;\n");
%!     fputs(fd, "set: real shaft_G = shaft_E / (2 * (1 + shaft_nu));\n");
%!     fputs(fd, "set: real shaft_rho = rho;\n");
%!     fputs(fd, "set: real shaft_damping_ratio = beta;\n");
%!     fputs(fd, "set: real initial_time = 0.;\n");
%!     fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: initial_time;\n");
%!     fputs(fd, "        final time: final_time;\n");
%!     fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!     fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!     fputs(fd, "        method: ms, 0.6;\n");
%!     fputs(fd, "        tolerance: 1e-6;\n");
%!     fputs(fd, "        max iterations: 100;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: auto;\n");
%!     fputs(fd, "        output: iterations;\n");
%!     fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always, max iterations, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "        structural nodes: 5;\n");
%!     fputs(fd, "        joints: 3;\n");
%!     fputs(fd, "        beams: 2;\n");
%!     fputs(fd, "        rigid bodies: 5;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_shaft,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!     fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null;\n");
%!     fputs(fd, "reference: ref_id_shaft_section,\n");
%!     fputs(fd, "           reference, ref_id_shaft, null,\n");
%!     fputs(fd, "           reference, ref_id_shaft, 1, 0., 0., 1.,\n");
%!     fputs(fd, "                                    2, 0., 1., 0.,\n");
%!     fputs(fd, "           reference, ref_id_shaft, null,\n");
%!     fputs(fd, "           reference, ref_id_shaft, null;\n");
%!     fputs(fd, "reference: ref_id_bearing1,\n");
%!     fputs(fd, "        reference, ref_id_shaft,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   -0.5 * l,\n");
%!     fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null;\n");
%!     fputs(fd, "reference: ref_id_bearing1_center,\n");
%!     fputs(fd, "        reference, ref_id_bearing1,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.5 * l1,\n");
%!     fputs(fd, "        reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "        reference, ref_id_rotor, null,\n");
%!     fputs(fd, "        reference, ref_id_rotor, null;\n");
%!     fputs(fd, "reference: ref_id_bearing2,\n");
%!     fputs(fd, "        reference, ref_id_shaft,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.5 * l,\n");
%!     fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null,\n");
%!     fputs(fd, "        reference, ref_id_shaft, null;\n");
%!     fputs(fd, "reference: ref_id_bearing2_center,\n");
%!     fputs(fd, "        reference, ref_id_bearing2,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   0.,\n");
%!     fputs(fd, "                   -0.5 * l2,\n");
%!     fputs(fd, "        reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "        reference, ref_id_rotor, null,\n");
%!     fputs(fd, "        reference, ref_id_rotor, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, dynamic,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_bearing1, dynamic,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing1, null;\n");
%!     fputs(fd, "        structural: node_id_bearing1_center, dynamic,\n");
%!     fputs(fd, "                reference, ref_id_bearing1_center, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing1_center, eye,\n");
%!     fputs(fd, "                reference, ref_id_bearing1_center, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing1_center, null;\n");
%!     fputs(fd, "        structural: node_id_bearing2, dynamic,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing2, null;\n");
%!     fputs(fd, "        structural: node_id_bearing2_center, dynamic,\n");
%!     fputs(fd, "                reference, ref_id_bearing2_center, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing2_center, eye,\n");
%!     fputs(fd, "                reference, ref_id_bearing2_center, null,\n");
%!     fputs(fd, "                reference, ref_id_bearing2_center, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\";\n");
%!     fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!     fputs(fd, "       body: body_id_rotor,\n");
%!     fputs(fd, "             node_id_rotor,\n");
%!     fputs(fd, "             condense, 4,\n");
%!     fputs(fd, "                rotor_m,\n");
%!     fputs(fd, "                  reference, node, null,\n");
%!     fputs(fd, "                  diag, rotor_Jx, rotor_Jy, rotor_Jz,\n");
%!     fputs(fd, "                shaft_rotor1_m,\n");
%!     fputs(fd, "                  reference, node, 0., 0., -w / 2. - l1 / 8.,\n");
%!     fputs(fd, "                  diag, shaft_rotor1_Jx, shaft_rotor1_Jy, shaft_rotor1_Jz,\n");
%!     fputs(fd, "                shaft_rotor2_m,\n");
%!     fputs(fd, "                  reference, node, 0., 0., w / 2. + l2 / 8.,\n");
%!     fputs(fd, "                  diag, shaft_rotor2_Jx, shaft_rotor2_Jy, shaft_rotor2_Jz,                \n");
%!     fputs(fd, "                dm,\n");
%!     fputs(fd, "                  reference, node, dr, 0., 0.,\n");
%!     fputs(fd, "                  diag, 0., 0., 0.;\n");
%!     fputs(fd, "        body: body_id_bearing1,\n");
%!     fputs(fd, "              node_id_bearing1,\n");
%!     fputs(fd, "              shaft_bearing1_m,\n");
%!     fputs(fd, "              reference, node, 0., 0., l1 / 8.,\n");
%!     fputs(fd, "              diag, shaft_bearing1_Jx, shaft_bearing1_Jy, shaft_bearing1_Jz;\n");
%!     fputs(fd, "        body: body_id_bearing1_center,\n");
%!     fputs(fd, "              node_id_bearing1_center,\n");
%!     fputs(fd, "              shaft_bearing1_center_m,\n");
%!     fputs(fd, "              reference, node, null,\n");
%!     fputs(fd, "              diag, shaft_bearing1_center_Jx, shaft_bearing1_center_Jy, shaft_bearing1_center_Jz;\n");
%!     fputs(fd, "        body: body_id_bearing2,\n");
%!     fputs(fd, "              node_id_bearing2,\n");
%!     fputs(fd, "              shaft_bearing2_m,\n");
%!     fputs(fd, "              reference, node, 0., 0., -l2 / 8.,\n");
%!     fputs(fd, "              diag, shaft_bearing2_Jx, shaft_bearing2_Jy, shaft_bearing2_Jz;\n");
%!     fputs(fd, "        body: body_id_bearing2_center,\n");
%!     fputs(fd, "              node_id_bearing2_center,\n");
%!     fputs(fd, "              shaft_bearing2_center_m,\n");
%!     fputs(fd, "              reference, node, null,\n");
%!     fputs(fd, "              diag, shaft_bearing2_center_Jx, shaft_bearing2_center_Jy, shaft_bearing2_center_Jz;\n");
%!     fputs(fd, "    beam3: beam_id_shaft1,\n");
%!     fputs(fd, "                # node 1\n");
%!     fputs(fd, "                node_id_bearing1, position, reference, node, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # node 2\n");
%!     fputs(fd, "                node_id_bearing1_center, position, reference, node, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # node 3,\n");
%!     fputs(fd, "                node_id_rotor, position, reference, node, 0., 0., -0.5 * w,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # orientation matrix section I\n");
%!     fputs(fd, "                reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # constitutive law section I\n");
%!     fputs(fd, "                linear viscoelastic generic,\n");
%!     fputs(fd, "                diag, shaft_E * shaft_A , shaft_G * shaft_As, shaft_G * shaft_As,\n");
%!     fputs(fd, "                      shaft_G * shaft_It, shaft_E * shaft_Iy, shaft_E * shaft_Iz,\n");
%!     fputs(fd, "                proportional, shaft_damping_ratio,\n");
%!     fputs(fd, "                # orientation matrix section II\n");
%!     fputs(fd, "                same,\n");
%!     fputs(fd, "                # constitutive law section II\n");
%!     fputs(fd, "                same;\n");
%!     fputs(fd, "    beam3: beam_id_shaft2,\n");
%!     fputs(fd, "                # node 1\n");
%!     fputs(fd, "                node_id_rotor, position, reference, node, 0., 0., 0.5 * w,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,                \n");
%!     fputs(fd, "                # node 2\n");
%!     fputs(fd, "                node_id_bearing2_center, position, reference, node, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # node 3,\n");
%!     fputs(fd, "                node_id_bearing2, position, reference, node, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # orientation matrix section I\n");
%!     fputs(fd, "                reference, ref_id_shaft_section, eye,\n");
%!     fputs(fd, "                # constitutive law section I\n");
%!     fputs(fd, "                linear viscoelastic generic,\n");
%!     fputs(fd, "                diag, shaft_E * shaft_A , shaft_G * shaft_As, shaft_G * shaft_As,\n");
%!     fputs(fd, "                      shaft_G * shaft_It, shaft_E * shaft_Iy, shaft_E * shaft_Iz,\n");
%!     fputs(fd, "                proportional, shaft_damping_ratio,\n");
%!     fputs(fd, "                # orientation matrix section II\n");
%!     fputs(fd, "                same,\n");
%!     fputs(fd, "                # constitutive law section II\n");
%!     fputs(fd, "                same;\n");
%!     fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!     fputs(fd, "                node_id_bearing1,\n");
%!     fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!     fputs(fd, "               position constraint,\n");
%!     fputs(fd, "                        active, active, active,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "               orientation constraint,\n");
%!     fputs(fd, "                        inactive, inactive, inactive,\n");
%!     fputs(fd, "                        null;\n");
%!     fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!     fputs(fd, "               node_id_bearing2,\n");
%!     fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!     fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!     fputs(fd, "               position constraint,\n");
%!     fputs(fd, "                        active, active, inactive,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "               orientation constraint,\n");
%!     fputs(fd, "                        inactive, inactive, inactive,\n");
%!     fputs(fd, "                        null;\n");
%!     fputs(fd, "	joint: joint_id_drive, angular velocity,\n");
%!     fputs(fd, "		# node label\n");
%!     fputs(fd, "		node_id_rotor, \n");
%!     fputs(fd, "		# direction\n");
%!     fputs(fd, "		0.,0.,1.,\n");
%!     fputs(fd, "		# angular velocity\n");
%!     fputs(fd, "		reference, drive_id_rotor_speed;\n");
%!     fputs(fd, "        inertia: elem_id_inertia,\n");
%!     fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!     fputs(fd, "                 output, both;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   geometry_file = [filename, ".geo"];
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");

%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif

%!     fn = fieldnames(param);

%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %g;\n", fn{i}, getfield(param, fn{i}));
%!     endfor

%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, -0.5 * l};\n");
%!     fputs(fd, "Point(2) = {0.5 * d, 0, -0.5 * l};\n");
%!     fputs(fd, "Point(3) = {0.5 * d, 0, -0.5 * w + o};\n");
%!     fputs(fd, "Point(4) = {0.5 * D, 0, -0.5 * w + o};\n");
%!     fputs(fd, "Point(5) = {0.5 * D, 0, 0.5 * w + o};\n");
%!     fputs(fd, "Point(6) = {0.5 * d, 0, 0.5 * w + o};\n");
%!     fputs(fd, "Point(7) = {0.5 * d, 0, 0.5 * l};\n");
%!     fputs(fd, "Point(8) = {0, 0, 0.5 * l};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "vol1[] = Extrude{{0, 0, 1},{0, 0, 0}, 2 * Pi}{ Surface{1}; };\n");
%!     fputs(fd, "B[] = Unique(Abs(Boundary{Volume{vol1};}));\n");
%!     fputs(fd, "Physical Volume(\"V\", 1) = {vol1};\n");
%!     fputs(fd, "For i In {0:#B[] - 1}\n");
%!     fputs(fd, "    Physical Surface(i) = {B[i]};\n");
%!     fputs(fd, "EndFor\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        "-clmin", sprintf("%g", 0.75 * param.h), ...
%!                        "-clmax", sprintf("%g", 1.25 * param.h), ...
%!                        "-ho_min", "0.1", ...
%!                        "-ho_max", "3", ...
%!                        geometry_file, ...
%!                        "-o", [filename, ".msh"]});

%!   status = spawn_wait(pid);

%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif

%!   fprintf(stderr, "loading mesh ...\n");
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");

%!   group_defs(1).id = 1;
%!   group_defs(1).name = "bearing1";
%!   group_defs(1).R = eye(3);
%!   group_defs(1).X0 = [0; 0; -0.5 * param.l];
%!   group_defs(1).type = "box";
%!   group_defs(1).geometry.xmin = -0.5 * param.d;
%!   group_defs(1).geometry.xmax = 0.5 * param.d;
%!   group_defs(1).geometry.ymin = -0.5 * param.d;
%!   group_defs(1).geometry.ymax = 0.5 * param.d;
%!   group_defs(1).geometry.zmin = 0;
%!   group_defs(1).geometry.zmax = 0;

%!   group_defs(2).id = 2;
%!   group_defs(2).name = "bearing2";
%!   group_defs(2).R = eye(3);
%!   group_defs(2).X0 = [0; 0; 0.5 * param.l];
%!   group_defs(2).type = "box";
%!   group_defs(2).geometry.xmin = -0.5 * param.d;
%!   group_defs(2).geometry.xmax = 0.5 * param.d;
%!   group_defs(2).geometry.ymin = -0.5 * param.d;
%!   group_defs(2).geometry.ymax = 0.5 * param.d;
%!   group_defs(2).geometry.zmin = 0;
%!   group_defs(2).geometry.zmax = 0;

%!   group_defs(3).id = 3;
%!   group_defs(3).name = "rotor";
%!   group_defs(3).R = eye(3);
%!   group_defs(3).X0 = [0; 0; param.o];
%!   group_defs(3).type = "cylinder";
%!   group_defs(3).geometry.rmin = 0.5 * param.D;
%!   group_defs(3).geometry.rmax = 0.5 * param.D;
%!   group_defs(3).geometry.zmin = -0.5 * param.w;
%!   group_defs(3).geometry.zmax = 0.5 * param.w;

%!   groups = fem_pre_mesh_groups_create(mesh, group_defs, options.geo_tol);

%!   mesh.groups.tria6 = groups.tria6;

%!   cms_opt.modes.number = int32(options.number_of_modes);
%!   cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!   cms_opt.nodes.interfaces(1).number = int32(rows(mesh.nodes) + 2);
%!   cms_opt.nodes.interfaces(2).number = int32(rows(mesh.nodes) + 3);
%!   cms_opt.invariants = true;
%!   mesh.nodes(cms_opt.nodes.modal.number, :) = [0, 0, param.o, 0, 0, 0];
%!   mesh.nodes([cms_opt.nodes.interfaces.number], :) = [0, 0, -0.5 * param.l, 0, 0, 0;
%!                                                       0, 0,  0.5 * param.l, 0, 0, 0];

%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 1, cms_opt.nodes.interfaces(1).number);
%!   mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces(2).number);
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 3, cms_opt.nodes.modal.number);

%!   load_case.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));

%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(find([mesh.groups.tet10.id == 1])).elements) = 1;
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(param.E, param.nu);

%!   fprintf(stderr, "building cms element ...\n");
%!   [mesh, mat_ass, dof_map, sol_eig] = fem_cms_create(mesh, load_case, cms_opt);
%!   for i=1:numel(sol_eig.f)
%!     if (isfinite(sol_eig.f(i)))
%!       fprintf(stderr, "mode %d: %.1fHz\n", i, sol_eig.f(i) * (1 / SI_unit_second));
%!     endif
%!   endfor

%!   mesh_post_pro_file = sprintf("%s_post.msh", filename);

%!   fem_post_mesh_export(mesh_post_pro_file, mesh);

%!   for j=1:numel(sol_eig.f)
%!     eig_post_pro_file_mode{j} = sprintf("%s_eig_def_%03d.msh", filename, j);
%!     fem_post_sol_step_export(eig_post_pro_file_mode{j}, sol_eig, j, j, sol_eig.f(j), options.scale_def / max(norm(sol_eig.def(:, 1:3, j), "rows")));
%!   endfor

%!   eig_post_pro_file = sprintf("%s_modes_post.geo", filename);
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(eig_post_pro_file, "w");

%!     if (fd == -1)
%!       error("failed to open file \"%s\"", eig_post_pro_file);
%!     endif

%!     fprintf(fd, "Merge \"%s\";\n", mesh_post_pro_file);

%!     for j=1:numel(eig_post_pro_file_mode)
%!       fprintf(fd, "Merge \"%s\";\n", eig_post_pro_file_mode{j});
%!     endfor

%!     fputs(fd, "View.Type = 1;\n");
%!     fputs(fd, "View.VectorType = 5;\n");
%!     fputs(fd, "View.Visible = 1;\n");
%!     fputs(fd, "View.DisplacementFactor = 1;\n");
%!     fputs(fd, "View.ShowTime = 6;\n");
%!     fputs(fd, "View.ShowElement = 1;\n");
%!     fputs(fd, "View.IntervalsType = 3;\n");
%!     fputs(fd, "View.NbIso = 20;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect

%!   mat_ass.Dred = param.alpha * mat_ass.Mred + param.beta * mat_ass.Kred;
%!   fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);

%!   for i=1:numel(mbdyn_filenames)
%!     options_mbd(i).output_file = sprintf("%s_%d", filename, i);
%!     options_mbd(i).mbdyn_command = options.mbdyn_command;
%!     if (~options.verbose)
%!       options_mbd(i).logfile = sprintf("%s_%d.stdout", filename, i);
%!     endif
%!     options_mbd(i).f_run_mbdyn2easyanim = false;
%!     param_file = sprintf("%s_%d.set", filename, i);
%!     putenv("MBDYN_ROTOR_DYN_CMS_ELEM_FILE", [filename, "_cms.elm"]);
%!     putenv("MBDYN_ROTOR_DYN_CMS_PARAM_FILE", param_file);
%!     mbdyn_pre_write_param_file(param_file, param);
%!     mbdyn_solver_run(mbdyn_filenames{i}, options_mbd(i));
%!     res(i).log_dat = mbdyn_post_load_log(options_mbd(i).output_file);
%!     [res(i).t, res(i).trajectory, res(i).deformation, res(i).velocity, res(i).acceleration, res(i).node_id] = mbdyn_post_load_output_struct(options_mbd(i).output_file);
%!     res(i).log_dat.vars = mbdyn_post_id_to_index(res(i), res(i).log_dat.vars);
%!   endfor

%!   omega_crit = zeros(1, numel(mbdyn_filenames));

%!   for i=1:numel(mbdyn_filenames)
%!     dm = res(i).log_dat.vars.dm;
%!     dr = res(i).log_dat.vars.dr;
%!     m1 = res(i).log_dat.vars.D^2 * pi / 4 * res(i).log_dat.vars.w * res(i).log_dat.vars.rho;
%!     ecg = dm * dr / (m1 + dm);
%!     d = res(i).log_dat.vars.d;
%!     E = res(i).log_dat.vars.E;
%!     l = res(i).log_dat.vars.l;
%!     g = 9.81 / (SI_unit_meter / SI_unit_second^2);
%!     Iy = d^4 * pi / 64;
%!     f = m1 * g * l^3 / (48 * E * Iy);
%!     omega_crit(i) = sqrt(g / f);
%!   endfor

%!   figure("visible", "off");
%!   hold on;

%!   for i=1:numel(mbdyn_filenames)
%!     omega = res(i).velocity{res(i).log_dat.vars.node_idx_rotor}(:, 6);
%!     r = norm(res(i).deformation{res(i).log_dat.vars.node_idx_rotor}(:, 1:2), "rows");

%!     plot(omega * 30 / pi * (1 / SI_unit_second), 1e3 * r * SI_unit_meter, sprintf("-;%s;%d", printable_title(mbdyn_filename_suffix{i}), i));
%!   endfor

%!   xlabel("n [rpm]");
%!   ylabel("r [mm]");
%!   grid on;
%!   grid minor on;
%!   title("resonance curve center of disk versus speed - magnitude");
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 2
%! close all;
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!             0,  0.5 * b,  0.5 * c;  ##  2
%!             0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!             0,  0.5 * b, -0.5 * c;  ##  6
%!             0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!             a,  0.5 * b,  0.5 * c;  ##  9
%!             a, -0.5 * b,  0.5 * c;  ## 10
%!             a,  0.5 * b, -0.5 * c;  ## 11
%!             a, -0.5 * b, -0.5 * c,  ## 12
%!         a + d,        0,        0;  ## 13
%!            -e,        0,        0]; ## 14
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.element.name = "elem_id_cube";
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.modal.name = "node_id_cube_modal";
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.nodes.interfaces.name = "node_id_cube_interface1";
%! cms_opt.number_of_threads = 1;
%! cms_opt.algorithm = "eliminate";
%! cms_opt.invariants = true;
%! [mesh_cms, ...
%!  mat_ass_cms, ...
%!  dof_map_cms, ...
%!  sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 3
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E12 = 80000e6 / SI_unit_pascal;
%! param.E3 = 1000 * param.E12;
%! param.nu = 0.3;
%! param.rho = 1700 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = 15e-3 / SI_unit_meter;
%! param.d3 = 150e-3 / SI_unit_meter;
%! param.l1 = 120e-3 / SI_unit_meter;
%! param.l2 = 100e-3 / SI_unit_meter;
%! param.l3 = 15e-3 / SI_unit_meter;
%! param.h = 1.5e-3 / SI_unit_meter;
%! param.h3 = 10e-3 / SI_unit_meter;
%! OMEGA = 2 * pi * [0, 50,   0,   0;
%!                   0,  0, 300,   0;
%!                   0,  0,   0,   0] / (1 / SI_unit_second);
%! OMEGADOT = [0,  0,  0,  2e5;
%!             0,  0,  0,    0;
%!             0,  0,  0,    0];
%! idx = 1:columns(OMEGA);
%! OMEGA = OMEGA(:, idx);
%! OMEGADOT = OMEGADOT(:, idx);
%! opt_solver.pre_scaling = true;
%! opt_solver.refine_max_iter = int32(100);
%! opt_solver.solver = "pardiso";
%! opt_solver.number_of_threads = int32(4);
%! opt_solver.symmetric = false; ## FEM_MAT_STIFFNESS_OMEGA_DOT makes it unsymmetric
%! I1 = param.d1^4 * pi / 64;
%! I2 = param.d2^4 * pi / 64;

%! alpha = d11 = param.l1 * param.l2^2 / (3 * param.E12 * I1) + (param.l2^3 - param.l3^3) / (3 * param.E12 * I2);
%! delta = gamma = d12 = param.l1 * param.l2 / (3 * param.E12 * I1) + (param.l2^2 - param.l3^2) / (2 * param.E12 * I2);
%! beta = d22 = param.l1 / (3 * param.E12 * I1) + (param.l2 - param.l3) / (param.E12 * I2);
%! m = param.rho * pi / 4 * param.d3^2 * 2 * param.l3;
%! Ja = m * (3 * param.d3^2 + 4 * (2 * param.l3)^2) / 48;
%! Jp = m * param.d3^2 / 8;
%! lambda0_ref = sqrt((alpha * m + beta * Ja) / (2 * m * Ja * (alpha * beta - gamma^2)) * (1 + [-1, 1] * sqrt(1 - (4 * m * Ja * (alpha * beta - gamma^2))/(alpha * m + beta * Ja)^2)));
%! lambda_ref = zeros(4, columns(OMEGA));
%! for i=1:columns(OMEGA)
%!   lambda_ref(:, i) = roots([m * Ja * (alpha * beta - gamma^2);
%!                             -m * Jp * OMEGA(1, i) * (alpha * beta - gamma^2);
%!                             -(alpha * m + beta * Ja);
%!                             beta * Jp * OMEGA(1, i);
%!                             1]);
%! endfor
%! lambda_ref = real(lambda_ref);
%! options.number_of_modes = int32(10);
%! fref = zeros(5, columns(OMEGA));
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   modal = struct("mbdyn", cell(1, columns(OMEGA)));
%!   for i=1:columns(OMEGA)
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real E = %.16e;\n", param.E12);
%!       fprintf(fd, " set: real nu = %.16e;\n", param.nu);
%!       fprintf(fd, " set: real rho = %.16e;\n", param.rho);
%!       fprintf(fd, " set: real d1 = %.16e;\n", param.d1);
%!       fprintf(fd, " set: real d2 = %.16e;\n", param.d2);
%!       fprintf(fd, " set: real d3 = %.16e;\n", param.d3);
%!       fprintf(fd, " set: real l1 = %.16e;\n", param.l1);
%!       fprintf(fd, " set: real l2 = %.16e;\n", param.l2);
%!       fprintf(fd, " set: real l3 = %.16e;\n", param.l3);
%!       fputs(fd, " set: real G = E / (2. * (1 + nu));\n");
%!       fputs(fd, " set: real A1 = d1^2 * pi / 4.;\n");
%!       fputs(fd, " set: real A2 = d2^2 * pi / 4.;\n");
%!       fputs(fd, " set: real As1 = 9./10. * A1;\n");
%!       fputs(fd, " set: real As2 = 9./10. * A2;\n");
%!       fputs(fd, " set: real Iy1 = d1^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iy2 = d2^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iz1 = Iy1;\n");
%!       fputs(fd, " set: real Iz2 = Iy2;\n");
%!       fputs(fd, " set: real It1 = Iy1 + Iz1;\n");
%!       fputs(fd, " set: real It2 = Iy2 + Iz2;\n");
%!       fputs(fd, " set: real m = rho * pi / 4 * d3^2 * 2 * l3;\n");
%!       fputs(fd, " set: real Ja = m * (3 * d3^2 + 4 * (2 * l3)^2) / 48;\n");
%!       fputs(fd, " set: real Jp = m * d3^2 / 8;\n");
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGA%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGA(j, i));
%!       endfor
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGADOT%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGADOT(j, i));
%!       endfor
%!       fputs(fd, " set: real t1 = 1;\n");
%!       fputs(fd, " set: real N = 20000;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / N;\n");
%!       fputs(fd, "         threads: disable;\n");
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-6;\n");
%!       fputs(fd, "    linear solver: umfpack, scale, row max column max, always, max iterations, 3;\n");
%!       fputs(fd, "    method: implicit euler;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations;\n");
%!       fputs(fd, "    eigenanalysis: t1,\n");
%!       fputs(fd, "    output matrices, \n");
%!       fputs(fd, "    parameter, 1e-3,\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "         output geometry,\n");
%!       fputs(fd, "         lower frequency limit, 0.0001, upper frequency limit, 1000,\n");
%!       fputs(fd, "    use lapack,balance,permute,suffix format, \"%02d\";\n");
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    output meter: closest next, 0., forever, t1 / 100;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fputs(fd, "        rigid body kinematics: const, angular velocity, OMEGAx, OMEGAy, OMEGAz;\n");
%!       fputs(fd, "        rigid body kinematics: const, angular acceleration, OMEGADOTx, OMEGADOTy, OMEGADOTz;\n");
%!       fputs(fd, "    structural nodes: 9;\n");
%!       fputs(fd, "    rigid bodies: 1;\n");
%!       fputs(fd, "    beams: 4;\n");
%!       fputs(fd, "    joints: 2;\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "    structural: 1, static, \n");
%!       fputs(fd, "            reference, global, null, \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 2, static, \n");
%!       fputs(fd, "            reference, global, 0.25 * l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 3, static, \n");
%!       fputs(fd, "            reference, global, 0.5 * l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 4, static, \n");
%!       fputs(fd, "            reference, global, 0.75 * l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 5, static, \n");
%!       fputs(fd, "            reference, global, l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 6, static, \n");
%!       fputs(fd, "            reference, global, l1 + 0.25 * (l2 - l3), 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 7, static, \n");
%!       fputs(fd, "            reference, global, l1 + 0.5 * (l2 - l3), 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 8, static, \n");
%!       fputs(fd, "            reference, global, l1 + 0.75 * (l2 - l3), 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 9, dynamic, \n");
%!       fputs(fd, "            reference, global, l1 + l2, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null; \n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    body: 1, \n");
%!       fputs(fd, "            9,\n");
%!       fputs(fd, "            m, \n");
%!       fputs(fd, "            null, \n");
%!       fputs(fd, "            diag,   Jp, \n");
%!       fputs(fd, "                    Ja, \n");
%!       fputs(fd, "                    Ja,\n");
%!       fputs(fd, "            orientation, reference, global, eye;\n");
%!       fputs(fd, "        beam3: 1,\n");
%!       fputs(fd, "            1, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            2, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            3, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A1 , G * As1, G * As1, \n");
%!       fputs(fd, "                  G * It1, E * Iy1, E * Iz1,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "        beam3: 2,\n");
%!       fputs(fd, "            3, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            4, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            5, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A1 , G * As1, G * As1, \n");
%!       fputs(fd, "                  G * It1, E * Iy1, E * Iz1,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "        beam3: 3,\n");
%!       fputs(fd, "            5, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            6, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            7, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A2 , G * As2, G * As2, \n");
%!       fputs(fd, "                  G * It2, E * Iy2, E * Iz2,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "       beam3: 4,\n");
%!       fputs(fd, "            7, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            8, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            9, position, reference, global, l1 + l2 - l3, 0., 0.,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A2 , G * As2, G * As2, \n");
%!       fputs(fd, "                  G * It2, E * Iy2, E * Iz2,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "    joint: 1, total pin joint,\n");
%!       fputs(fd, "                    1,\n");
%!       fputs(fd, "                            position,		reference, global, null,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            position,		reference, global, null,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, "    joint: 2, total pin joint,\n");
%!       fputs(fd, "                    5,\n");
%!       fputs(fd, "                            position,		reference, global, l1, 0., 0.,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            position,		reference, global, l1, 0., 0.,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    inactive, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     modal(i).mbdyn = mbdyn_post_load_output_eig(opt_mbdyn.output_file);
%!     for j=1:numel(modal(i).mbdyn.f)
%!       opt_modal.mode_index = j;
%!       opt_modal.scale = 100;
%!       mode_file = [opt_mbdyn.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(opt_mbdyn.output_file, [mode_file, ".mov"], opt_modal, modal(i).mbdyn);
%!       [err, msg] = symlink([opt_mbdyn.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = true;
%!       opt_post.f_runEasyAnim = false;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!     fref(:, i) = modal(i).mbdyn.f(:);
%!   endfor
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.Tolerance = 1e-6;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "s2[] = Extrude {l2 - l3, 0, 0}{ Line{s1[0]}; Layers{Ceil((l2 - l3) / h)}; Recombine; };\n");
%!     fputs(fd, "s3[] = Extrude {2 * l3, 0, 0}{ Line{s2[0]}; Layers{Ceil(2 * l3 / h3)}; Recombine; };\n");
%!     fputs(fd, "s4[] = Extrude {0, 0.5 * (d2 - d1), 0}{ Line{s2[3],s3[3]}; Layers{Ceil(0.5 * (d2 - d1) / h)}; Recombine; };\n");
%!     fputs(fd, "s6[] = Extrude {0, 0.5 * (d3 - d2), 0}{ Line{s4[4]}; Layers{Ceil(0.5 * (d3 - d1) / h3)}; Recombine; };\n");
%!     fputs(fd, "se0[] = {s1[1], s2[1], s3[1], s4[1], s4[5], s6[1]};\n");
%!     fputs(fd, "v1[] = {};\n");
%!     fputs(fd, "se1[] = {};\n");
%!     fputs(fd, "For i In {0:#se0[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se0[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v1[i] = vtmp[1];\n");
%!     fputs(fd, "  se1[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v2[] = {};\n");
%!     fputs(fd, "se2[] = {};\n");
%!     fputs(fd, "For i In {0:#se1[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se1[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v2[i] = vtmp[1];\n");
%!     fputs(fd, "  se2[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v3[] = {};\n");
%!     fputs(fd, "se3[] = {};\n");
%!     fputs(fd, "For i In {0:#se2[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se2[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v3[i] = vtmp[1];\n");
%!     fputs(fd, "  se3[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v4[] = {};\n");
%!     fputs(fd, "se4[] = {};\n");
%!     fputs(fd, "For i In {0:#se3[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se3[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v4[i] = vtmp[1];\n");
%!     fputs(fd, "  se4[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v[] = {v1[], v2[], v3[], v4[]};\n");
%!     fputs(fd, "Physical Volume(1) = {19, 1, 13, 7};\n");
%!     fputs(fd, "Physical Volume(2) = {22, 4, 16, 10, 14, 8, 20, 2};\n");
%!     fputs(fd, "Physical Volume(3) = {18, 12, 24, 6, 23, 5, 17, 11, 21, 3, 15, 9};\n");
%!     fputs(fd, "Physical Surface(1) = {89, 8, 62, 35};\n");
%!     fputs(fd, "Physical Surface(2) = {19, 100, 73, 46};\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   e1 = [1; 0.5; 0.4];
%!   e2 = [0; 1; 0];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   R *= diag(1 ./ norm(R, "cols"));
%!   T = [R.', zeros(3, 3);
%!        zeros(3, 3), R.'];
%!   mesh.nodes *= T;
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, "iso4");
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, "iso4");
%!   mesh.elements.joints(2).nodes = node_idx_bearing2;
%!   mesh.elements.joints(2).C = [0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0] * T;
%!   mesh.elements.joints(1).nodes = node_idx_bearing1;
%!   mesh.elements.joints(1).C = [1, 0, 0, 0, 0, 0;
%!                                0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0;
%!                                0, 0, 0, 1, 0, 0] * T;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   for i=1:3
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == i])).elements]) = i;
%!   endfor
%!   mesh.material_data(1).rho = 0;
%!   mesh.material_data(1).E = param.E12;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(2).rho = 0;
%!   mesh.material_data(2).E = param.E12;
%!   mesh.material_data(2).nu = param.nu;
%!   mesh.material_data(3).rho = param.rho;
%!   mesh.material_data(3).E = param.E3;
%!   mesh.material_data(3).nu = param.nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = opt_solver.number_of_threads;
%!   load_case = fem_pre_load_case_create_empty(columns(OMEGA));
%!   empty_cell = cell(1, columns(OMEGA));
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   sol_eig = struct("lambda", empty_cell, "f", empty_cell, "def", empty_cell);
%!   mat_red = struct("Mred", empty_cell, "Kred", empty_cell, "Dred", empty_cell);
%!   for i=1:columns(OMEGA)
%!     fprintf(stderr, "n=%.0frpm\n", norm(OMEGA(:, i)) * 30 / pi / SI_unit_second);
%!     load_case(i).omega = R * OMEGA(:, i);
%!     load_case(i).omegadot = R * OMEGADOT(:, i);
%!     [mat_ass.M, ...
%!      mat_ass.K, ...
%!      mat_ass.KOMEGA, ...
%!      mat_ass.KOMEGA_DOT, ...
%!      mat_ass.DOMEGA, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_MASS, ...
%!                                   FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                   FEM_MAT_DAMPING_OMEGA, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS, ...
%!                                   FEM_VEC_INERTIA_M1, ...
%!                                   FEM_MAT_INERTIA_J], ...
%!                                  load_case(i));
%!     Tred = zeros(rows(mesh.nodes) * 3, 6);
%!     for j=1:rows(mesh.nodes)
%!       Tred((j - 1) * 3 + (1:3), 1:3) = eye(3);
%!       Tred((j - 1) * 3 + (1:3), 4:6) = -skew(mesh.nodes(j, 1:3));
%!     endfor
%!     idx = dof_map.ndof(:, 1:3).'(:);
%!     mat_red(i).Mred = Tred.' * mat_ass.M(idx, idx) * Tred;
%!     mat_red(i).Kred = Tred.' * mat_ass.K(idx, idx) * Tred;
%!     mat_red(i).Dred = Tred.' * mat_ass.DOMEGA(idx, idx) * Tred;
%!     sol_stat(i).def = fem_sol_static(mesh, dof_map, mat_ass).def;
%!     sol_stat(i).stress = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         load_case(i), ...
%!                                         sol_stat(i));
%!     load_case_pre_stress.tau0 = sol_stat(i).stress.tau;
%!     mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress);
%!     mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0 + mat_ass.KOMEGA_DOT;
%!     mat_ass.D = mat_ass.DOMEGA;
%!     sol_eig(i) = fem_sol_modal_damped(mesh, ...
%!                                       dof_map, ...
%!                                       mat_ass, ...
%!                                       options.number_of_modes, ...
%!                                       opt_solver);
%!   endfor
%!   f = zeros(options.number_of_modes, numel(sol_eig));
%!   for i=1:columns(f)
%!     f(:, i) = sort(sol_eig(i).f(:));
%!   endfor
%!   tol = 2e-2;
%!   assert(f(floor(end/2+1):end,:), fref, tol * max(max(abs(fref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 4
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E12 = 80000e6 / SI_unit_pascal;
%! param.E3 =  210000e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.rho = 1700 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = 15e-3 / SI_unit_meter;
%! param.d3 = 150e-3 / SI_unit_meter;
%! param.l1 = 120e-3 / SI_unit_meter;
%! param.l2 = 100e-3 / SI_unit_meter;
%! param.l3 = 15e-3 / SI_unit_meter;
%! param.h = 10e-3 / SI_unit_meter;
%! param.h3 = 10e-3 / SI_unit_meter;
%! OMEGA = 2 * pi * [0, 50,   0,   0;
%!                   0,  0, 150,   0;
%!                   0,  0,   0,   0] / (1 / SI_unit_second);
%! OMEGADOT = [0,  0,  0,  2e5;
%!             0,  0,  0,    0;
%!             0,  0,  0,    0] / (1 / SI_unit_second^2);
%! idx = 1:columns(OMEGA);
%! OMEGA = OMEGA(:, idx);
%! OMEGADOT = OMEGADOT(:, idx);
%! options.post_proc_modes = false;
%! options.verbose = false;
%! opt_solver.pre_scaling = true;
%! opt_solver.refine_max_iter = int32(100);
%! opt_solver.solver = "pardiso";
%! opt_solver.number_of_threads = int32(4);
%! opt_solver.symmetric = false; ## FEM_MAT_STIFFNESS_OMEGA_DOT makes it unsymmetric
%! I1 = param.d1^4 * pi / 64;
%! I2 = param.d2^4 * pi / 64;
%! alpha = d11 = param.l1 * param.l2^2 / (3 * param.E12 * I1) + (param.l2^3 - param.l3^3) / (3 * param.E12 * I2);
%! delta = gamma = d12 = param.l1 * param.l2 / (3 * param.E12 * I1) + (param.l2^2 - param.l3^2) / (2 * param.E12 * I2);
%! beta = d22 = param.l1 / (3 * param.E12 * I1) + (param.l2 - param.l3) / (param.E12 * I2);
%! m = param.rho * pi / 4 * param.d3^2 * 2 * param.l3;
%! Ja = m * (3 * param.d3^2 + 4 * (2 * param.l3)^2) / 48;
%! Jp = m * param.d3^2 / 8;
%! lambda0_ref = sqrt((alpha * m + beta * Ja) / (2 * m * Ja * (alpha * beta - gamma^2)) * (1 + [-1, 1] * sqrt(1 - (4 * m * Ja * (alpha * beta - gamma^2))/(alpha * m + beta * Ja)^2)));
%! lambda_ref = zeros(4, columns(OMEGA));
%! for i=1:columns(OMEGA)
%!   lambda_ref(:, i) = roots([m * Ja * (alpha * beta - gamma^2);
%!                             -m * Jp * OMEGA(1, i) * (alpha * beta - gamma^2);
%!                             -(alpha * m + beta * Ja);
%!                             beta * Jp * OMEGA(1, i);
%!                             1]);
%! endfor
%! lambda_ref = real(lambda_ref);
%! options.number_of_modes = int32(10);
%! fref = zeros(5, columns(OMEGA));
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   modal = struct("mbdyn", cell(1, columns(OMEGA)));
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.Tolerance = 1e-6;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "s2[] = Extrude {l2 - l3, 0, 0}{ Line{s1[0]}; Layers{Ceil((l2 - l3) / h)}; Recombine; };\n");
%!     fputs(fd, "s3[] = Extrude {2 * l3, 0, 0}{ Line{s2[0]}; Layers{Ceil(2 * l3 / h3)}; Recombine; };\n");
%!     fputs(fd, "s4[] = Extrude {0, 0.5 * (d2 - d1), 0}{ Line{s2[3],s3[3]}; Layers{Ceil(0.5 * (d2 - d1) / h)}; Recombine; };\n");
%!     fputs(fd, "s6[] = Extrude {0, 0.5 * (d3 - d2), 0}{ Line{s4[4]}; Layers{Ceil(0.5 * (d3 - d1) / h3)}; Recombine; };\n");
%!     fputs(fd, "se0[] = {s1[1], s2[1], s3[1], s4[1], s4[5], s6[1]};\n");
%!     fputs(fd, "v1[] = {};\n");
%!     fputs(fd, "se1[] = {};\n");
%!     fputs(fd, "For i In {0:#se0[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se0[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v1[i] = vtmp[1];\n");
%!     fputs(fd, "  se1[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v2[] = {};\n");
%!     fputs(fd, "se2[] = {};\n");
%!     fputs(fd, "For i In {0:#se1[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se1[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v2[i] = vtmp[1];\n");
%!     fputs(fd, "  se2[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v3[] = {};\n");
%!     fputs(fd, "se3[] = {};\n");
%!     fputs(fd, "For i In {0:#se2[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se2[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v3[i] = vtmp[1];\n");
%!     fputs(fd, "  se3[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v4[] = {};\n");
%!     fputs(fd, "se4[] = {};\n");
%!     fputs(fd, "For i In {0:#se3[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se3[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v4[i] = vtmp[1];\n");
%!     fputs(fd, "  se4[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v[] = {v1[], v2[], v3[], v4[]};\n");
%!     fputs(fd, "Physical Volume(1) = {19, 1, 13, 7};\n");
%!     fputs(fd, "Physical Volume(2) = {22, 4, 16, 10, 14, 8, 20, 2};\n");
%!     fputs(fd, "Physical Volume(3) = {18, 12, 24, 6, 23, 5, 17, 11, 21, 3, 15, 9};\n");
%!     fputs(fd, "Physical Surface(1) = {89, 8, 62, 35};\n");
%!     fputs(fd, "Physical Surface(2) = {19, 100, 73, 46};\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   e1 = [1; 0.5; 0.4];
%!   e2 = [0; 1; 0];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   R *= diag(1 ./ norm(R, "cols"));
%!   R = eye(3);
%!   T = [R.', zeros(3, 3);
%!        zeros(3, 3), R.'];
%!   mesh.nodes *= T;
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, "iso4");
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, "iso4");
%!   mesh.elements.joints(2).nodes = node_idx_bearing2;
%!   mesh.elements.joints(2).C = [0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0] * T;
%!   mesh.elements.joints(1).nodes = node_idx_bearing1;
%!   mesh.elements.joints(1).C = [1, 0, 0, 0, 0, 0;
%!                                0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0;
%!                                0, 0, 0, 1, 0, 0] * T;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   for i=1:3
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == i])).elements]) = i;
%!   endfor
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E12;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(2).rho = param.rho;
%!   mesh.material_data(2).E = param.E12;
%!   mesh.material_data(2).nu = param.nu;
%!   mesh.material_data(3).rho = param.rho;
%!   mesh.material_data(3).E = param.E3;
%!   mesh.material_data(3).nu = param.nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = opt_solver.number_of_threads;
%!   load_case = fem_pre_load_case_create_empty(columns(OMEGA));
%!   empty_cell = cell(1, columns(OMEGA));
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   sol_eig = struct("lambda", empty_cell, "f", empty_cell, "def", empty_cell);
%!   mat_red = struct("Mred", empty_cell, "Kred", empty_cell, "Dred", empty_cell);
%!   for i=1:columns(OMEGA)
%!     fprintf(stderr, "n=%.0frpm\n", norm(OMEGA(:, i)) * 30 / pi / SI_unit_second);
%!     load_case(i).omega = R * OMEGA(:, i);
%!     load_case(i).omegadot = R * OMEGADOT(:, i);
%!     [mat_ass.M, ...
%!      mat_ass.K, ...
%!      mat_ass.KOMEGA, ...
%!      mat_ass.KOMEGA_DOT, ...
%!      mat_ass.DOMEGA, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_MASS, ...
%!                                   FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                   FEM_MAT_DAMPING_OMEGA, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS, ...
%!                                   FEM_VEC_INERTIA_M1, ...
%!                                   FEM_MAT_INERTIA_J], ...
%!                                  load_case(i));
%!     Tred = zeros(rows(mesh.nodes) * 3, 6);
%!     for j=1:rows(mesh.nodes)
%!       Tred((j - 1) * 3 + (1:3), 1:3) = eye(3);
%!       Tred((j - 1) * 3 + (1:3), 4:6) = -skew(mesh.nodes(j, 1:3));
%!     endfor
%!     idx = dof_map.ndof(:, 1:3).'(:);
%!     mat_red(i).Mred = Tred.' * mat_ass.M(idx, idx) * Tred;
%!     mat_red(i).Kred = Tred.' * mat_ass.K(idx, idx) * Tred;
%!     mat_red(i).Dred = Tred.' * mat_ass.DOMEGA(idx, idx) * Tred;
%!     sol_stat(i).def = fem_sol_static(mesh, dof_map, mat_ass).def;
%!     sol_stat(i).stress = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         load_case(i), ...
%!                                         sol_stat(i));
%!     load_case_pre_stress.tau0 = sol_stat(i).stress.tau;
%!     mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress);
%!     mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0 + mat_ass.KOMEGA_DOT;
%!     mat_ass.D = mat_ass.DOMEGA;
%!     sol_eig(i) = fem_sol_modal_damped(mesh, ...
%!                                       dof_map, ...
%!                                       mat_ass, ...
%!                                       options.number_of_modes, ...
%!                                       opt_solver);
%!   endfor
%!   for i=1:columns(OMEGA)
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.joints.number = 2;
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     load_case_empty = struct();
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real E = %.16e;\n", param.E12);
%!       fprintf(fd, " set: real nu = %.16e;\n", param.nu);
%!       fprintf(fd, " set: real rho = %.16e;\n", param.rho);
%!       fprintf(fd, " set: real d1 = %.16e;\n", param.d1);
%!       fprintf(fd, " set: real d2 = %.16e;\n", param.d2);
%!       fprintf(fd, " set: real d3 = %.16e;\n", param.d3);
%!       fprintf(fd, " set: real l1 = %.16e;\n", param.l1);
%!       fprintf(fd, " set: real l2 = %.16e;\n", param.l2);
%!       fprintf(fd, " set: real l3 = %.16e;\n", param.l3);
%!       fputs(fd, " set: real G = E / (2. * (1 + nu));\n");
%!       fputs(fd, " set: real A1 = d1^2 * pi / 4.;\n");
%!       fputs(fd, " set: real A2 = d2^2 * pi / 4.;\n");
%!       fputs(fd, " set: real As1 = 9./10. * A1;\n");
%!       fputs(fd, " set: real As2 = 9./10. * A2;\n");
%!       fputs(fd, " set: real Iy1 = d1^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iy2 = d2^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iz1 = Iy1;\n");
%!       fputs(fd, " set: real Iz2 = Iy2;\n");
%!       fputs(fd, " set: real It1 = Iy1 + Iz1;\n");
%!       fputs(fd, " set: real It2 = Iy2 + Iz2;\n");
%!       fputs(fd, " set: real m = rho * pi / 4 * d3^2 * 2 * l3;\n");
%!       fputs(fd, " set: real Ja = m * (3 * d3^2 + 4 * (2 * l3)^2) / 48;\n");
%!       fputs(fd, " set: real Jp = m * d3^2 / 8;\n");
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGA%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGA(j, i));
%!       endfor
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGADOT%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGADOT(j, i));
%!       endfor
%!       fprintf(fd, " set: real t1 = %g;\n", 1000 / SI_unit_second);
%!       fputs(fd, " set: real N = 1000;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / N;\n");
%!       fputs(fd, "         threads: assembly, 4;\n");
%!       fputs(fd, "         threads: solver, 4;\n");
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-3, 1e-3;\n");
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!       fputs(fd, "    method: bdf;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    eigenanalysis: t1,\n");
%!       fputs(fd, "    # output matrices, \n");
%!       fputs(fd, "    # parameter, 1e-3, ## use default estimate\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "         # output geometry,\n");
%!       fprintf(fd, "         lower frequency limit, %e,\n", 0.01 / (SI_unit_second^-1));
%!       fprintf(fd, "         upper frequency limit, %e,\n", 1000 / (SI_unit_second^-1));
%!       fprintf(fd, "    use arpack,%d,%d,0.,suffix format, \"%%02d\";\n", 2 * options.number_of_modes, options.number_of_modes * 10);
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-8, forcing term max tolerance, 1e-3;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "        angular acceleration,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGADOTx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    joint: 1, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, "    joint: 2, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    inactive, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     if (~options.verbose)
%!       opt_mbdyn.logfile = [filename, ".stdout"];
%!     endif
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     [mesh_sol(i), sol(i)] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!     modal(i).mbdyn = mbdyn_post_load_output_eig(opt_mbdyn.output_file);
%!     if (options.post_proc_modes)
%!     for j=1:numel(modal(i).mbdyn.f)
%!       opt_modal.mode_index = j;
%!       opt_modal.scale = 100;
%!       mode_file = [opt_mbdyn.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(opt_mbdyn.output_file, [mode_file, ".mov"], opt_modal, modal(i).mbdyn);
%!       [err, msg] = symlink([opt_mbdyn.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = true;
%!       opt_post.f_runEasyAnim = false;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!     endif
%!     fref(:, i) = modal(i).mbdyn.f(1:rows(fref));
%!   endfor
%!   f = zeros(options.number_of_modes, numel(sol_eig));
%!   for i=1:columns(f)
%!     f(:, i) = sort(sol_eig(i).f(:));
%!   endfor
%!   tol = 2e-2;
%!   assert(f(floor(end/2+1):end,:), fref, tol * max(max(abs(fref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 5
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E12 = 80000e6 / SI_unit_pascal;
%! param.E3 =  210000e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.rho = 1700 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = 15e-3 / SI_unit_meter;
%! param.d3 = 150e-3 / SI_unit_meter;
%! param.l1 = 120e-3 / SI_unit_meter;
%! param.l2 = 100e-3 / SI_unit_meter;
%! param.l3 = 15e-3 / SI_unit_meter;
%! param.h = 40e-3 / SI_unit_meter;
%! param.h3 = 10e-3 / SI_unit_meter;
%! OMEGA = 2 * pi * [0, 50,   0,   0;
%!                   0,  0, 150,   0;
%!                   0,  0,   0,   0] / (1 / SI_unit_second);
%! OMEGADOT = [0,  0,  0,  2e5;
%!             0,  0,  0,    0;
%!             0,  0,  0,    0] / (1 / SI_unit_second^2);
%! idx = 1:columns(OMEGA);
%! OMEGA = OMEGA(:, idx);
%! OMEGADOT = OMEGADOT(:, idx);
%! options.post_proc_modes = false;
%! options.verbose = false;
%! opt_solver.pre_scaling = true;
%! opt_solver.refine_max_iter = int32(100);
%! opt_solver.solver = "pardiso";
%! opt_solver.number_of_threads = int32(4);
%! opt_solver.symmetric = false; ## FEM_MAT_STIFFNESS_OMEGA_DOT makes it unsymmetric
%! I1 = param.d1^4 * pi / 64;
%! I2 = param.d2^4 * pi / 64;
%! alpha = d11 = param.l1 * param.l2^2 / (3 * param.E12 * I1) + (param.l2^3 - param.l3^3) / (3 * param.E12 * I2);
%! delta = gamma = d12 = param.l1 * param.l2 / (3 * param.E12 * I1) + (param.l2^2 - param.l3^2) / (2 * param.E12 * I2);
%! beta = d22 = param.l1 / (3 * param.E12 * I1) + (param.l2 - param.l3) / (param.E12 * I2);
%! m = param.rho * pi / 4 * param.d3^2 * 2 * param.l3;
%! Ja = m * (3 * param.d3^2 + 4 * (2 * param.l3)^2) / 48;
%! Jp = m * param.d3^2 / 8;
%! lambda0_ref = sqrt((alpha * m + beta * Ja) / (2 * m * Ja * (alpha * beta - gamma^2)) * (1 + [-1, 1] * sqrt(1 - (4 * m * Ja * (alpha * beta - gamma^2))/(alpha * m + beta * Ja)^2)));
%! lambda_ref = zeros(4, columns(OMEGA));
%! for i=1:columns(OMEGA)
%!   lambda_ref(:, i) = roots([m * Ja * (alpha * beta - gamma^2);
%!                             -m * Jp * OMEGA(1, i) * (alpha * beta - gamma^2);
%!                             -(alpha * m + beta * Ja);
%!                             beta * Jp * OMEGA(1, i);
%!                             1]);
%! endfor
%! lambda_ref = real(lambda_ref);
%! options.number_of_modes = int32(10);
%! fref = zeros(5, columns(OMEGA));
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   modal = struct("mbdyn", cell(1, columns(OMEGA)));
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.Tolerance = 1e-6;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "s2[] = Extrude {l2 - l3, 0, 0}{ Line{s1[0]}; Layers{Ceil((l2 - l3) / h)}; Recombine; };\n");
%!     fputs(fd, "s3[] = Extrude {2 * l3, 0, 0}{ Line{s2[0]}; Layers{Ceil(2 * l3 / h3)}; Recombine; };\n");
%!     fputs(fd, "s4[] = Extrude {0, 0.5 * (d2 - d1), 0}{ Line{s2[3],s3[3]}; Layers{Ceil(0.5 * (d2 - d1) / h)}; Recombine; };\n");
%!     fputs(fd, "s6[] = Extrude {0, 0.5 * (d3 - d2), 0}{ Line{s4[4]}; Layers{Ceil(0.5 * (d3 - d1) / h3)}; Recombine; };\n");
%!     fputs(fd, "se0[] = {s1[1], s2[1], s3[1], s4[1], s4[5], s6[1]};\n");
%!     fputs(fd, "v1[] = {};\n");
%!     fputs(fd, "se1[] = {};\n");
%!     fputs(fd, "For i In {0:#se0[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se0[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v1[i] = vtmp[1];\n");
%!     fputs(fd, "  se1[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v2[] = {};\n");
%!     fputs(fd, "se2[] = {};\n");
%!     fputs(fd, "For i In {0:#se1[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se1[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v2[i] = vtmp[1];\n");
%!     fputs(fd, "  se2[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v3[] = {};\n");
%!     fputs(fd, "se3[] = {};\n");
%!     fputs(fd, "For i In {0:#se2[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se2[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v3[i] = vtmp[1];\n");
%!     fputs(fd, "  se3[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v4[] = {};\n");
%!     fputs(fd, "se4[] = {};\n");
%!     fputs(fd, "For i In {0:#se3[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se3[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v4[i] = vtmp[1];\n");
%!     fputs(fd, "  se4[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v[] = {v1[], v2[], v3[], v4[]};\n");
%!     fputs(fd, "Physical Volume(1) = {19, 1, 13, 7};\n");
%!     fputs(fd, "Physical Volume(2) = {22, 4, 16, 10, 14, 8, 20, 2};\n");
%!     fputs(fd, "Physical Volume(3) = {18, 12, 24, 6, 23, 5, 17, 11, 21, 3, 15, 9};\n");
%!     fputs(fd, "Physical Surface(1) = {89, 8, 62, 35};\n");
%!     fputs(fd, "Physical Surface(2) = {19, 100, 73, 46};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   e1 = [1; 0.5; 0.4];
%!   e2 = [0; 1; 0];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   R *= diag(1 ./ norm(R, "cols"));
%!   R = eye(3);
%!   T = [R.', zeros(3, 3);
%!        zeros(3, 3), R.'];
%!   mesh.nodes *= T;
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, "quad8");
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, "tria6h");
%!   mesh.elements.joints(2).nodes = node_idx_bearing2;
%!   mesh.elements.joints(2).C = [0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0] * T;
%!   mesh.elements.joints(1).nodes = node_idx_bearing1;
%!   mesh.elements.joints(1).C = [1, 0, 0, 0, 0, 0;
%!                                0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0;
%!                                0, 0, 0, 1, 0, 0] * T;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   for i=1:3
%!     mesh.materials.iso20([mesh.groups.iso20(find([[mesh.groups.iso20.id] == i])).elements]) = i;
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == i])).elements]) = i;
%!   endfor
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E12;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(2).rho = param.rho;
%!   mesh.material_data(2).E = param.E12;
%!   mesh.material_data(2).nu = param.nu;
%!   mesh.material_data(3).rho = param.rho;
%!   mesh.material_data(3).E = param.E3;
%!   mesh.material_data(3).nu = param.nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = opt_solver.number_of_threads;
%!   load_case = fem_pre_load_case_create_empty(columns(OMEGA));
%!   empty_cell = cell(1, columns(OMEGA));
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   sol_eig = struct("lambda", empty_cell, "f", empty_cell, "def", empty_cell);
%!   mat_red = struct("Mred", empty_cell, "Kred", empty_cell, "Dred", empty_cell);
%!   for i=1:columns(OMEGA)
%!     fprintf(stderr, "n=%.0frpm\n", norm(OMEGA(:, i)) * 30 / pi / SI_unit_second);
%!     load_case(i).omega = R * OMEGA(:, i);
%!     load_case(i).omegadot = R * OMEGADOT(:, i);
%!     [mat_ass.M, ...
%!      mat_ass.K, ...
%!      mat_ass.KOMEGA, ...
%!      mat_ass.KOMEGA_DOT, ...
%!      mat_ass.DOMEGA, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_MASS, ...
%!                                   FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                   FEM_MAT_DAMPING_OMEGA, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS, ...
%!                                   FEM_VEC_INERTIA_M1, ...
%!                                   FEM_MAT_INERTIA_J], ...
%!                                  load_case(i));
%!     Tred = zeros(rows(mesh.nodes) * 3, 6);
%!     for j=1:rows(mesh.nodes)
%!       Tred((j - 1) * 3 + (1:3), 1:3) = eye(3);
%!       Tred((j - 1) * 3 + (1:3), 4:6) = -skew(mesh.nodes(j, 1:3));
%!     endfor
%!     idx = dof_map.ndof(:, 1:3).'(:);
%!     mat_red(i).Mred = Tred.' * mat_ass.M(idx, idx) * Tred;
%!     mat_red(i).Kred = Tred.' * mat_ass.K(idx, idx) * Tred;
%!     mat_red(i).Dred = Tred.' * mat_ass.DOMEGA(idx, idx) * Tred;
%!     sol_stat(i).def = fem_sol_static(mesh, dof_map, mat_ass).def;
%!     sol_stat(i).stress = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         load_case(i), ...
%!                                         sol_stat(i));
%!     load_case_pre_stress.tau0 = sol_stat(i).stress.tau;
%!     mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress);
%!     mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0 + mat_ass.KOMEGA_DOT;
%!     mat_ass.D = mat_ass.DOMEGA;
%!     sol_eig(i) = fem_sol_modal_damped(mesh, ...
%!                                       dof_map, ...
%!                                       mat_ass, ...
%!                                       options.number_of_modes, ...
%!                                       opt_solver);
%!   endfor
%!   for i=1:columns(OMEGA)
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.joints.number = 2;
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     load_case_empty = struct();
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real E = %.16e;\n", param.E12);
%!       fprintf(fd, " set: real nu = %.16e;\n", param.nu);
%!       fprintf(fd, " set: real rho = %.16e;\n", param.rho);
%!       fprintf(fd, " set: real d1 = %.16e;\n", param.d1);
%!       fprintf(fd, " set: real d2 = %.16e;\n", param.d2);
%!       fprintf(fd, " set: real d3 = %.16e;\n", param.d3);
%!       fprintf(fd, " set: real l1 = %.16e;\n", param.l1);
%!       fprintf(fd, " set: real l2 = %.16e;\n", param.l2);
%!       fprintf(fd, " set: real l3 = %.16e;\n", param.l3);
%!       fputs(fd, " set: real G = E / (2. * (1 + nu));\n");
%!       fputs(fd, " set: real A1 = d1^2 * pi / 4.;\n");
%!       fputs(fd, " set: real A2 = d2^2 * pi / 4.;\n");
%!       fputs(fd, " set: real As1 = 9./10. * A1;\n");
%!       fputs(fd, " set: real As2 = 9./10. * A2;\n");
%!       fputs(fd, " set: real Iy1 = d1^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iy2 = d2^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iz1 = Iy1;\n");
%!       fputs(fd, " set: real Iz2 = Iy2;\n");
%!       fputs(fd, " set: real It1 = Iy1 + Iz1;\n");
%!       fputs(fd, " set: real It2 = Iy2 + Iz2;\n");
%!       fputs(fd, " set: real m = rho * pi / 4 * d3^2 * 2 * l3;\n");
%!       fputs(fd, " set: real Ja = m * (3 * d3^2 + 4 * (2 * l3)^2) / 48;\n");
%!       fputs(fd, " set: real Jp = m * d3^2 / 8;\n");
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGA%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGA(j, i));
%!       endfor
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGADOT%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGADOT(j, i));
%!       endfor
%!       fprintf(fd, " set: real t1 = %g;\n", 1000 / SI_unit_second);
%!       fputs(fd, " set: real N = 1000;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / N;\n");
%!       fputs(fd, "         threads: assembly, 4;\n");
%!       fputs(fd, "         threads: solver, 4;\n");
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-3, 1e-3;\n");
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!       fputs(fd, "    method: bdf;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    eigenanalysis: t1,\n");
%!       fputs(fd, "    # output matrices, \n");
%!       fputs(fd, "    # parameter, 1e-3, ## use default estimate\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "         # output geometry,\n");
%!       fprintf(fd, "         lower frequency limit, %e,\n", 0.01 / (SI_unit_second^-1));
%!       fprintf(fd, "         upper frequency limit, %e,\n", 1000 / (SI_unit_second^-1));
%!       fprintf(fd, "    use arpack,%d,%d,0.,suffix format, \"%%02d\";\n", 2 * options.number_of_modes, options.number_of_modes * 10);
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-8, forcing term max tolerance, 1e-3, inner iterations before assembly, 30;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "        angular acceleration,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGADOTx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    joint: 1, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, "    joint: 2, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    inactive, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     if (~options.verbose)
%!       opt_mbdyn.logfile = [filename, ".stdout"];
%!     endif
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     [mesh_sol(i), sol(i)] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!     modal(i).mbdyn = mbdyn_post_load_output_eig(opt_mbdyn.output_file);
%!     if (options.post_proc_modes)
%!     for j=1:numel(modal(i).mbdyn.f)
%!       opt_modal.mode_index = j;
%!       opt_modal.scale = 100;
%!       mode_file = [opt_mbdyn.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(opt_mbdyn.output_file, [mode_file, ".mov"], opt_modal, modal(i).mbdyn);
%!       [err, msg] = symlink([opt_mbdyn.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = true;
%!       opt_post.f_runEasyAnim = false;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!     endif
%!     fref(:, i) = modal(i).mbdyn.f(1:rows(fref));
%!   endfor
%!   f = zeros(options.number_of_modes, numel(sol_eig));
%!   for i=1:columns(f)
%!     f(:, i) = sort(sol_eig(i).f(:));
%!   endfor
%!   tol = 2e-2;
%!   assert(f(floor(end/2+1):end,:), fref, tol * max(max(abs(fref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 6
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20"))
%!     mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!     mesh.materials.iso20([mesh.groups.iso20(find([[mesh.groups.iso20.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fputs(fd, "         threads: assembly, 4;\n");
%!     fputs(fd, "         threads: solver, 4;\n");
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert(N, Nref, tol * norm(Nref));
%!   assert(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 7
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 0.8333e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso4", "iso8", "tria3", "penta6"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"iso4"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"iso4"});
%!   if (isfield(mesh.elements, "iso8"))
%!     mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fputs(fd, "         threads: assembly, 4;\n");
%!     fputs(fd, "         threads: solver, 4;\n");
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 1.5e-2;
%!   assert(N, Nref, tol * norm(Nref));
%!   assert(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 8
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 1.25e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fprintf(fd, "h = %e;\n", param.h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "Point(2) = {l1, 0, 0};\n");
%!     fputs(fd, "Point(3) = {l1, 0.5 * d1, 0};\n");
%!     fputs(fd, "Point(4) = {0, 0.5 * d1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "v[] = Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi}{ Surface{1}; };\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v[1]}; } } = h;\n");
%!     fputs(fd, "Physical Volume(1) = {v[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4};\n");
%!     fputs(fd, "Physical Surface(2) = {2};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"tria6h"});
%!   if (isfield(mesh.elements, "tet10h"))
%!     mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!     mesh.materials.tet10h([mesh.groups.tet10h(find([[mesh.groups.tet10h.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 50;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fputs(fd, "         threads: assembly, 4;\n");
%!     fputs(fd, "         threads: solver, 4;\n");
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"0.5 * pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert(N, Nref, tol * norm(Nref));
%!   assert(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 9
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20r", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20r"))
%!     mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!     mesh.materials.iso20r([mesh.groups.iso20r(find([[mesh.groups.iso20r.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fputs(fd, "         threads: assembly, 4;\n");
%!     fputs(fd, "         threads: solver, 4;\n");
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert(N, Nref, tol * norm(Nref));
%!   assert(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 10
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 1e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.h1 = 5e-3 / SI_unit_meter;
%! param.h2 = 5e-3 / SI_unit_meter;
%! param.z0 = 10e-3 / SI_unit_meter;
%! param.vz0 = -30 / (SI_unit_meter / SI_unit_second);
%! param.gz = -9.81 / (SI_unit_meter / SI_unit_second^2);
%! param.t1 = 0.001 / SI_unit_second;
%! param.N = 1000;
%! options.verbose = true;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Sphere(1) = {0, 0, z0, 0.5 * d1};\n");
%!     fputs(fd, "Physical Volume(1) = {1};\n");
%!     fputs(fd, "Physical Surface(1) = {1};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1};}} = h1;\n");
%!     fputs(fd, "MeshSize{PointsOf{Surface{1};}} = h2;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   node_id_ground = rows(mesh.nodes) + 1;
%!   mesh.nodes(node_id_ground, 1:3) = zeros(1, 3);
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.tet10h([mesh.groups.tet10h(find([[mesh.groups.tet10h.id] == 1])).elements]) = 1;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_id_ground) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_sphere";
%!   opt_mbd_mesh.joints.number = 1;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   node_id_cont = mesh.groups.tria6h(1).nodes;
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_sphere = 1000;\n");
%!     fputs(fd, "set: integer drive_id_E = 2000;\n");
%!     fprintf(fd, " set: real t1 = %g;\n", param.t1);
%!     fprintf(fd, " set: real gz = %g;\n", param.gz);
%!     fprintf(fd, " set: real z0 = %g;\n", param.z0);
%!     fprintf(fd, " set: real vz0 = %g;\n", param.vz0);
%!     fprintf(fd, " set: integer N = %d;\n", param.N);
%!     fprintf(fd, " set: integer node_id_ground = %d;\n", node_id_ground);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: pardiso, grad, max iterations, 100;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: mcp newton min fb;\n");
%!     fputs(fd, "         threads: assembly, 4;\n");
%!     fputs(fd, "         threads: solver, 4;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 100;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    loadable elements: %d;\n", numel(node_id_cont));
%!     fputs(fd, "      gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, "reference: ref_id_sphere,\n");
%!     fputs(fd, "  position, reference, global, null,\n");
%!     fputs(fd, "  orientation, reference, global, eye,\n");
%!     fputs(fd, "  velocity, reference, global, 0., 0., vz0,\n");
%!     fputs(fd, "  angular velocity, reference, global, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     for i=1:numel(node_id_cont)
%!       fprintf(fd, "user defined: %d, unilateral disp in plane,\n", i);
%!       fprintf(fd, "node1, %d,\n", node_id_cont(i));
%!       fputs(fd, "    offset, 1,\n");
%!       fputs(fd, "    reference, node, null,\n");
%!       fputs(fd, "    enable mcp, yes,\n");
%!       fputs(fd, " node2, node_id_ground,\n");
%!       fputs(fd, " offset, reference, node, null,\n");
%!       fputs(fd, " hinge, 3, 0, 0, 1,\n");
%!       fputs(fd, "        1, 1, 0, 0;\n");
%!     endfor
%!     fputs(fd, "  joint: 1, clamp, node_id_ground, node, node;\n");
%!     fputs(fd, "  gravity: uniform, component, 0., 0., gz;\n");
%!     fprintf(fd, "drive caller: drive_id_E, array, %d", rows(mesh.elements.tet10h));
%!     for i=1:rows(mesh.elements.tet10h)
%!       fprintf(fd, ",\n    element, %d, solid, string, \"E\", direct", i);
%!     endfor
%!     fputs(fd, ",\n  output, yes;\n\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   #shell(sprintf("cat \"%s\"", nodes_file));
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [drive_id, drive_data] = mbdyn_post_load_output_drv(opt_mbdyn.output_file);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! close all;
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!             0,  0.5 * b,  0.5 * c;  ##  2
%!             0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!             0,  0.5 * b, -0.5 * c;  ##  6
%!             0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!             a,  0.5 * b,  0.5 * c;  ##  9
%!             a, -0.5 * b,  0.5 * c;  ## 10
%!             a,  0.5 * b, -0.5 * c;  ## 11
%!             a, -0.5 * b, -0.5 * c,  ## 12
%!         a + d,        0,        0;  ## 13
%!            -e,        0,        0]; ## 14
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.element.name = "elem_id_cube";
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.modal.name = "node_id_cube_modal";
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.nodes.interfaces.name = "node_id_cube_interface1";
%! cms_opt.number_of_threads = 1;
%! cms_opt.algorithm = "eliminate";
%! cms_opt.invariants = true;
%! [mesh_cms, ...
%!  mat_ass_cms, ...
%!  dof_map_cms, ...
%!  sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
