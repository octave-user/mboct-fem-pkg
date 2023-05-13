## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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

%!test
%! ## Rotating thin ring
%! ## W.Beitz, K.-H.Grote, 1997
%! ## Dubbel page C41, chapter 6.2
%! close all;
%! ## Define the unit system for the solver
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second); ## Rayleigh damping mass coefficient
%! param.beta = 1e-7 / (SI_unit_second); ## Rayleigh damping stiffness coefficient
%! param.h = 1e-3 / SI_unit_meter; ## mesh size
%! param.w = 2e-3 / SI_unit_meter; ## width of disk
%! param.d = 48e-3 / SI_unit_meter; ## inner diameter disk
%! param.D = 50e-3 / SI_unit_meter; ## outer diameter disk
%! param.E = 70000e6 / SI_unit_pascal; ## Young's modulus
%! param.nu = 0.3; ## Poisson number
%! param.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3); ## density
%! param.omega0 = 15000 * pi / 30 / (1 / SI_unit_second); ## angular velocity
%! param.omega1 = 15000 * pi / 30 / (1 / SI_unit_second);
%! param.n = 1; ## number of revolutions

%! cms_opt.algorithm = "shift-invert";
%! cms_opt.refine_max_iter = int32(10);
%! cms_opt.element.name = "elem_id_disk";
%! cms_opt.nodes.modal.name = "node_id_inner_diameter";
%! cms_opt.nodes.interfaces(1).name = "node_id_outer_diameter";
%! cms_opt.number_of_threads = int32(4);
%! cms_opt.verbose = false;

%! ## It is very important to have at least one mode shape which accounts for centrifugal loads!
%! cms_opt.load_cases = "all";
%! options.number_of_modes = int32(10);
%! options.scale_def = 10e-3;
%! options.geo_tol = sqrt(eps);
%! options.mbdyn_command = "mbdyn";
%! options.f_run_mbdyn = true;
%! options.verbose = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     geometry_file = [filename, ".geo"];
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %g;\n", fn{i}, getfield(param, fn{i}));
%!     endfor

%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0.5 * d, 0, 0.5 * w};\n");
%!     fputs(fd, "Point(2) = {0.5 * D, 0, 0.5 * w};\n");
%!     fputs(fd, "Point(3) = {0.5 * D, 0, -0.5 * w};\n");
%!     fputs(fd, "Point(4) = {0.5 * d, 0, -0.5 * w};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "vol1[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{1}; };\n");
%!     fputs(fd, "vol2[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol1[0]}; };\n");
%!     fputs(fd, "vol3[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol2[0]}; };\n");
%!     fputs(fd, "vol4[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol3[0]}; };\n");
%!     fputs(fd, "vol6 = newv;\n");
%!     fputs(fd, "BooleanUnion(vol6) = { Volume{vol1[1]}; Delete; }{ Volume{vol2[1], vol3[1], vol4[1]}; Delete; };\n");
%!     fputs(fd, "B[] = Unique(Abs(Boundary{Volume{vol6};}));\n");
%!     fputs(fd, "Physical Volume(\"V\", 1) = {vol6};\n");
%!     fputs(fd, "For i In {0:#B[] - 1}\n");
%!     fputs(fd, "  Physical Surface(i) = {B[i]};\n");
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
%!                        "-ho_min", "0.5", ...
%!                        "-ho_max", "1.5", ...
%!                        geometry_file, ...
%!                        "-o", [filename, ".msh"]});

%!   status = spawn_wait(pid);

%!   if status ~= 0
%!     warning("gmsh failed with status %d", status);
%!   endif

%!   fprintf(stderr, "loading mesh ...\n");
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));

%!   group_defs(1).id = 1;
%!   group_defs(1).name = "inner_diameter";
%!   group_defs(1).R = eye(3);
%!   group_defs(1).X0 = zeros(3, 1);
%!   group_defs(1).type = "cylinder";
%!   group_defs(1).geometry.rmin = 0.5 * param.d;
%!   group_defs(1).geometry.rmax = 0.5 * param.d;
%!   group_defs(1).geometry.zmin = -0.5 * param.w;
%!   group_defs(1).geometry.zmax = 0.5 * param.w;

%!   group_defs(2).id = 2;
%!   group_defs(2).name = "outer_diameter";
%!   group_defs(2).R = eye(3);
%!   group_defs(2).X0 = zeros(3, 1);
%!   group_defs(2).type = "cylinder";
%!   group_defs(2).geometry.rmin = 0.5 * param.D;
%!   group_defs(2).geometry.rmax = 0.5 * param.D;
%!   group_defs(2).geometry.zmin = -0.5 * param.w;
%!   group_defs(2).geometry.zmax = 0.5 * param.w;

%!   groups = fem_pre_mesh_groups_create(mesh, group_defs, options.geo_tol);

%!   mesh.groups.tria6 = groups.tria6;

%!   cms_opt.modes.number = int32(options.number_of_modes);
%!   cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!   cms_opt.nodes.interfaces(1).number = int32(rows(mesh.nodes) + 2);
%!   cms_opt.invariants = true;
%!   mesh.nodes(cms_opt.nodes.modal.number, :) = [0, 0, 0, 0, 0, 0];
%!   mesh.nodes([cms_opt.nodes.interfaces.number], :) = [0, 0, 0, 0, 0, 0];
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, cms_opt.nodes.modal.number);
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces(1).number);

%!   load_case(1).locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));

%!   for i=1:numel(mesh.groups.tria6)
%!     ## In order to account for the centrifugal load in this demo,
%!     ## we need at least one load case which applies a force in radial direction.
%!     ## This is done here by applying a virtual pressure load.
%!     load_case(end + 1).pressure.tria6.elements = mesh.elements.tria6(mesh.groups.tria6(i).elements, :);
%!     load_case(end).pressure.tria6.p = ones(size(load_case(end).pressure.tria6.elements));
%!   endfor

%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(find([mesh.groups.tet10.id == 1])).elements) = 1;
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(param.E, param.nu);

%!   fprintf(stderr, "building cms element ...\n");
%!   [cms_data.mesh, cms_data.mat_ass, cms_data.dof_map, cms_data.sol_eig, cms_data.cms_opt] = fem_cms_create(mesh, load_case, cms_opt);

%!   for i=1:numel(cms_data.sol_eig.f)
%!     fprintf(stderr, "mode %d: %.1fHz\n", i, cms_data.sol_eig.f(i));
%!   endfor

%!   mesh_post_pro_file = sprintf("%s_post.msh", filename);

%!   fem_post_mesh_export(mesh_post_pro_file, cms_data.mesh);

%!   for j=1:numel(cms_data.sol_eig.f)
%!     eig_post_pro_file_mode{j} = sprintf("%s_eig_def_%03d.msh", filename, j);
%!     fem_post_sol_step_export(eig_post_pro_file_mode{j}, cms_data.sol_eig, j, j, cms_data.sol_eig.f(j), options.scale_def / max(norm(cms_data.sol_eig.def(:, 1:3, j), "rows")));
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

%!   cms_data.mat_ass.Dred = param.alpha * cms_data.mat_ass.Mred + param.beta * cms_data.mat_ass.Kred;
%!   fem_cms_export([filename, "_cms"], cms_data.mesh, cms_data.dof_map, cms_data.mat_ass, cms_opt);
%!   options_mbd.output_file = filename;
%!   options_mbd.mbdyn_command = options.mbdyn_command;

%!   options_mbd.logfile = [filename, ".stdout"];
%!   options_mbd.f_run_mbdyn2easyanim = false;
%!   param_file = [filename, ".set"];
%!   putenv("MBDYN_ROT_DISK_CMS_ELEM_FILE", [filename, "_cms.elm"]);
%!   putenv("MBDYN_ROT_DISK_CMS_PARAM_FILE", param_file);
%!   mbdyn_pre_write_param_file(param_file, param);
%!   mbdyn_filename = [filename, ".mbdyn"];

%!   fd = -1;

%!   unwind_protect
%!     fd = fopen(mbdyn_filename, "w");

%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_filename);
%!     endif

%!     fputs(fd, "include: \"${MBDYN_ROT_DISK_CMS_PARAM_FILE}\";\n");

%!     fputs(fd, "set: integer ref_id_disk = 1001;\n");
%!     fputs(fd, "set: integer node_id_inner_diameter = 2001;\n");
%!     fputs(fd, "set: integer node_id_outer_diameter = 2002;\n");
%!     fputs(fd, "set: integer elem_id_disk = 3001;\n");
%!     fputs(fd, "set: integer joint_id_drive = 3002;\n");
%!     fputs(fd, "set: integer joint_id_bearing = 3003;\n");
%!     fputs(fd, "set: integer elem_id_inertia = 3004;\n");
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
%!     fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always;\n");
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
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output frequency: 10;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");

%!     fputs(fd, "        structural nodes:\n");
%!     fputs(fd, "                +1		# modal\n");
%!     fputs(fd, "                +1		# interface 1\n");
%!     fputs(fd, "        ;\n");
%!     fputs(fd, "        joints:\n");
%!     fputs(fd, "                +1		# modal\n");
%!     fputs(fd, "                +1              # drive\n");
%!     fputs(fd, "                +1              # pin\n");
%!     fputs(fd, "        ;\n");
%!     fputs(fd, "end: control data;\n");

%!     fputs(fd, "reference: ref_id_disk,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, 0., 0., omega0;\n");

%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_inner_diameter, modal,\n");
%!     fputs(fd, "                reference, ref_id_disk, null,\n");
%!     fputs(fd, "                reference, ref_id_disk, eye,\n");
%!     fputs(fd, "                reference, ref_id_disk, null,\n");
%!     fputs(fd, "                reference, ref_id_disk, null;\n");

%!     fputs(fd, "        structural: node_id_outer_diameter, static,\n");
%!     fputs(fd, "                reference, ref_id_disk, null,\n");
%!     fputs(fd, "                reference, ref_id_disk, eye,\n");
%!     fputs(fd, "                reference, ref_id_disk, null,\n");
%!     fputs(fd, "                reference, ref_id_disk, null;\n");
%!     fputs(fd, "end: nodes;\n");

%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\";\n");
%!     fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");

%!     fputs(fd, "	joint: joint_id_drive, angular velocity,\n");
%!     fputs(fd, "		# node label\n");
%!     fputs(fd, "		node_id_inner_diameter,\n");
%!     fputs(fd, "		# direction\n");
%!     fputs(fd, "		0.,0.,1.,\n");
%!     fputs(fd, "		# angular velocity\n");
%!     fputs(fd, "		reference, drive_id_rotor_speed;\n");

%!     fputs(fd, "        joint: joint_id_bearing, revolute pin,\n");
%!     fputs(fd, "               node_id_inner_diameter,\n");
%!     fputs(fd, "               position, reference, ref_id_disk, null,\n");
%!     fputs(fd, "               orientation, reference, ref_id_disk, eye,\n");
%!     fputs(fd, "               position, reference, ref_id_disk, null,\n");
%!     fputs(fd, "               orientation, reference, ref_id_disk, eye;\n");

%!     fputs(fd, "        include: \"${MBDYN_ROT_DISK_CMS_ELEM_FILE}\";\n");

%!     fputs(fd, "        inertia: elem_id_inertia,\n");
%!     fputs(fd, "                 position, reference, ref_id_disk, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_disk, eye,\n");
%!     fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!     fputs(fd, "                 output, both;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect

%!   mbdyn_solver_run(mbdyn_filename, options_mbd);

%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);

%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(options_mbd.output_file);

%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);

%!   res.sol_dyn = fem_post_cms_sol_import(options_mbd.output_file, cms_data);

%!   cms_data_p.sol_tran.def = fem_post_cms_expand_body(cms_data.mesh, ...
%!                                                      cms_data.dof_map, ...
%!                                                      cms_data.mat_ass, ...
%!                                                      res.sol_dyn.bodies(1).q);

%!   cms_data_p.sol_tran.stress = fem_ass_matrix(cms_data.mesh, ...
%!                                               cms_data.dof_map, ...
%!                                               [FEM_VEC_STRESS_CAUCH], ...
%!                                               load_case, ...
%!                                               cms_data_p.sol_tran);

%!   R = 0.5 * param.D;
%!   omega = param.omega1;
%!   rho = param.rho;
%!   E = param.E;
%!   TAU_t_ref = diag([0; rho * omega^2 * R^2; 0]);
%!   U_r_ref = [rho * omega^2 * R^3 / E; 0; 0];
%!   idxtens = int32([1, 4, 6;
%!                    4, 2, 5;
%!                    6, 5, 3]);
%!   tol = 7e-2;
%!
%!   inode = cms_data.mesh.elements.tet10(1, 1);
%!   x = cms_data.mesh.nodes(inode, 1:3).';
%!   e1 = x / norm(x(1:2));
%!   e3 = [0; 0; 1];
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   tau = reshape(cms_data_p.sol_tran.stress.taum.tet10(1, 1, :, :), 6, numel(res.t));
%!   U = reshape(cms_data_p.sol_tran.def(inode, 1:3, :), 3, numel(res.t));
%!   u_r = zeros(1, numel(res.t));
%!   tau_t = zeros(1, numel(res.t));
%!   for i=1:numel(res.t)
%!     TAU = tau(:, i)(idxtens);
%!     TAU_t = R.' * TAU * R;
%!     tau_t(i) = TAU_t(2, 2);
%!     U_t = R.' * U(:, i);
%!     u_r(i) = U_t(1);
%!   endfor

%!   figure("visible", "off");
%!   hold on;
%!   plot(res.t * SI_unit_second, tau_t * SI_unit_pascal, "-;FEM tau11;1");
%!   plot(res.t([1, end]) * SI_unit_second, TAU_t_ref(2, 2)([1, end]) * SI_unit_pascal, "-;reference tau11;0");
%!   xlabel("t [s]");
%!   ylabel("stress [Pa]");
%!   grid on;
%!   grid minor on;
%!   title("rotating thin ring - stress versus time");

%!   figure("visible", "off");
%!   hold on;
%!   plot(res.t * SI_unit_second, u_r * SI_unit_meter, "-;FEM tau11;1");
%!   plot(res.t([1, end]) * SI_unit_second, U_r_ref(1)([1, end]) * SI_unit_meter, "-;reference tau11;0");
%!   xlabel("t [s]");
%!   ylabel("radial deformation [m]");
%!   grid on;
%!   grid minor on;
%!   title("rotating thin ring - deformation versus time");

%!   for j=1:rows(cms_data_p.sol_tran.stress.taum.tet10)
%!     for k=1:columns(cms_data_p.sol_tran.stress.taum.tet10)
%!       inode = cms_data.mesh.elements.tet10(j, k);
%!       x = cms_data.mesh.nodes(inode, 1:3).';
%!       e1 = x / norm(x(1:2));
%!       e3 = [0; 0; 1];
%!       e2 = cross(e3, e1);
%!       R = [e1, e2, e3];
%!       tau = cms_data_p.sol_tran.stress.taum.tet10(j, k, :, end)(:);
%!       U = cms_data_p.sol_tran.def(inode, 1:3, end)(:);
%!       TAU = tau(idxtens);
%!       TAU_t = R.' * TAU * R;
%!       U_t = R.' * U;
%!       assert(TAU_t, TAU_t_ref, tol * norm(TAU_t_ref));
%!       assert(U_t, U_r_ref, tol * norm(U_r_ref));
%!     endfor
%!   endfor
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(fn))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       rc = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (rc ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
