## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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

  if (ispc())
    filename(find(filename == '\')) = '/'; ## Required for MBDyn's parser
  endif
  
  fd = -1;
  
  unwind_protect
    [fd, msg] = fopen([filename, ".elm"], "wt");

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
                             mat_ass.Inv9);
  else  
    mbdyn_pre_write_fem_data([filename, ".fem"], ...
                             mat_ass.Mred, ...
                             mat_ass.Dred, ...
                             mat_ass.Kred, ...
                             Phi, ...
                             mesh.nodes(:, 1:3).', ...
                             zeros(columns(mat_ass.Mred), 1), ...
                             zeros(columns(mat_ass.Mred), 1), ...
                             mat_ass.diagM);
  endif
  
  if (cms_opt.verbose)
    fprintf(stderr, "file \"%s\" created ...\n", [filename, ".fem"]);
    toc();
  endif
endfunction

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
%! mesh.material_data.E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
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
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
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
%!     fd = fopen(mbdyn_filenames{1}, "wt");
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
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
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
%!     fd = fopen(mbdyn_filenames{2}, "wt");
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
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
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
%!     [fd, msg] = fopen(geometry_file, "wt");

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

%!   group_defs(end + 1).id = 1;
%!   group_defs(end).name = "bearing1";
%!   group_defs(end).R = eye(3);
%!   group_defs(end).X0 = [0; 0; -0.5 * param.l];
%!   group_defs(end).type = "box";
%!   group_defs(end).geometry.xmin = -0.5 * param.d;
%!   group_defs(end).geometry.xmax = 0.5 * param.d;
%!   group_defs(end).geometry.ymin = -0.5 * param.d;
%!   group_defs(end).geometry.ymax = 0.5 * param.d;
%!   group_defs(end).geometry.zmin = 0;
%!   group_defs(end).geometry.zmax = 0;

%!   group_defs(end + 1).id = 2;
%!   group_defs(end).name = "bearing2";
%!   group_defs(end).R = eye(3);
%!   group_defs(end).X0 = [0; 0; 0.5 * param.l];
%!   group_defs(end).type = "box";
%!   group_defs(end).geometry.xmin = -0.5 * param.d;
%!   group_defs(end).geometry.xmax = 0.5 * param.d;
%!   group_defs(end).geometry.ymin = -0.5 * param.d;
%!   group_defs(end).geometry.ymax = 0.5 * param.d;
%!   group_defs(end).geometry.zmin = 0;
%!   group_defs(end).geometry.zmax = 0;

%!   group_defs(end + 1).id = 3;
%!   group_defs(end).name = "rotor";
%!   group_defs(end).R = eye(3);
%!   group_defs(end).X0 = [0; 0; param.o];
%!   group_defs(end).type = "cylinder";
%!   group_defs(end).geometry.rmin = 0.5 * param.D;
%!   group_defs(end).geometry.rmax = 0.5 * param.D;
%!   group_defs(end).geometry.zmin = -0.5 * param.w;
%!   group_defs(end).geometry.zmax = 0.5 * param.w;
  
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
%!   mesh.material_data.E = param.E;
%!   mesh.material_data.nu = param.nu;
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);

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
%!     [fd, msg] = fopen(eig_post_pro_file, "wt");

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
