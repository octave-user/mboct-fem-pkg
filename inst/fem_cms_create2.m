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
## @deftypefn {Function File} [@var{mesh}, @var{mat_ass}, @var{dof_map}, @var{sol_eig}, @var{cms_opt}, @var{sol_tau}] = fem_cms_create2(@var{mesh}, @var{load_case}, @var{cms_opt})
## Build a reduced order model by using the Craig Bampton approach. That model can be exported to MBDyn by means of fem_cms_export.
##
## @var{mesh} @dots{} Finite element mesh data structure containing constraints for the modal node
##
## @var{load_case}.locked_dof @dots{} Boolean matrix of constrained degrees of freedom
##
## @var{cms_opt}.nodes.modal @dots{} Struct containing node number and node name of the modal node. Appropriate constraints must be defined for that node in @var{mesh} or @var{load_case}
##
## @var{cms_opt}.nodes.interfaces @dots{} Struct array containing node numbers and node names for interface nodes accessible to MBDyn. Static mode shapes will be generated for those nodes.
##
## @var{cms_opt}.verbose @dots{} Enable verbose output
##
## @var{cms_opt}.algorithm @dots{} Algorithm used for eigenanalysis
##
## @var{cms_opt}.number_of_threads @dots{} Number of threads used for linear solver and assembly
##
## @var{cms_opt}.refine_max_iter @dots{} Maximum number of refinement iterations for the linear solver
##
## @seealso{fem_cms_export}
## @end deftypefn

function [mesh, mat_ass_itf, dof_map_itf, sol_eig_itf, cms_opt, sol_tau] = fem_cms_create2(mesh, load_case_dof, cms_opt)
  if (~isfield(cms_opt, "modes"))
    cms_opt.modes = struct();
  endif

  if (~isfield(cms_opt.modes, "number"))
    cms_opt.modes.number = int32(0);
  endif

  if (~isfield(cms_opt, "number_of_threads"))
    cms_opt.number_of_threads = int32(1);
  endif

  if (~isfield(cms_opt, "threshold_elem"))
    cms_opt.threshold_elem = int32(10000);
  endif

  if (~isfield(cms_opt, "tolerance_tau"))
    cms_opt.tolerance_tau = -1;
  endif

  if (~isfield(cms_opt, "algorithm"))
    cms_opt.algorithm = "shift-invert";
  endif

  if (~isfield(cms_opt, "enable_KTAU0W"))
    cms_opt.enable_KTAU0W = true(6, 1);
  endif

  if (~isfield(cms_opt, "enable_KTAU0WP"))
    cms_opt.enable_KTAU0WP = true(3, 1);
  endif

  if (~isfield(cms_opt, "enable_KTAU0VP"))
    cms_opt.enable_KTAU0VP = true(3, 1);
  endif

  for i=1:numel(cms_opt.nodes.interfaces)
    if (~isfield(cms_opt.nodes.interfaces(i), "enable_KTAU0") || isempty(cms_opt.nodes.interfaces(i).enable_KTAU0))
      cms_opt.nodes.interfaces(i).enable_KTAU0 = true(6, 1);
    endif
  endfor

  ## required for angular velocity and angular acceleration loads
  mesh.nodes -= mesh.nodes(cms_opt.nodes.modal.number, :);

  if (~isfield(cms_opt, "invariants"))
    cms_opt.invariants = true;
  endif

  joint_modal.C = eye(6)(~load_case_dof.locked_dof(cms_opt.nodes.modal.number, :), :);
  joint_modal.nodes = cms_opt.nodes.modal.number;

  if (~isfield(mesh.elements, "joints"))
    joint_offset = int32(0);
    mesh.elements.joints = struct("C", cell(0, 1), "nodes", cell(0, 1));
  else
    joint_offset = numel(mesh.elements.joints);
  endif

  if (rows(joint_modal.C))
    mesh.elements.joints(++joint_offset) = joint_modal;
  endif

  dof_map_tau = fem_ass_dof_map(mesh, load_case_dof);
  dof_map_tau.parallel.threads_ass = cms_opt.number_of_threads;
  dof_map_tau.parallel.threshold_elem = cms_opt.threshold_elem;

  nomega = sum(cms_opt.enable_KTAU0W);
  nomega_dot = sum(cms_opt.enable_KTAU0WP);
  ng = sum(cms_opt.enable_KTAU0VP);

  num_loads_KTAU0 = nomega + nomega_dot + ng;

  for i=1:numel(cms_opt.nodes.interfaces)
    num_loads_KTAU0 += sum(dof_map_tau.ndof(cms_opt.nodes.interfaces(i).number, cms_opt.nodes.interfaces(i).enable_KTAU0) > 0);
  endfor

  cms_opt.index_KTAU0red = zeros(num_loads_KTAU0, 1, "int32");

  load_case_tau = fem_pre_load_case_create_empty(num_loads_KTAU0);

  for i=1:numel(load_case_tau)
    load_case_tau(i).omegaq = zeros(6, 1);
    load_case_tau(i).omegadot = zeros(3, 1);
    load_case_tau(i).g = zeros(3, 1);
  endfor

  idx = int32(0);
  id = int32(0);

  for i=1:6
    ++id;

    if (cms_opt.enable_KTAU0W(i))
      load_case_tau(++idx).omegaq(i) = 1;
      cms_opt.index_KTAU0red(idx) = id;
    endif
  endfor

  for i=1:3
    ++id;

    if (cms_opt.enable_KTAU0WP(i))
      load_case_tau(++idx).omegadot(i) = 1;
      cms_opt.index_KTAU0red(idx) = id;
    endif
  endfor

  for i=1:3
    ++id;

    if (cms_opt.enable_KTAU0VP(i))
      load_case_tau(++idx).g(i) = -1; ## By convention MBDyn uses (vP - g)
      cms_opt.index_KTAU0red(idx) = id;
    endif
  endfor

  for i=1:numel(cms_opt.nodes.interfaces)
    for j=1:columns(dof_map_tau.ndof)
      ++id;

      if (~(dof_map_tau.ndof(cms_opt.nodes.interfaces(i).number, j) > 0))
        continue;
      endif

      if (cms_opt.nodes.interfaces(i).enable_KTAU0(j))
        ++idx;
        load_case_tau(idx).loaded_nodes = cms_opt.nodes.interfaces(i).number;
        load_case_tau(idx).loads = zeros(1, 6);
        load_case_tau(idx).loads(j) = 1;
        cms_opt.index_KTAU0red(idx) = id;
      endif
    endfor
  endfor

  [mat_ass_tau.K, ...
   mat_ass_tau.R, ...
   mat_ass_tau.mat_info, ...
   mat_ass_tau.mesh_info] = fem_ass_matrix(mesh, ...
                                           dof_map_tau, ...
                                           [FEM_MAT_STIFFNESS_SYM_L, ...
                                            FEM_VEC_LOAD_CONSISTENT], ...
                                           load_case_tau);

  [sol_tau, U_tau] = fem_sol_static(mesh, dof_map_tau, mat_ass_tau, cms_opt);

  KTAU0 = cell(1, size(sol_tau.def, 3));

  for i=1:size(sol_tau.def, 3)
    sol_tau.stress = fem_ass_matrix(mesh, ...
                                    dof_map_tau, ...
                                    [FEM_VEC_STRESS_CAUCH], ...
                                    load_case_tau, ...
                                    struct("def", sol_tau.def(:, :, i)));
    KTAU0{i} = fem_ass_matrix(mesh, ...
                              dof_map_tau, ...
                              [FEM_MAT_STIFFNESS_TAU0], ...
                              struct("tau0", sol_tau.stress.tau));
  endfor

  num_joints_itf = numel(cms_opt.nodes.interfaces);

  joints_itf = struct("C", cell(1, num_joints_itf), "nodes", cell(1, num_joints_itf));
  nodes_itf = int32([[cms_opt.nodes.interfaces].number]);
  idx_joint_itf = int32(0);

  for i=1:num_joints_itf
    f_active_dof = dof_map_tau.ndof(nodes_itf(i), :) > 0;

    if (~any(f_active_dof))
      continue;
    endif

    ++idx_joint_itf;
    joints_itf(idx_joint_itf).C = eye(6)(f_active_dof, :);
    joints_itf(idx_joint_itf).nodes = nodes_itf(i);
  endfor

  mesh.elements.joints(joint_offset + (1:idx_joint_itf)) = joints_itf(1:idx_joint_itf);

  dof_map_itf = fem_ass_dof_map(mesh, load_case_dof);
  dof_map_itf.parallel.threads_ass = cms_opt.number_of_threads;
  dof_map_itf.parallel.threshold_elem = cms_opt.threshold_elem;

  num_loads_itf = int32(0);

  for i=1:numel(cms_opt.nodes.interfaces)
    num_loads_itf += sum(dof_map_tau.ndof(cms_opt.nodes.interfaces(i).number, :) > 0);
  endfor

  load_case_itf = fem_pre_load_case_create_empty(num_loads_itf);

  for i=1:numel(load_case_itf)
    for j=1:numel(mesh.elements.joints)
      load_case_itf(i).joints(j).U = zeros(rows(mesh.elements.joints(j).C), 1);
    endfor
  endfor

  idx = int32(0);

  for i=1:numel(nodes_itf)
    k = int32(0);

    for j=1:6
      if (~(dof_map_itf.ndof(nodes_itf(i), j) > 0))
        continue;
      endif

      load_case_itf(++idx).joints(i + joint_offset).U(++k) = 1;
    endfor
  endfor

  if (cms_opt.invariants)
    [mat_ass_itf.M, ...
     mat_ass_itf.D, ...
     mat_ass_itf.K, ...
     mat_ass_itf.R, ...
     mat_ass_itf.dm, ...
     mat_ass_itf.S, ...
     mat_ass_itf.J, ...
     mat_ass_itf.mat_info, ...
     mat_ass_itf.mesh_info] = fem_ass_matrix(mesh, ...
                                             dof_map_itf, ...
                                             [FEM_MAT_MASS_SYM_L, ...
                                              FEM_MAT_DAMPING_SYM_L, ...
                                              FEM_MAT_STIFFNESS_SYM_L, ...
                                              FEM_VEC_LOAD_CONSISTENT, ...
                                              FEM_SCA_TOT_MASS, ...
                                              FEM_VEC_INERTIA_M1, ...
                                              FEM_MAT_INERTIA_J], ...
                                             load_case_itf);
  else
    [mat_ass_itf.M, ...
     mat_ass_itf.D, ...
     mat_ass_itf.K, ...
     mat_ass_itf.R, ...
     mat_ass_itf.Mlumped, ...
     mat_ass_itf.mat_info, ...
     mat_ass_itf.mesh_info] = fem_ass_matrix(mesh, ...
                                             dof_map_itf, ...
                                             [FEM_MAT_MASS_SYM_L, ...
                                              FEM_MAT_DAMPING_SYM_L, ...
                                              FEM_MAT_STIFFNESS_SYM_L, ...
                                              FEM_VEC_LOAD_CONSISTENT, ...
                                              FEM_MAT_MASS_LUMPED], ...
                                             load_case_itf);

    diagMlumped = full(diag(mat_ass_itf.Mlumped));
    mat_ass_itf.diagM = zeros(6 * rows(mesh.nodes), 1);

    ridx_node(dof_map_itf.idx_node) = 1:numel(dof_map_itf.idx_node);

    for i=1:columns(dof_map_itf.ndof)
      idx_act_dof = find(dof_map_itf.ndof(:, i) > 0);
      idx_glob_dof = dof_map_itf.ndof(idx_act_dof, i);
      mat_ass_itf.diagM((idx_act_dof - 1) * 6 + i) = diagMlumped(idx_glob_dof);
    endfor
  endif

  Msym = fem_mat_sym(mat_ass_itf.M);
  Ksym = fem_mat_sym(mat_ass_itf.K);
  Dsym = fem_mat_sym(mat_ass_itf.D);

  mat_ass_itf.Tred =  zeros(numel(dof_map_itf.idx_node), columns(mat_ass_itf.R) + cms_opt.modes.number);

  Kfact = fem_sol_factor(mat_ass_itf.K, cms_opt);

  if (columns(mat_ass_itf.R))
    mat_ass_itf.Tred(:, cms_opt.modes.number + (1:columns(mat_ass_itf.R))) = (Kfact \ mat_ass_itf.R)(dof_map_itf.idx_node, :);
  endif

  switch (cms_opt.algorithm)
    case {"shift-invert", "diag-shift-invert"}
    otherwise
      error("algorithm \"%s\" not supported", cms_opt.algorithm);
  endswitch

  if (cms_opt.modes.number)
    rndstate = rand("state");

    unwind_protect
      rand("seed", 0);

      SIGMA = 0;
      op{1} = @(x) Msym * x;
      op{2} = @(x) Kfact \ x;

      opt_eig.disp = 0;
      opt_eig.maxit = 50000;
      opt_eig.tol = 0;

      [PHI_d, mu] = eig_sym(op, columns(Msym), cms_opt.modes.number, SIGMA, opt_eig);
      mat_ass_itf.Tred(:, 1:cms_opt.modes.number) = PHI_d(dof_map_itf.idx_node, :);

      lambda = sqrt(-diag(mu)).';
      clear op PHI_d;
    unwind_protect_cleanup
      rand("state", rndstate);
    end_unwind_protect
  else
    lambda = [];
  endif

  if (cms_opt.tolerance_tau >= 0)
    Phi_tau = U_tau(dof_map_tau.idx_node, 1:(nomega + nomega_dot + ng));
    Phi_tau_res = Phi_tau - mat_ass_itf.Tred * (mat_ass_itf.Tred \ Phi_tau);
    mat_ass_itf.Tred = [mat_ass_itf.Tred, Phi_tau_res(:, norm(Phi_tau_res, "cols") > cms_opt.tolerance_tau * norm(Phi_tau, "cols"))];
  endif

  mat_ass_itf.Mred = fem_cms_matrix_trans(mat_ass_itf.Tred, Msym(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");
  mat_ass_itf.Kred = fem_cms_matrix_trans(mat_ass_itf.Tred, Ksym(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");

  switch (cms_opt.algorithm)
    case "diag-shift-invert"
      [PHI_diag, lambda_diag] = eig(mat_ass_itf.Kred, mat_ass_itf.Mred, "chol", "vector");
      mat_ass_itf.Mred = PHI_diag.' * mat_ass_itf.Mred * PHI_diag;
      mat_ass_itf.Kred = PHI_diag.' * mat_ass_itf.Kred * PHI_diag;
      mat_ass_itf.Tred *= PHI_diag;

      clear PHI_diag lambda_diag;
  endswitch

  mat_ass_itf.Dred = fem_cms_matrix_trans(mat_ass_itf.Tred, Dsym(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");

  mat_ass_itf.KTAU0red = zeros(columns(mat_ass_itf.Tred), columns(mat_ass_itf.Tred), numel(KTAU0));

  for i=1:numel(KTAU0)
    mat_ass_itf.KTAU0red(:, :, i) = fem_cms_matrix_trans(mat_ass_itf.Tred, KTAU0{i}(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");
  endfor

  S = diag(1 ./ sqrt(abs(diag(mat_ass_itf.Mred))));

  mat_ass_itf.Mred = S * mat_ass_itf.Mred * S;
  mat_ass_itf.Kred = S * mat_ass_itf.Kred * S;
  mat_ass_itf.Dred = S * mat_ass_itf.Dred * S;
  mat_ass_itf.Tred *= S;

  for i=1:numel(KTAU0)
    mat_ass_itf.KTAU0red(:, :, i) = S * mat_ass_itf.KTAU0red(:, :, i) * S;
  endfor

  if (nargout >= 4)
    sol_eig_itf.def = zeros(rows(mesh.nodes), columns(mesh.nodes), columns(mat_ass_itf.Tred));
    PHI = zeros(dof_map_itf.totdof, 1);

    for i=1:columns(mat_ass_itf.Tred)
      PHI(dof_map_itf.idx_node) = mat_ass_itf.Tred(:, i);
      sol_eig_itf.def(:, :, i) = fem_post_def_nodal(mesh, dof_map_itf, PHI);
    endfor

    clear PHI;
    sol_eig_itf.f = [imag(lambda) / (2 * pi), repmat(-inf, 1, columns(mat_ass_itf.R) + columns(U_tau))];
  endif

  if (cms_opt.invariants)
    [mat_ass_itf.Inv3, ...
     mat_ass_itf.Inv4, ...
     mat_ass_itf.Inv5, ...
     mat_ass_itf.Inv8, ...
     mat_ass_itf.Inv9] = fem_ass_matrix(mesh, ...
                                        dof_map_itf, ...
                                        [FEM_MAT_INERTIA_INV3, ...
                                         FEM_MAT_INERTIA_INV4, ...
                                         FEM_MAT_INERTIA_INV5, ...
                                         FEM_MAT_INERTIA_INV8, ...
                                         FEM_MAT_INERTIA_INV9], ...
                                        load_case_itf, ...
                                        sol_eig_itf);
  endif
endfunction

%!test
%! ## TEST1
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   a = 800e-3 / SI_unit_meter;
%!   b = 40e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   h = c;
%!   options.interactive = false;
%!   options.plot = true;
%!   options.verbose = false;
%!   options.number_of_beams = int32(40);
%!   options.number_of_threads = int32(4);
%!   if (options.plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fprintf(fd, "a = %.16e;\n", a);
%!     fprintf(fd, "b = %.16e;\n", b);
%!     fprintf(fd, "c = %.16e;\n", c);
%!     fprintf(fd, "h = %.16e;\n", h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(3) = {a,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(4) = {a, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,c}{Surface{1}; Layers{Max(1, Round(c / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {v1[2]};\n");
%!     fputs(fd, "Physical Surface(2) = {v1[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   mesh.material_data.E = 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data.alpha = 0e-5 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0e-5 / (SI_unit_second);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "diag-shift-invert";
%!   cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!   cms_opt.nodes.modal.name = "node_id_modal";
%!   cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!   cms_opt.nodes.interfaces.name = "node_id_interface1";
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [0, 0, 0];
%!   mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces.number, "quad8");
%!   cms_opt.refine_max_iter = 30;
%!   cms_opt.pre_scaling = false;
%!   cms_opt.solver = "pardiso";
%!   cms_opt.modes.number = 20;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   #cms_opt.enable_KTAU0WP = [true; false(2, 1)];
%!   #cms_opt.enable_KTAU0VP = [true; false(2, 1)];
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, :) = true;
%!   [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms, cms_opt, sol_tau_cms] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%!   pert.omega = [1e2; 3e2; 2e2] / (SI_unit_rad / SI_unit_second);
%!   pert.omegadot = [1e5; 1e3; 3e3] / (SI_unit_rad / SI_unit_second^2);
%!   pert.loads = [[1e4; 1e3; 1e2] / (SI_unit_newton);
%!                 [1e2; 1e1; 1e1] / (SI_unit_newton * SI_unit_meter)];
%!   pert.g = [1e4; -1e3; -1e2] / (SI_unit_meter / SI_unit_second^2);
%!   pert.a = [-1e4; 1e3; 1e2] / (SI_unit_meter / SI_unit_second^2);
%!   empty_cell = cell(7, 3, 2);
%!   res = struct("info", empty_cell, ...
%!                "t", empty_cell, ...
%!                "trajectory", empty_cell, ...
%!                "deformation", empty_cell, ...
%!                "velocity", empty_cell, ...
%!                "acceleration", empty_cell, ...
%!                "node_id", empty_cell, ...
%!                "force", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "orientation_description", empty_cell, ...
%!                "drive_id", empty_cell, ...
%!                "drive_value", empty_cell, ...
%!                "modal", empty_cell);
%!   empty_cell = cell(7, 3);
%!   param = struct("omega", empty_cell, ...
%!                  "omegadot", empty_cell, ...
%!                  "F1", empty_cell, ...
%!                  "M1", empty_cell, ...
%!                  "a", empty_cell, ...
%!                  "g", empty_cell, ...
%!                  "t1", empty_cell, ...
%!                  "holonomic", empty_cell);
%!   idx_j = 1:rows(param);
%!   idx_k = 1:columns(param);
%!   for j=idx_j
%!     for k=idx_k
%!       param(j, k).omega = zeros(3, 1);
%!       param(j, k).omegadot = zeros(3, 1);
%!       param(j, k).F1 = zeros(3, 1);
%!       param(j, k).M1 = zeros(3, 1);
%!       param(j, k).a = zeros(3, 1);
%!       param(j, k).g = zeros(3, 1);
%!       param(j, k).holonomic = false;
%!       param(j, k).gamma = zeros(3, 1);
%!       param(j, k).N = 50;
%!       switch (j)
%!         case 1
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).gamma = [20; 45; 30] * pi / 180;
%!         case 2
%!           param(j, k).omega(k) = pert.omega(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           param(j, k).holonomic = true;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 2000;
%!         case 3
%!           param(j, k).omegadot(k) = pert.omegadot(k);
%!           param(j, k).t1 = 1e-2 / SI_unit_second;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 200;
%!         case 4
%!           param(j, k).F1(k) = pert.loads(k);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 5
%!           param(j, k).M1(k) = pert.loads(k + 3);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           param(j, k).gamma(1) = 30 * pi / 180;
%!         case 6
%!           param(j, k).a(k) = pert.a(k);
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).N = 200;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 7
%!           param(j, k).g(k) = pert.g(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!       endswitch
%!       for l=1:2
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%d_%d_%d.mbdyn", filename, j, k, l);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fprintf(fd, "set: real a = %.16e;\n", a);
%!           fprintf(fd, "set: real b = %.16e;\n", b);
%!           fprintf(fd, "set: real c = %.16e;\n", c);
%!           fprintf(fd, "set: real d = %.16e;\n", d);
%!           for i=1:3
%!             fprintf(fd, "set: real gamma%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).gamma(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGA%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omega(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGAP%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omegadot(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real F1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).F1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real M1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).M1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real a%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).a(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).g(i));
%!           endfor
%!           fprintf(fd, "set: real t1 = %.16e;\n", param(j, k).t1);
%!           fprintf(fd, "set: integer N = %d;\n", param(j, k).N);
%!           fprintf(fd, "set: integer M = %d;\n", options.number_of_beams);
%!           fputs(fd, "set: integer ref_id_ground = 1;\n");
%!           fputs(fd, "set: integer ref_id_tilt = 2;\n");
%!           fputs(fd, "set: integer joint_id_ground = 1;\n");
%!           fputs(fd, "set: integer force_id1;\n");
%!           fputs(fd, "set: integer torque_id1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_PHI1 = 1;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP1 = 3;\n");
%!           fputs(fd, "set: integer drive_id_PHI2 = 4;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA2 = 5;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP2 = 6;\n");
%!           fputs(fd, "set: integer drive_id_PHIx = 7;\n");
%!           fputs(fd, "set: integer drive_id_PHIy = 8;\n");
%!           fputs(fd, "set: integer drive_id_PHIz = 9;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAx = 10;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAy = 11;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAz = 12;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPx = 13;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPy = 14;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPz = 15;\n");
%!           fputs(fd, "set: integer drive_id_Xx = 16;\n");
%!           fputs(fd, "set: integer drive_id_Xy = 17;\n");
%!           fputs(fd, "set: integer drive_id_Xz = 18;\n");
%!           fputs(fd, "set: integer drive_id_XPx = 19;\n");
%!           fputs(fd, "set: integer drive_id_XPy = 20;\n");
%!           fputs(fd, "set: integer drive_id_XPz = 21;\n");
%!           fputs(fd, "set: integer drive_id_XPPx = 22;\n");
%!           fputs(fd, "set: integer drive_id_XPPy = 23;\n");
%!           fputs(fd, "set: integer drive_id_XPPz = 24;\n");
%!           fputs(fd, "set: integer drive_id_gx = 25;\n");
%!           fputs(fd, "set: integer drive_id_gy = 26;\n");
%!           fputs(fd, "set: integer drive_id_gz = 27;\n");
%!           fputs(fd, "set: integer drive_id_F1x = 28;\n");
%!           fputs(fd, "set: integer drive_id_F1y = 29;\n");
%!           fputs(fd, "set: integer drive_id_F1z = 30;\n");
%!           fputs(fd, "set: integer drive_id_M1x = 31;\n");
%!           fputs(fd, "set: integer drive_id_M1y = 32;\n");
%!           fputs(fd, "set: integer drive_id_M1z = 33;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "set: integer ref_id_modal = 3;\n");
%!               fputs(fd, "set: integer node_id_modal = 1;\n");
%!               fputs(fd, "set: integer ref_id_interface1 = 4;\n");
%!               fputs(fd, "set: integer node_id_interface1 = 2;\n");
%!               fputs(fd, "set: integer elem_id_modal = 2;\n");
%!             case 2
%!               fputs(fd, "set: integer ref_id_beam1 = 3;\n");
%!               fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!               fprintf(fd, "set: real E = %.16e;\n", mesh.material_data.E);
%!               fprintf(fd, "set: real nu = %.16e;\n", mesh.material_data.nu);
%!               fprintf(fd, "set: real rho = %.16e;\n", mesh.material_data.rho);
%!               fprintf(fd, "set: real alpha = %.16e;\n", mesh.material_data.alpha);
%!               fprintf(fd, "set: real beta = %.16e;\n", mesh.material_data.beta);
%!               fputs(fd, "set: real G = E / (2. * (1. + nu));\n");
%!               fputs(fd, "set: real A = b * c;\n");
%!               fputs(fd, "set: real As = 9. / 10. * A;\n");
%!               fputs(fd, "set: real Iy = b * c^3 / 12.;\n");
%!               fputs(fd, "set: real Iz = c * b^3 / 12.;\n");
%!               fputs(fd, "set: real Ip = Iy + Iz;\n");
%!               c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!               w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!               fprintf(fd, "set: real It = %.16e;\n", interp1(w_h, c2, max(c, b) / min(c, b)) * max(c, b) * min(c, b)^3);
%!           endswitch
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / N;\n");
%!           fputs(fd, "        max time step: t1 / N;\n");
%!           fputs(fd, "        min time step: t1 / N;\n");
%!           fputs(fd, "        method: ss4, 0.;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!           fputs(fd, "        strategy: factor, 0.8, 3, 1.25, 3, 3, 6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!           fputs(fd, "             keep jacobian matrix,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 300,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 300;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           fputs(fd, "        threads: assembly, 1;\n");
%!           fputs(fd, "    eigenanalysis: list, 1, t1,\n");
%!           fputs(fd, "    output matrices, \n");
%!           fprintf(fd, "          parameter, %.16e,\n", 1);
%!           fputs(fd, "    output eigenvectors,\n");
%!           fputs(fd, "        output geometry,\n");
%!           fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1 / SI_unit_second^-1, 100000 / SI_unit_second^-1);
%!           switch (l)
%!           case 1
%!             fputs(fd, "    use lapack, balance, permute, suffix format, \"%02d\";\n");
%!           case 2
%!             fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", 2 * cms_opt.modes.number, 4 * cms_opt.modes.number + 1);
%!           endswitch
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "        output meter: closest next, 0., forever, t1 / 100.;\n");
%!           switch (l)
%!           case 2
%!             fputs(fd, "        rigid body kinematics: drive,\n");
%!             fputs(fd, "            angular velocity,\n");
%!             fputs(fd, "                   component,\n");
%!             for i=1:3
%!               fprintf(fd, "                postponed, drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             endfor
%!             fputs(fd, "            acceleration,\n");
%!             fputs(fd, "                   component,\n");
%!             for i=1:3
%!               fprintf(fd, "               postponed, drive_id_XPP%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!             fputs(fd, "            angular acceleration,\n");
%!             fputs(fd, "                   component");
%!             for i=1:3
%!               fprintf(fd, ",\n               postponed, drive_id_OMEGAP%s", {"x","y","z"}{i});
%!             endfor
%!             fputs(fd, ";\n");
%!           endswitch
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: none, structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural nodes:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# interface 1\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        joints:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# ground\n");
%!               fputs(fd, "        ;\n");
%!             case 2
%!               fputs(fd, "       structural nodes: 2 * M + 1;\n");
%!               fputs(fd, "       rigid bodies: 2 * M + 1;\n");
%!               fputs(fd, "       beams: M;\n");
%!               fputs(fd, "       joints: 1;\n");
%!           endswitch
%!           fputs(fd, "        forces: 2;\n");
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "drive caller: drive_id_PHI1, string, \"(((pi*Time)/(2*t1)-sin((pi*Time)/t1)/2)*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA1, string, \"sin((pi*Time)/(2*t1))^2\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP1, string, \"(pi*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1)))/t1\";\n");
%!           fputs(fd, "drive caller: drive_id_PHI2, string, \"-(4*sin((pi*Time)/(2*t1))^3*t1^2)/(3*pi^2)\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA2, string, \"-(2*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1))^2*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP2, string, \"-(3*cos((pi*Time)/(2*t1))^2-1)*sin((pi*Time)/(2*t1))\";\n");
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_PHI%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_PHI1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGA1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGAP%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGAP1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_X%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XPP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_g%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, g%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_F1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, F1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_M1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, M1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_tilt,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, euler123, gammax, gammay, gammaz,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "reference: ref_id_modal,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fputs(fd, "reference: ref_id_interface1,\n");
%!               fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!               fputs(fd, "        reference, ref_id_modal, eye,\n");
%!               fputs(fd, "        reference, ref_id_modal, null,\n");
%!               fputs(fd, "        reference, ref_id_modal, null;\n");
%!             case 2
%!               fputs(fd, "reference: ref_id_beam1,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!           endswitch
%!           fputs(fd, "begin: nodes;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural: node_id_modal, modal,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, eye,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!                 fprintf(fd, "                reference, ref_id_beam1, 0.5 * a / M * %d, 0., 0.,\n", i - 1);
%!                 fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null, accelerations, yes;\n");
%!               endfor
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           switch (l)
%!           case 1
%!           fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!           fprintf(fd, "                %s,\n", {"node_id_modal", "node_id_beam1"}{l});
%!           fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fputs(fd, "               position constraint,\n");
%!           if (~param(j, k).holonomic)
%!             fputs(fd, "                        velocity, velocity, velocity,\n");
%!             fputs(fd, "                        component,\n");
%!             for i=1:3
%!               fprintf(fd, "                      reference, drive_id_XP%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!           else
%!             fputs(fd, "                        active, active, active,\n");
%!             fputs(fd, "                        component,\n");
%!             for i=1:3
%!               fprintf(fd, "                      reference, drive_id_X%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!           endif
%!           fputs(fd, "               orientation constraint,\n");
%!           if (~param(j, k).holonomic)
%!             fputs(fd, "                        angular velocity, angular velocity, angular velocity,\n");
%!             fputs(fd, "                        component");
%!             for i=1:3
%!               fprintf(fd, ",\n                   reference, drive_id_OMEGA%s", {"x","y","z"}{i});
%!             endfor
%!           else
%!             fputs(fd, "                        active, active, active,\n");
%!             fputs(fd, "                        component");
%!             for i=1:3
%!               fprintf(fd, ",\n                   reference, drive_id_PHI%s", {"x","y","z"}{i});
%!             endfor
%!           endif
%!           fputs(fd, ";\n");
%!           case 2
%!             fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!           endswitch
%!           switch (l)
%!             case 1
%!               fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!                 fputs(fd, "               rho * A * a / (2 * M + 1),\n");
%!                 fputs(fd, "               reference, node, null, \n");
%!                 fputs(fd, "               diag,   rho * Ip * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iy * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iz * a / (2 * M + 1),\n");
%!                 fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!               endfor
%!               for i=1:options.number_of_beams
%!                 fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "               linear elastic generic, \n");
%!                 fputs(fd, "               diag, E * A , G * As, G * As, \n");
%!                 fputs(fd, "                     G * It, E * Iy, E * Iz,\n");
%!                 fputs(fd, "               same,\n");
%!                 fputs(fd, "               same;\n");
%!               endfor
%!           endswitch
%!           fprintf(fd, "        force: force_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M"}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_F1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fprintf(fd, "        couple: torque_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M"}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_M1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fputs(fd, "        gravity: uniform, component");
%!           for i=1:3
%!             fprintf(fd, ",\n       reference, drive_id_g%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd,";\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         options_mbd.output_file = sprintf("%s_%d_%d_%d_mbd", filename, j, k, l);
%!         if (~options.verbose)
%!           options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!         endif
%!         options_mbd.mbdyn_command = "mbdyn";
%!         options_eig.positive_frequencies = false;
%!         if (options.verbose)
%!           shell(sprintf("cat %s | nl", filename_mbdyn));
%!         endif
%!         res(j, k, l).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!         output_file_rel_frame = [options_mbd.output_file, "_rel"];
%!         mbdyn_post_abs_to_rel(1, options_mbd.output_file, output_file_rel_frame, 0);
%!         exts = {".log", ".out"};
%!         for i=1:numel(exts)
%!           [err, msg] = symlink([options_mbd.output_file, exts{i}], [output_file_rel_frame, exts{i}]);
%!           if (err ~= 0)
%!             error("failed to create symlink: %s", msg);
%!           endif
%!         endfor
%!         [res(j, k, l).t, ...
%!          res(j, k, l).trajectory, ...
%!          res(j, k, l).deformation, ...
%!          res(j, k, l).velocity, ...
%!          res(j, k, l).acceleration, ...
%!          res(j, k, l).node_id, ...
%!          res(j, k, l).force, ...
%!          res(j, k, l).force_id, ...
%!          res(j, k, l).force_node_id, ...
%!          res(j, k, l).orientation_description] = mbdyn_post_load_output_struct(output_file_rel_frame);
%!         res(j, k, l).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!         [res(j, k, l).drive_id, ...
%!         res(j, k, l).drive_value] = mbdyn_post_load_output_drv(options_mbd.output_file, [], numel(res(j, k, l).t));
%!         res(j, k, l).modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, 0);
%!       endfor
%!     endfor
%!   endfor
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = options.number_of_threads;
%!   load_case = struct("omega", empty_cell, ...
%!                      "omegadot", empty_cell, ...
%!                      "loads", empty_cell, ...
%!                      "loaded_nodes", empty_cell, ...
%!                      "joints", empty_cell, ...
%!                      "g", empty_cell, ...
%!                      "tau0", empty_cell);
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   for i=1:numel(load_case)
%!     load_case(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!     load_case(i).loads = zeros(1, 6);
%!     load_case(i).omega = zeros(3, 1);
%!     load_case(i).omegadot = zeros(3, 1);
%!     load_case(i).g = zeros(3, 1);
%!     load_case(i).tau0.iso20 = zeros(rows(mesh.elements.iso20), columns(mesh.elements.iso20), 6);
%!   endfor
%!   sol_eig = struct("def", empty_cell, "lambda", empty_cell, "f", empty_cell);
%!   sol_eig_red = struct("lambda_red", empty_cell, "Ured", empty_cell);
%!   for j=idx_j
%!     for k=idx_k
%!       R = euler123_to_rotation_matrix(param(j, k).gamma);
%!       load_case(j, k).omega = R.' * param(j, k).omega;
%!       load_case(j, k).omegadot = R.' * param(j, k).omegadot;
%!       load_case(j, k).loads = [(R.' * param(j, k).F1).', (R.' * param(j, k).M1).'];
%!       load_case(j, k).g = R.' * (param(j, k).g - param(j, k).a);
%!       [mat_ass.M, ...
%!        mat_ass.D, ...
%!        mat_ass.K, ...
%!        mat_ass.KOMEGA, ...
%!        mat_ass.KOMEGA_DOT, ...
%!        mat_ass.DOMEGA, ...
%!        mat_ass.R, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                             FEM_MAT_DAMPING_OMEGA, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case(j, k));
%!       cms_opt.symmetric = false;
%!       sol_statjk = fem_sol_static(mesh, dof_map, mat_ass, cms_opt);
%!       sol_statjk.stress = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_VEC_STRESS_CAUCH], ...
%!                                          load_case(j, k), ...
%!                                          sol_statjk);
%!       sol_stat(j, k) = sol_statjk;
%!       load_case(j, k).tau0 = sol_stat(j, k).stress.tau;
%!       mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case(j, k));
%!       mat_ass.K += mat_ass.KOMEGA + mat_ass.KOMEGA_DOT + mat_ass.KTAU0;
%!       mat_ass.D += mat_ass.DOMEGA;
%!       sol_eig(j, k) = fem_sol_modal_damped(mesh, ...
%!                                            dof_map, ...
%!                                            mat_ass, ...
%!                                            cms_opt.modes.number, ...
%!                                            cms_opt);
%!       Mred = mat_ass_cms.Mred;
%!       Dred = mat_ass_cms.Dred;
%!       Kred = mat_ass_cms.Kred;
%!       Dred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.DOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA_DOT(dof_map.idx_node, dof_map.idx_node), "Full");
%!       omegaq = [load_case(j, k).omega.^2;
%!                 load_case(j, k).omega(1) * load_case(j, k).omega(2);
%!                 load_case(j, k).omega(2) * load_case(j, k).omega(3);
%!                 load_case(j, k).omega(3) * load_case(j, k).omega(1)];
%!       idx = int32(0);
%!       for i=1:numel(omegaq)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * omegaq(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).omegadot)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).omegadot(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).g)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred -= mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).g(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).loads)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).loads(i);
%!       endfor
%!       [sol_eig_red(j, k).Ured, sol_eig_red(j, k).lambda_red] = fem_sol_eigsd(Kred, Dred, Mred, cms_opt.modes.number, cms_opt);
%!     endfor
%!   endfor
%! tol_abs = [0, 0] / SI_unit_second^-1;
%! tol_rel = [0.3e-2, 3e-2];
%! tol_disp_rel = 3e-2;
%! err_u_modal = err_v_modal = zeros(size(param));
%! printf("deformation/velocity:\n");
%! colors = rainbow(3);
%! width = 1:size(res, 3);
%! linestyle = {"-", "--"};
%! for i=idx_j
%!   for j=idx_k
%!      u_modal = res(i, j, 1).trajectory{end} - res(i, j, 1).trajectory{end}(1, :);
%!      u_beam = res(i, j, 2).trajectory{end} - res(i, j, 2).trajectory{end}(1, :);
%!      v_modal = res(i, j, 1).velocity{end};
%!      v_beam = res(i, j, 2).velocity{end};
%!      if (options.plot)
%!      figure("visible", "off");
%!      hold on;
%!      for k=1:size(res, 3)
%!        for l=1:3
%!          hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l) - res(i, j, k).trajectory{end}(1, l)) * SI_unit_meter);
%!          set(hnd, "color", colors(l, :));
%!          set(hnd, "linewidth", width(k));
%!          set(hnd, "linestyle", linestyle{k});
%!        endfor
%!      endfor
%!      xlabel("t [s]");
%!      ylabel("u [m]");
%!      grid on;
%!      grid minor on;
%!      title(sprintf("linear displacement %d:%d", i, j));
%!      figure("visible", "off");
%!      hold on;
%!      for k=1:size(res, 3)
%!        for l=1:3
%!          hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l + 3) - res(i, j, k).trajectory{end}(1, l + 3)) * 180 / pi);
%!          set(hnd, "color", colors(l, :));
%!          set(hnd, "linewidth", width(k));
%!          set(hnd, "linestyle", linestyle{k});
%!        endfor
%!      endfor
%!      xlabel("t [s]");
%!      ylabel("Phi [deg]");
%!      grid on;
%!      grid minor on;
%!      title(sprintf("angular displacement %d:%d", i, j));
%!      endif
%!      err_u_modal(i, j) = max(max(abs(u_modal - u_beam))) / max(1, max(max(abs(u_beam))));
%!      err_v_modal(i, j) = max(max(abs(v_modal - v_beam))) / max(1, max(max(abs(v_beam))));
%!      printf("%d:%d %.1f%%/%.1f%%\n", i, j, 100 * err_u_modal(i, j), 100 * err_v_modal(i, j));
%!   endfor
%! endfor
%! printf("natural frequencies:\n");
%! MACR = cell(size(param));
%! result_data = struct("f_mbd", cell(size(param)), "f_fem", cell(size(param)));
%! for i=idx_j
%!   for j=idx_k
%!     f_fem = sort(sol_eig(i, j).f(:));
%!     f_fem = f_fem(f_fem > 0);
%!     f_mbd = zeros(rows(f_fem), size(res, 3));
%!     PhiR = zeros(6, rows(f_fem), size(res, 3));
%!     for k=1:size(res, 3)
%!       [f_mbd_k, idx_mbd_k] = sort(res(i, j, k).modal.f(:));
%!       D_mbd_k = res(i, j, k).modal.D(idx_mbd_k);
%!       idx_mbd_k = idx_mbd_k(f_mbd_k > 0);
%!       f_mbd_k = f_mbd_k(f_mbd_k > 0);
%!       idx_mbd_k = idx_mbd_k(1:rows(f_fem));
%!       f_mbd(:, k) = f_mbd_k(1:rows(f_fem));
%!       PhiR(:, :, k) = res(i, j, k).modal.VR(res(i, j, k).modal.idx(end) + (1:6), idx_mbd_k);
%!     endfor
%!     result_data(i, j).f_fem = f_fem;
%!     result_data(i, j).f_mbd = f_mbd;
%!     MACR{i, j} = MACL{i, j} = zeros(rows(f_fem), rows(f_fem));
%!     for k=1:rows(f_fem)
%!       for l=1:rows(f_fem)
%!         MACR{i, j}(k, l) = (PhiR(:, k, 1)' * PhiR(:, k, 2)) * conj(PhiR(:, k, 1)' * PhiR(:, k, 2)) / ((PhiR(:, k, 1)' * PhiR(:, k, 1)) * (PhiR(:, k, 2)' * PhiR(:, k, 2)));
%!       endfor
%!     endfor
%!     printf("%d:%d\n", i, j);
%!     for k=1:rows(f_fem)
%!       printf("%10.2f", f_fem(k) * SI_unit_second^-1);
%!       for l=1:columns(f_mbd)
%!         printf("\t%10.2f", f_mbd(k, l) * SI_unit_second^-1);
%!       endfor
%!       for l=1:columns(f_mbd)
%!         printf("\t%.1f%%", 100 * (f_mbd(k, l) / f_fem(k) - 1));
%!       endfor
%!       printf("\t%.3f", MACR{i, j}(k, k));
%!       fputs(stdout, "\n");
%!     endfor
%!    fputs(stdout, "\n\n");
%!   endfor
%! endfor
%! for i=idx_j
%!   for j=idx_k
%!     for k=1:rows(result_data(i, j).f_fem)
%!       for l=1:columns(result_data(i, j).f_mbd)
%!         assert(result_data(i, j).f_mbd(k, l), result_data(i, j).f_fem(k), tol_abs(l) + tol_rel(l) * abs(result_data(i, j).f_fem(k)));
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! assert(all(all(err_u_modal < tol_disp_rel)));
%! for j=idx_j
%!   for k=idx_k
%!       tol = 2e-2;
%!       [lambda_s] = sortrows([imag(sol_eig(j, k).lambda), real(sol_eig(j, k).lambda)],[1,2]);
%!       [lambda_red_s] = sortrows([imag(sol_eig_red(j, k).lambda_red), real(sol_eig_red(j, k).lambda_red)],[1,2]);
%!       K = min(20, rows(lambda_s));
%!       lambda_s = 1j * lambda_s(:,1) + lambda_s(:, 2);
%!       lambda_red_s = 1j * lambda_red_s(:, 1) + lambda_red_s(:, 2);
%!       assert(lambda_red_s(1:K), lambda_s(1:K), tol * norm(lambda_s(1:K)));
%!   endfor
%! endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST2
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 5;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round((l1 - 0.5 * w2) / hw1)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round((l2 - 0.5 * w1) / hw2)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(w2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round((l2 - 0.5 * w1) / hw2)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Max(1, Round(w1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Max(1, Round(w2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Max(1, Round((l1 - 0.5 * w2) / hw1)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Max(1, Round(w1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Max(1, Round(w1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Max(1, Round(w2 / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Max(1, Round(h1 / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{2,12};\n");
%!     fputs(fd, "Recombine Surface{3,16};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "1", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.iso8].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.iso8].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.iso4].id] == 1);
%!   mesh.materials.iso8(mesh.groups.iso8(grp_idx_beam1).elements) = 1;
%!   mesh.materials.iso8(mesh.groups.iso8(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.iso4(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = int32(4);
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = int32(4);
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.03);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST3
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 3;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2, param.h1]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Round(h1 / h)}; Recombine;};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.penta15].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.penta15].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_beam1).elements) = 1;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = int32(4);
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = int32(4);
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.025);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST4
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 3;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Round(h1 / h)}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{2,12};\n");
%!     fputs(fd, "Recombine Surface{3,16};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.iso20].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam1).elements) = 1;
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = int32(4);
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = int32(4);
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.04);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST5
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 2;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3};};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1, 2, 3}; } } = h;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMin=0.9;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMax=1.1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh);
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.tet10h].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.tet10h].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.tria6h].id] == 1);
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_beam1).elements) = 1;
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = int32(4);
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = int32(4);
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                              dof_map, ...
%!                              mat_ass, ...
%!                              options.number_of_modes, ...
%!                              0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST6
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 2;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3};};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1, 2, 3}; } } = h;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=1;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMin=0.9;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMax=1.1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "3", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   opt_mesh.elem_type = {"tet20", "tria10"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh);
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.tet20 = zeros(rows(mesh.elements.tet20), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.tet20].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.tet20].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.tria10].id] == 1);
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_beam1).elements) = 1;
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.tria10(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = int32(4);
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = int32(4);
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert(max(abs(sol_eig_diag.f / sol_eig(1).f - 1)) < 0.03);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST7
%! printf("fem_cms_create2: test7\n");
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   a = 800e-3 / SI_unit_meter;
%!   b = 40e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   h = 2 * c;
%!   options.interactive = false;
%!   options.plot = true;
%!   options.verbose = false;
%!   options.number_of_beams = int32(40);
%!   options.number_of_threads = int32(4);
%!   if (options.plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fprintf(fd, "a = %.16e;\n", a);
%!     fprintf(fd, "b = %.16e;\n", b);
%!     fprintf(fd, "c = %.16e;\n", c);
%!     fprintf(fd, "h = %.16e;\n", h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(3) = {a,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(4) = {a, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,c}{Surface{1}; Layers{Max(1, Round(c / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {v1[2]};\n");
%!     fputs(fd, "Physical Surface(2) = {v1[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   mesh.material_data.E = 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data.alpha = 0e-5 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0e-5 / (SI_unit_second);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "diag-shift-invert";
%!   cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!   cms_opt.nodes.modal.name = "node_id_modal";
%!   cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!   cms_opt.nodes.interfaces.name = "node_id_interface1";
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [0, 0, 0];
%!   mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [2], [cms_opt.nodes.interfaces.number], "quad8");
%!   cms_opt.refine_max_iter = 30;
%!   cms_opt.pre_scaling = false;
%!   cms_opt.solver = "pardiso";
%!   cms_opt.modes.number = 20;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, 1:3) = true;
%!   [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms, cms_opt, sol_tau_cms] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%!   pert.omega = [1e2; 3e2; 2e2] / (SI_unit_rad / SI_unit_second);
%!   pert.omegadot = [1e5; 1e3; 3e3] / (SI_unit_rad / SI_unit_second^2);
%!   pert.loads = [[1e4; 1e3; 1e2] / (SI_unit_newton);
%!                 [1e2; 1e1; 1e1] / (SI_unit_newton * SI_unit_meter)];
%!   pert.g = [1e4; -1e3; -1e2] / (SI_unit_meter / SI_unit_second^2);
%!   pert.a = [-1e4; 1e3; 1e2] / (SI_unit_meter / SI_unit_second^2);
%!   empty_cell = cell(7, 3, 2);
%!   res = struct("info", empty_cell, ...
%!                "t", empty_cell, ...
%!                "trajectory", empty_cell, ...
%!                "deformation", empty_cell, ...
%!                "velocity", empty_cell, ...
%!                "acceleration", empty_cell, ...
%!                "node_id", empty_cell, ...
%!                "force", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "orientation_description", empty_cell, ...
%!                "drive_id", empty_cell, ...
%!                "drive_value", empty_cell, ...
%!                "modal", empty_cell);
%!   empty_cell = cell(7, 3);
%!   param = struct("omega", empty_cell, ...
%!                  "omegadot", empty_cell, ...
%!                  "F1", empty_cell, ...
%!                  "M1", empty_cell, ...
%!                  "a", empty_cell, ...
%!                  "g", empty_cell, ...
%!                  "t1", empty_cell, ...
%!                  "holonomic", empty_cell);
%!   idx_j = 1:rows(param);
%!   idx_k = 1:columns(param);
%!   for j=idx_j
%!     for k=idx_k
%!       param(j, k).omega = zeros(3, 1);
%!       param(j, k).omegadot = zeros(3, 1);
%!       param(j, k).F1 = zeros(3, 1);
%!       param(j, k).M1 = zeros(3, 1);
%!       param(j, k).a = zeros(3, 1);
%!       param(j, k).g = zeros(3, 1);
%!       param(j, k).holonomic = false;
%!       param(j, k).gamma = zeros(3, 1);
%!       param(j, k).N = 50;
%!       switch (j)
%!         case 1
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).gamma = [20; 45; 30] * pi / 180;
%!         case 2
%!           param(j, k).omega(k) = pert.omega(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           param(j, k).holonomic = true;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 2000;
%!         case 3
%!           param(j, k).omegadot(k) = pert.omegadot(k);
%!           param(j, k).t1 = 1e-2 / SI_unit_second;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 200;
%!         case 4
%!           param(j, k).F1(k) = pert.loads(k);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 5
%!           param(j, k).M1(k) = pert.loads(k + 3);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           param(j, k).gamma(1) = 30 * pi / 180;
%!         case 6
%!           param(j, k).a(k) = pert.a(k);
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).N = 200;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 7
%!           param(j, k).g(k) = pert.g(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!       endswitch
%!       for l=1:3
%!         switch (l)
%!           case 3
%!             nodes_file = sprintf("%s_%d_%d_%d.nod", filename, j, k, l);
%!             csl_file = sprintf("%s_%d_%d_%d.csl", filename, j, k, l);
%!             elem_file = sprintf("%s_%d_%d_%d.elm", filename, j, k, l);
%!             opt_mbd_mesh = struct();
%!             opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_clamp";
%!             opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!             opt_mbd_mesh.struct_nodes.type(cms_opt.nodes.interfaces.number) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!             opt_mbd_mesh.struct_nodes.type(cms_opt.nodes.modal.number) = MBDYN_NODE_TYPE_STATIC_STRUCT_DISP;
%!             opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!             opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!             load_case_empty = struct();
%!             opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!         endswitch
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%d_%d_%d.mbdyn", filename, j, k, l);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fprintf(fd, "set: real a = %.16e;\n", a);
%!           fprintf(fd, "set: real b = %.16e;\n", b);
%!           fprintf(fd, "set: real c = %.16e;\n", c);
%!           fprintf(fd, "set: real d = %.16e;\n", d);
%!           for i=1:3
%!             fprintf(fd, "set: real gamma%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).gamma(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGA%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omega(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGAP%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omegadot(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real F1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).F1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real M1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).M1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real a%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).a(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).g(i));
%!           endfor
%!           fprintf(fd, "set: real t1 = %.16e;\n", param(j, k).t1);
%!           fprintf(fd, "set: integer N = %d;\n", param(j, k).N);
%!           fprintf(fd, "set: integer M = %d;\n", options.number_of_beams);
%!           fputs(fd, "set: integer ref_id_ground = 1;\n");
%!           fputs(fd, "set: integer ref_id_tilt = 2;\n");
%!           fputs(fd, "set: integer joint_id_ground = 1;\n");
%!           fputs(fd, "set: integer force_id1;\n");
%!           fputs(fd, "set: integer torque_id1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_PHI1 = 1;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP1 = 3;\n");
%!           fputs(fd, "set: integer drive_id_PHI2 = 4;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA2 = 5;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP2 = 6;\n");
%!           fputs(fd, "set: integer drive_id_PHIx = 7;\n");
%!           fputs(fd, "set: integer drive_id_PHIy = 8;\n");
%!           fputs(fd, "set: integer drive_id_PHIz = 9;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAx = 10;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAy = 11;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAz = 12;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPx = 13;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPy = 14;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPz = 15;\n");
%!           fputs(fd, "set: integer drive_id_Xx = 16;\n");
%!           fputs(fd, "set: integer drive_id_Xy = 17;\n");
%!           fputs(fd, "set: integer drive_id_Xz = 18;\n");
%!           fputs(fd, "set: integer drive_id_XPx = 19;\n");
%!           fputs(fd, "set: integer drive_id_XPy = 20;\n");
%!           fputs(fd, "set: integer drive_id_XPz = 21;\n");
%!           fputs(fd, "set: integer drive_id_XPPx = 22;\n");
%!           fputs(fd, "set: integer drive_id_XPPy = 23;\n");
%!           fputs(fd, "set: integer drive_id_XPPz = 24;\n");
%!           fputs(fd, "set: integer drive_id_gx = 25;\n");
%!           fputs(fd, "set: integer drive_id_gy = 26;\n");
%!           fputs(fd, "set: integer drive_id_gz = 27;\n");
%!           fputs(fd, "set: integer drive_id_F1x = 28;\n");
%!           fputs(fd, "set: integer drive_id_F1y = 29;\n");
%!           fputs(fd, "set: integer drive_id_F1z = 30;\n");
%!           fputs(fd, "set: integer drive_id_M1x = 31;\n");
%!           fputs(fd, "set: integer drive_id_M1y = 32;\n");
%!           fputs(fd, "set: integer drive_id_M1z = 33;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "set: integer ref_id_modal = 3;\n");
%!               fputs(fd, "set: integer node_id_modal = 1;\n");
%!               fputs(fd, "set: integer ref_id_interface1 = 4;\n");
%!               fputs(fd, "set: integer node_id_interface1 = 2;\n");
%!               fputs(fd, "set: integer elem_id_modal = 2;\n");
%!             case 2
%!               fputs(fd, "set: integer ref_id_beam1 = 3;\n");
%!               fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!               fprintf(fd, "set: real E = %.16e;\n", mesh.material_data.E);
%!               fprintf(fd, "set: real nu = %.16e;\n", mesh.material_data.nu);
%!               fprintf(fd, "set: real rho = %.16e;\n", mesh.material_data.rho);
%!               fprintf(fd, "set: real alpha = %.16e;\n", mesh.material_data.alpha);
%!               fprintf(fd, "set: real beta = %.16e;\n", mesh.material_data.beta);
%!               fputs(fd, "set: real G = E / (2. * (1. + nu));\n");
%!               fputs(fd, "set: real A = b * c;\n");
%!               fputs(fd, "set: real As = 9. / 10. * A;\n");
%!               fputs(fd, "set: real Iy = b * c^3 / 12.;\n");
%!               fputs(fd, "set: real Iz = c * b^3 / 12.;\n");
%!               fputs(fd, "set: real Ip = Iy + Iz;\n");
%!               c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!               w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!               fprintf(fd, "set: real It = %.16e;\n", interp1(w_h, c2, max(c, b) / min(c, b)) * max(c, b) * min(c, b)^3);
%!             case 3
%!               fputs(fd, "set: integer ref_id_clamp = 3;\n");
%!               fprintf(fd, "set: integer node_id_interface1 = %d;\n", cms_opt.nodes.interfaces.number);
%!           endswitch
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / N;\n");
%!           fputs(fd, "        max time step: t1 / N;\n");
%!           fputs(fd, "        min time step: t1 / N;\n");
%!           fputs(fd, "        method: ss4, 0.;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1;\n");
%!           fputs(fd, "        strategy: factor, 0.8, 3, 1.25, 3, 3, 6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes, cpu time;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!           fputs(fd, "             keep jacobian matrix,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 30,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 30;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           switch (l)
%!           case 3
%!             fprintf(fd, "        threads: assembly, %d;\n", cms_opt.number_of_threads);
%!           otherwise
%!             fputs(fd, "        threads: assembly, 1;\n");
%!           endswitch
%!           fputs(fd, "    eigenanalysis: list, 1, t1,\n");
%!           fputs(fd, "    output eigenvectors,\n");
%!           fputs(fd, "    # output matrices, \n");
%!           fputs(fd, "        output geometry,\n");
%!           fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1. / SI_unit_second^-1, 10000 / SI_unit_second^-1);
%!           switch (l)
%!             case 1
%!               fputs(fd, "    use lapack, balance, permute, suffix format, \"%02d\";\n");
%!             case {2, 3}
%!               fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", cms_opt.modes.number, 2 * cms_opt.modes.number + 1);
%!           endswitch
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "        output meter: closest next, 0., forever, t1 / 20.;\n");
%!           switch (l)
%!             case {2, 3}
%!               fputs(fd, "        rigid body kinematics: drive,\n");
%!               fputs(fd, "            angular velocity,\n");
%!               fputs(fd, "                   component,\n");
%!               for i=1:3
%!                 fprintf(fd, "                postponed, drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!               endfor
%!               fputs(fd, "            acceleration,\n");
%!               fputs(fd, "                   component,\n");
%!               for i=1:3
%!                 fprintf(fd, "               postponed, drive_id_XPP%s,\n", {"x", "y", "z"}{i});
%!               endfor
%!               fputs(fd, "            angular acceleration,\n");
%!               fputs(fd, "                   component");
%!               for i=1:3
%!                 fprintf(fd, ",\n               postponed, drive_id_OMEGAP%s", {"x","y","z"}{i});
%!               endfor
%!               fputs(fd, ";\n");
%!           endswitch
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: none, structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural nodes:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# interface 1\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        joints:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# ground\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        forces: 2;\n");
%!             case 2
%!               fputs(fd, "       structural nodes: 2 * M + 1;\n");
%!               fputs(fd, "       rigid bodies: 2 * M + 1;\n");
%!               fputs(fd, "       beams: M;\n");
%!               fputs(fd, "       joints: 1;\n");
%!               fputs(fd, "       forces: 2;\n");
%!             case 3
%!               fprintf(fd, "     structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!               fprintf(fd, "     solids: %d;\n", opt_mbd_mesh.solids.number);
%!               fprintf(fd, "     genels: %d;\n", opt_mbd_mesh.genels.number);
%!               fprintf(fd, "     joints: %d;\n", opt_mbd_mesh.joints.number);
%!               fprintf(fd, "     forces: %d;\n", opt_mbd_mesh.forces.number + 2);
%!           endswitch
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "drive caller: drive_id_PHI1, string, \"(((pi*Time)/(2*t1)-sin((pi*Time)/t1)/2)*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA1, string, \"sin((pi*Time)/(2*t1))^2\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP1, string, \"(pi*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1)))/t1\";\n");
%!           fputs(fd, "drive caller: drive_id_PHI2, string, \"-(4*sin((pi*Time)/(2*t1))^3*t1^2)/(3*pi^2)\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA2, string, \"-(2*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1))^2*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP2, string, \"-(3*cos((pi*Time)/(2*t1))^2-1)*sin((pi*Time)/(2*t1))\";\n");
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_PHI%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_PHI1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGA1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGAP%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGAP1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_X%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XPP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_g%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, g%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_F1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, F1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_M1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, M1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_tilt,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, euler123, gammax, gammay, gammaz,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "reference: ref_id_modal,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fputs(fd, "reference: ref_id_interface1,\n");
%!               fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!               fputs(fd, "        reference, ref_id_modal, eye,\n");
%!               fputs(fd, "        reference, ref_id_modal, null,\n");
%!               fputs(fd, "        reference, ref_id_modal, null;\n");
%!             case 2
%!               fputs(fd, "reference: ref_id_beam1,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!             case 3
%!               fputs(fd, "reference: ref_id_clamp,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fprintf(fd, "include: \"%s\";\n", csl_file);
%!           endswitch
%!           fputs(fd, "begin: nodes;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural: node_id_modal, modal,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, eye,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!                 fprintf(fd, "                reference, ref_id_beam1, 0.5 * a / M * %d, 0., 0.,\n", i - 1);
%!                 fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null, accelerations, yes;\n");
%!               endfor
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", nodes_file);
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!               fprintf(fd, "                %s,\n", {"node_id_modal", "node_id_beam1"}{l});
%!               fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fputs(fd, "               position constraint,\n");
%!               if (~param(j, k).holonomic)
%!                 fputs(fd, "                        velocity, velocity, velocity,\n");
%!                 fputs(fd, "                        component,\n");
%!                 for i=1:3
%!                   fprintf(fd, "                      reference, drive_id_XP%s,\n", {"x", "y", "z"}{i});
%!                 endfor
%!               else
%!                 fputs(fd, "                        active, active, active,\n");
%!                 fputs(fd, "                        component,\n");
%!                 for i=1:3
%!                   fprintf(fd, "                      reference, drive_id_X%s,\n", {"x", "y", "z"}{i});
%!                 endfor
%!               endif
%!               fputs(fd, "               orientation constraint,\n");
%!               if (~param(j, k).holonomic)
%!                 fputs(fd, "                        angular velocity, angular velocity, angular velocity,\n");
%!                 fputs(fd, "                        component");
%!                 for i=1:3
%!                   fprintf(fd, ",\n                   reference, drive_id_OMEGA%s", {"x","y","z"}{i});
%!                 endfor
%!               else
%!                 fputs(fd, "                        active, active, active,\n");
%!                 fputs(fd, "                        component");
%!                 for i=1:3
%!                   fprintf(fd, ",\n                   reference, drive_id_PHI%s", {"x","y","z"}{i});
%!                 endfor
%!               endif
%!               fputs(fd, ";\n");
%!             case 2
%!               fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!           endswitch
%!           switch (l)
%!             case 1
%!               fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!                 fputs(fd, "               rho * A * a / (2 * M + 1),\n");
%!                 fputs(fd, "               reference, node, null, \n");
%!                 fputs(fd, "               diag,   rho * Ip * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iy * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iz * a / (2 * M + 1),\n");
%!                 fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!               endfor
%!               for i=1:options.number_of_beams
%!                 fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "               linear elastic generic, \n");
%!                 fputs(fd, "               diag, E * A , G * As, G * As, \n");
%!                 fputs(fd, "                     G * It, E * Iy, E * Iz,\n");
%!                 fputs(fd, "               same,\n");
%!                 fputs(fd, "               same;\n");
%!               endfor
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", elem_file);
%!           endswitch
%!           fprintf(fd, "        force: force_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M", sprintf("%d", cms_opt.nodes.interfaces.number)}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_F1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fprintf(fd, "        couple: torque_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M", sprintf("%d", cms_opt.nodes.interfaces.number)}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_M1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fputs(fd, "        gravity: uniform, component");
%!           for i=1:3
%!             fprintf(fd, ",\n       reference, drive_id_g%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd,";\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         options_mbd.output_file = sprintf("%s_%d_%d_%d_mbd", filename, j, k, l);
%!         if (~options.verbose)
%!           options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!         endif
%!         options_mbd.mbdyn_command = "mbdyn -C";
%!         options_eig.positive_frequencies = false;
%!         if (options.verbose)
%!           shell(sprintf("cat %s | nl", filename_mbdyn));
%!           switch (l)
%!           case 3
%!             shell(sprintf("cat %s | nl", csl_file));
%!             shell(sprintf("cat %s | nl", nodes_file));
%!             shell(sprintf("cat %s | nl", elem_file));
%!           endswitch
%!         endif
%!         res(j, k, l).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!         output_file_rel_frame = [options_mbd.output_file, "_rel"];
%!         mbdyn_post_abs_to_rel(1, options_mbd.output_file, output_file_rel_frame, 0);
%!         exts = {".log", ".out"};
%!         for i=1:numel(exts)
%!           [err, msg] = symlink([options_mbd.output_file, exts{i}], [output_file_rel_frame, exts{i}]);
%!           if (err ~= 0)
%!             error("failed to create symlink: %s", msg);
%!           endif
%!         endfor
%!         [res(j, k, l).t, ...
%!          res(j, k, l).trajectory, ...
%!          res(j, k, l).deformation, ...
%!          res(j, k, l).velocity, ...
%!          res(j, k, l).acceleration, ...
%!          res(j, k, l).node_id, ...
%!          res(j, k, l).force, ...
%!          res(j, k, l).force_id, ...
%!          res(j, k, l).force_node_id, ...
%!          res(j, k, l).orientation_description] = mbdyn_post_load_output_struct(output_file_rel_frame);
%!         res(j, k, l).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!         [res(j, k, l).drive_id, ...
%!          res(j, k, l).drive_value] = mbdyn_post_load_output_drv(options_mbd.output_file, [], numel(res(j, k, l).t));
%!         res(j, k, l).modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, 0);
%!       endfor
%!     endfor
%!   endfor
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, 1:3) = true;
%!   mesh.elements.joints.nodes = cms_opt.nodes.modal.number;
%!   mesh.elements.joints.C = eye(6);
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = options.number_of_threads;
%!   load_case = struct("omega", empty_cell, ...
%!                      "omegadot", empty_cell, ...
%!                      "loads", empty_cell, ...
%!                      "loaded_nodes", empty_cell, ...
%!                      "joints", empty_cell, ...
%!                      "g", empty_cell, ...
%!                      "tau0", empty_cell);
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   for i=1:numel(load_case)
%!     load_case(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!     load_case(i).loads = zeros(1, 6);
%!     load_case(i).omega = zeros(3, 1);
%!     load_case(i).omegadot = zeros(3, 1);
%!     load_case(i).g = zeros(3, 1);
%!     load_case(i).tau0.iso20 = zeros(rows(mesh.elements.iso20), columns(mesh.elements.iso20), 6);
%!   endfor
%!   sol_eig = struct("def", empty_cell, "lambda", empty_cell, "f", empty_cell);
%!   sol_eig_red = struct("lambda_red", empty_cell, "Ured", empty_cell);
%!   for j=idx_j
%!     for k=idx_k
%!       R = euler123_to_rotation_matrix(param(j, k).gamma);
%!       load_case(j, k).omega = R.' * param(j, k).omega;
%!       load_case(j, k).omegadot = R.' * param(j, k).omegadot;
%!       load_case(j, k).loads = [(R.' * param(j, k).F1).', (R.' * param(j, k).M1).'];
%!       load_case(j, k).g = R.' * (param(j, k).g - param(j, k).a);
%!       [mat_ass.M, ...
%!        mat_ass.D, ...
%!        mat_ass.K, ...
%!        mat_ass.KOMEGA, ...
%!        mat_ass.KOMEGA_DOT, ...
%!        mat_ass.DOMEGA, ...
%!        mat_ass.R, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                             FEM_MAT_DAMPING_OMEGA, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case(j, k));
%!       cms_opt.symmetric = false;
%!       sol_statjk = fem_sol_static(mesh, dof_map, mat_ass, cms_opt);
%!       sol_statjk.stress = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_VEC_STRESS_CAUCH], ...
%!                                          load_case(j, k), ...
%!                                          sol_statjk);
%!       sol_stat(j, k) = sol_statjk;
%!       load_case(j, k).tau0 = sol_stat(j, k).stress.tau;
%!       mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case(j, k));
%!       mat_ass.K += mat_ass.KOMEGA + mat_ass.KOMEGA_DOT + mat_ass.KTAU0;
%!       mat_ass.D += mat_ass.DOMEGA;
%!       sol_eig(j, k) = fem_sol_modal_damped(mesh, ...
%!                                            dof_map, ...
%!                                            mat_ass, ...
%!                                            cms_opt.modes.number, ...
%!                                            cms_opt);
%!       Mred = mat_ass_cms.Mred;
%!       Dred = mat_ass_cms.Dred;
%!       Kred = mat_ass_cms.Kred;
%!       Dred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.DOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA_DOT(dof_map.idx_node, dof_map.idx_node), "Full");
%!       omegaq = [load_case(j, k).omega.^2;
%!                 load_case(j, k).omega(1) * load_case(j, k).omega(2);
%!                 load_case(j, k).omega(2) * load_case(j, k).omega(3);
%!                 load_case(j, k).omega(3) * load_case(j, k).omega(1)];
%!       idx = int32(0);
%!       for i=1:numel(omegaq)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * omegaq(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).omegadot)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).omegadot(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).g)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred -= mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).g(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).loads)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).loads(i);
%!       endfor
%!       [sol_eig_red(j, k).Ured, sol_eig_red(j, k).lambda_red] = fem_sol_eigsd(Kred, Dred, Mred, cms_opt.modes.number, cms_opt);
%!     endfor
%!   endfor
%!   tol_abs = [0, 0, 0] / SI_unit_second^-1;
%!   tol_rel = [0.3e-2, 3e-2, 3e-2];
%!   tol_disp_rel = 3e-2;
%!   err_u_modal = err_v_modal = zeros(size(param));
%!   printf("deformation/velocity:\n");
%!   colors = rainbow(3);
%!   width = 1:size(res, 3);
%!   linestyle = {"-", "--", "-."};
%!   for i=idx_j
%!     for j=idx_k
%!       u_modal = res(i, j, 1).trajectory{end} - res(i, j, 1).trajectory{end}(1, :);
%!       u_beam = res(i, j, 2).trajectory{end} - res(i, j, 2).trajectory{end}(1, :);
%!       v_modal = res(i, j, 1).velocity{end};
%!       v_beam = res(i, j, 2).velocity{end};
%!       if (options.plot)
%!         figure("visible", "off");
%!         hold on;
%!         for k=1:size(res, 3)
%!           for l=1:3
%!             hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l) - res(i, j, k).trajectory{end}(1, l)) * SI_unit_meter);
%!             set(hnd, "color", colors(l, :));
%!             set(hnd, "linewidth", width(k));
%!             set(hnd, "linestyle", linestyle{k});
%!           endfor
%!         endfor
%!         xlabel("t [s]");
%!         ylabel("u [m]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("linear displacement %d:%d", i, j));
%!         figure("visible", "off");
%!         hold on;
%!         for k=1:size(res, 3)
%!           for l=1:3
%!             hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l + 3) - res(i, j, k).trajectory{end}(1, l + 3)) * 180 / pi);
%!             set(hnd, "color", colors(l, :));
%!             set(hnd, "linewidth", width(k));
%!             set(hnd, "linestyle", linestyle{k});
%!           endfor
%!         endfor
%!         xlabel("t [s]");
%!         ylabel("Phi [deg]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("angular displacement %d:%d", i, j));
%!       endif
%!       err_u_modal(i, j) = max(max(abs(u_modal - u_beam))) / max(1, max(max(abs(u_beam))));
%!       err_v_modal(i, j) = max(max(abs(v_modal - v_beam))) / max(1, max(max(abs(v_beam))));
%!       printf("%d:%d %.1f%%/%.1f%%\n", i, j, 100 * err_u_modal(i, j), 100 * err_v_modal(i, j));
%!     endfor
%!   endfor
%!   printf("natural frequencies:\n");
%!   MACR = cell(size(param));
%!   result_data = struct("f_mbd", cell(size(param)), "f_fem", cell(size(param)));
%!   for i=idx_j
%!     for j=idx_k
%!       f_fem = sort(sol_eig(i, j).f(:));
%!       f_fem = f_fem(f_fem > 0);
%!       f_mbd = zeros(rows(f_fem), size(res, 3));
%!       PhiR = zeros(6, rows(f_fem), size(res, 3));
%!       for k=1:size(res, 3)
%!         [f_mbd_k, idx_mbd_k] = sort(res(i, j, k).modal.f(:));
%!         D_mbd_k = res(i, j, k).modal.D(idx_mbd_k);
%!         idx_mbd_k = idx_mbd_k(f_mbd_k > 0);
%!         f_mbd_k = f_mbd_k(f_mbd_k > 0);
%!         idx_mbd_k = idx_mbd_k(1:rows(f_fem));
%!         f_mbd(:, k) = f_mbd_k(1:rows(f_fem));
%!         PhiR(:, :, k) = res(i, j, k).modal.VR(res(i, j, k).modal.idx(end) + (1:6), idx_mbd_k);
%!       endfor
%!       result_data(i, j).f_fem = f_fem;
%!       result_data(i, j).f_mbd = f_mbd;
%!       MACR{i, j} = MACL{i, j} = zeros(rows(f_fem), rows(f_fem));
%!       for k=1:rows(f_fem)
%!         for l=1:rows(f_fem)
%!           MACR{i, j}(k, l) = (PhiR(:, k, 1)' * PhiR(:, k, 2)) * conj(PhiR(:, k, 1)' * PhiR(:, k, 2)) / ((PhiR(:, k, 1)' * PhiR(:, k, 1)) * (PhiR(:, k, 2)' * PhiR(:, k, 2)));
%!         endfor
%!       endfor
%!       printf("%d:%d\n", i, j);
%!       for k=1:rows(f_fem)
%!         printf("%10.2f", f_fem(k) * SI_unit_second^-1);
%!         for l=1:columns(f_mbd)
%!           printf("\t%10.2f", f_mbd(k, l) * SI_unit_second^-1);
%!         endfor
%!         for l=1:columns(f_mbd)
%!           printf("\t%.1f%%", 100 * (f_mbd(k, l) / f_fem(k) - 1));
%!         endfor
%!         printf("\t%.3f", MACR{i, j}(k, k));
%!         fputs(stdout, "\n");
%!       endfor
%!       fputs(stdout, "\n\n");
%!     endfor
%!   endfor
%!   for i=idx_j
%!     for j=idx_k
%!       for k=1:rows(result_data(i, j).f_fem)
%!         for l=1:columns(result_data(i, j).f_mbd)
%!           assert(result_data(i, j).f_mbd(k, l), result_data(i, j).f_fem(k), tol_abs(l) + tol_rel(l) * abs(result_data(i, j).f_fem(k)));
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%!   assert(all(all(err_u_modal < tol_disp_rel)));
%!   for j=idx_j
%!     for k=idx_k
%!       tol = 2e-2;
%!       [lambda_s] = sortrows([imag(sol_eig(j, k).lambda), real(sol_eig(j, k).lambda)],[1,2]);
%!       [lambda_red_s] = sortrows([imag(sol_eig_red(j, k).lambda_red), real(sol_eig_red(j, k).lambda_red)],[1,2]);
%!       K = min(20, rows(lambda_s));
%!       lambda_s = 1j * lambda_s(:,1) + lambda_s(:, 2);
%!       lambda_red_s = 1j * lambda_red_s(:, 1) + lambda_red_s(:, 2);
%!       assert(lambda_red_s(1:K), lambda_s(1:K), tol * norm(lambda_s(1:K)));
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST8
%! printf("fem_cms_create2: test8\n");
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   h = [5e-3; 5e-3; 5e-3];
%!   geometry.zF = 0;
%!   geometry.l = 10e-3;
%!   geometry.w = 10e-3;
%!   geometry.h = 10e-3;
%!   geometry.a = 0;
%!   t1 = 1;
%!   t2 = 0;
%!   dt = t1 / 100;
%!   dto = t1 / 50;
%!   model = "static";
%!   method = "implicit euler";
%!   options.verbose = false;
%!   f_rigid_body_load = false;
%!   f_rigid_body_clamp = false;
%!   material.E = 210000e6;
%!   material.ET = 2100e6;
%!   material.sigmayv = 235000e6;
%!   material.nu = 0.3;
%!   material.rho = 7850;
%!   material.beta = 0;
%!   elem_types = {"iso8", "iso20", "iso20r", "penta15", "tet10h"};
%!   mat_types = {"linear elastic generic", "neo hookean", "bilinear isotropic hardening"};
%!   f_transfinite_mesh = [true, false];
%!   boundary_cond = {"symmetry", "three point"};
%!   load_type = {"traction", "pressure"};
%!   for idx_load_type=1:numel(load_type)
%!   for idx_boundary_cond=1:numel(boundary_cond)
%!     for idx_transfinite=1:numel(f_transfinite_mesh)
%!       for idx_mat_type=1:numel(mat_types)
%!         material.type = mat_types{idx_mat_type};
%!         for idx_elem_type=1:numel(elem_types)
%!           elem_type = elem_types{idx_elem_type};
%!           switch (elem_type)
%!             case "tet10h"
%!               if (f_transfinite_mesh(idx_transfinite))
%!                 continue;
%!               endif
%!           endswitch
%!           file_prefix = sprintf("%s_%d_%d_%d", filename, idx_transfinite, idx_mat_type, idx_elem_type);
%!           geo_file = [file_prefix, "_gmsh.geo"];
%!           mesh_file = [file_prefix, "_gmsh.msh"];
%!           nodes_file = [file_prefix, "_mbd.nod"];
%!           elem_file = [file_prefix, "_mbd.elm"];
%!           set_file = [file_prefix, "_mbd.set"];
%!           csl_file = [file_prefix, "_mbd.csl"];
%!           control_file = [file_prefix, "_mbd.con"];
%!           initial_value_file = [file_prefix, "_mbd.inv"];
%!           input_file = [file_prefix, "_mbd_inp.mbdyn"];
%!           output_file = [file_prefix, "_mbd_out"];
%!           opt_mbd.output_file = output_file;
%!           if (~options.verbose)
%!             opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!           endif
%!           opt_mbd.mbdyn_command = "mbdyn -C";
%!           opt_mbd.f_run_mbdyn = true;
%!           R = eye(3);
%!           F = R * [0.;
%!                    0.;
%!                    0];
%!           px = -material.sigmayv * 1.5;
%!           py = -material.sigmayv * 0.3;
%!           pz = -material.sigmayv * 0.2;
%!           M = R * [0;
%!                    0;
%!                    0];
%!           W = R * [0; 0; 0];
%!           WP = R * [0; 0; 0];
%!           XPP = R * [0; 0; 0];
%!           g = R * [0; 0; 0];
%!           switch (elem_type)
%!             case "iso8"
%!               mesh_order = 1;
%!               elem_type_solid = {elem_type};
%!               elem_type_surf = {"iso4"};
%!             case "iso20"
%!               mesh_order = 2;
%!               elem_type_solid = {"iso20", "penta15"};
%!               elem_type_surf = {"quad8", "tria6h"};
%!             case "iso20r"
%!               mesh_order = 2;
%!               elem_type_solid = {"iso20r", "penta15"}
%!               elem_type_surf = {"quad8r", "tria6h"};
%!             case "penta15"
%!               mesh_order = 2;
%!               elem_type_solid = {elem_type};
%!               elem_type_surf = {"quad8", "tria6h"};
%!             case "tet10h"
%!               mesh_order = 2;
%!               elem_type_solid = {elem_type};
%!               elem_type_surf = {"tria6h"};
%!             otherwise
%!               error("unknown element type \"%s\"", elem_type);
%!           endswitch
%!           fd = -1;
%!           unwind_protect
%!             [fd, msg] = fopen(geo_file, "w");
%!             if (fd == -1)
%!               error("failed to open file \"%s.geo\"", geo_file);
%!             endif
%!             fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!             fprintf(fd, "a=%g;\n", geometry.l);
%!             fprintf(fd, "b=%g;\n", geometry.w);
%!             fprintf(fd, "c=%g;\n", geometry.h);
%!             fprintf(fd, "hx = %g;\n", h(1));
%!             fprintf(fd, "hy = %g;\n", h(2));
%!             fprintf(fd, "hz = %g;\n", h(3));
%!             switch (elem_type)
%!               case {"iso20", "iso20r", "penta15"}
%!                 fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!             endswitch
%!             fprintf(fd, "Mesh.ElementOrder = %d;\n", mesh_order);
%!             fputs(fd, "Point(1) = { 0, -0.5 * b, -0.5 * c};\n");
%!             fputs(fd, "Point(2) = { a, -0.5 * b, -0.5 * c};\n");
%!             fputs(fd, "Point(3) = { a,  0.5 * b, -0.5 * c};\n");
%!             fputs(fd, "Point(4) = { 0,  0.5 * b, -0.5 * c};\n");
%!             fputs(fd, "Line(1) = {4,3};\n");
%!             fputs(fd, "Line(2) = {3,2};\n");
%!             fputs(fd, "Line(3) = {2,1};\n");
%!             fputs(fd, "Line(4) = {1,4};\n");
%!             if (f_transfinite_mesh(idx_transfinite))
%!               fprintf(fd, "Transfinite Curve(1) = Max(2, Round(a / hx) + 1);\n");
%!               fprintf(fd, "Transfinite Curve(2) = Max(2, Round(b / hy) + 1);\n");
%!               fprintf(fd, "Transfinite Curve(3) = Max(2, Round(a / hx) + 1);\n");
%!               fprintf(fd, "Transfinite Curve(4) = Max(2, Round(b / hy) + 1);\n");
%!             endif
%!             fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!             fputs(fd, "Plane Surface(6) = {5};\n");
%!             if (f_transfinite_mesh(idx_transfinite))
%!               fputs(fd, "Transfinite Surface(6) = {2,3,4,1};\n");
%!             endif
%!             fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!             switch (elem_type)
%!               case "tet10h"
%!                 fputs(fd, "  Surface{6};\n");
%!               otherwise
%!                 switch (elem_type)
%!                   case "iso20r"
%!                     num_layers = 2; ## because of hourglass instability
%!                   otherwise
%!                     num_layers = 1;
%!                 endswitch
%!                 fprintf(fd, "  Surface{6}; Layers{Max(%d, Round(c/hz))}; Recombine;\n", num_layers);
%!             endswitch
%!             fputs(fd, "};\n");
%!             f_unstruct_mesh_size = false;
%!             switch (elem_type)
%!               case {"iso8", "iso20", "iso20r"}
%!                 f_unstruct_mesh_size = ~f_transfinite_mesh(idx_transfinite);
%!                 fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!               otherwise
%!                 f_unstruct_mesh_size = true;
%!             endswitch
%!             if (f_unstruct_mesh_size)
%!               fprintf(fd, "MeshSize{PointsOf{Volume{tmp[1]};}} = %.16e;\n", mean(h));
%!             endif
%!             switch (elem_type)
%!               case "tet10h"
%!                 fputs(fd, "Mesh.HighOrderOptimize=1;\n");
%!                 fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!             endswitch
%!             fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!             fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[4]};\n");
%!             fputs(fd, "Physical Surface(\"load+x\",3) = {tmp[2]};\n");
%!             fputs(fd, "Physical Surface(\"load+y\",6) = {tmp[5]};\n");
%!             fputs(fd, "Physical Surface(\"load+z\",7) = {tmp[0]};\n");
%!             fputs(fd, "Physical Surface(\"load-x\",8) = {tmp[4]};\n");
%!             fputs(fd, "Physical Surface(\"load-y\",9) = {tmp[3]};\n");
%!             fputs(fd, "Physical Surface(\"load-z\",10) = {6};\n");
%!             fputs(fd, "Physical Surface(\"symmetry-xy\",4) = {6};\n");
%!             fputs(fd, "Physical Surface(\"symmetry-xz\",5) = {tmp[3]};\n");
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!           end_unwind_protect
%!           pid = spawn("gmsh", {"-format", "msh2", "-3", geo_file});
%!           status = spawn_wait(pid);
%!           if (status ~= 0)
%!             warning("gmsh failed with status %d", status);
%!           endif
%!           opt_msh.elem_type = {elem_type_solid{:}, elem_type_surf{:}};
%!           mesh = fem_pre_mesh_import(mesh_file, "gmsh", opt_msh);
%!           opt_mbd_mesh = struct();
%!           switch (model)
%!             case "dynamic"
%!               opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!             case "static"
%!               opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!           endswitch
%!           grp_idx_volume = zeros(1, numel(elem_type_solid), "int32");
%!           for i=1:numel(elem_type_solid)
%!             if (~isfield(mesh.groups, elem_type_solid{i}))
%!               continue;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_solid{i}).id] == 1);
%!             if (isempty(idx))
%!               continue;
%!             endif
%!             grp_idx_volume(i) = idx;
%!           endfor
%!           grp_idx_load_px = grp_idx_load_py = grp_idx_load_pz = grp_idx_clamp = grp_idx_symmetry_xy = grp_idx_symmetry_xz = zeros(1, numel(elem_type_surf), "int32");
%!           grp_idx_load_mx = grp_idx_load_my = grp_idx_load_mz = zeros(1, numel(elem_type_surf), "int32");
%!           for i=1:numel(elem_type_surf)
%!             if (~isfield(mesh.groups, elem_type_surf{i}))
%!               continue;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 2);
%!             if (~isempty(idx))
%!               grp_idx_clamp(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 3);
%!             if (~isempty(idx))
%!               grp_idx_load_px(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 4);
%!             if (~isempty(idx))
%!               grp_idx_symmetry_xy(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 5);
%!             if (~isempty(idx))
%!               grp_idx_symmetry_xz(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 6);
%!             if (~isempty(idx))
%!               grp_idx_load_py(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 7);
%!             if (~isempty(idx))
%!               grp_idx_load_pz(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 8);
%!             if (~isempty(idx))
%!               grp_idx_load_mx(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 9);
%!             if (~isempty(idx))
%!               grp_idx_load_my(i) = idx;
%!             endif
%!             idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 10);
%!             if (~isempty(idx))
%!               grp_idx_load_mz(i) = idx;
%!             endif
%!           endfor
%!           if (~f_rigid_body_load)
%!             if (norm(F) > 0 || norm(M) > 0)
%!               for i=1:numel(elem_type_surf)
%!                 if (grp_idx_load_px(i) > 0)
%!                   [Fpress, grp_idx_press, weight_press] = fem_pre_mesh_nodal_pressure_load(mesh, 3, elem_type_surf{i});
%!                   mesh_group_i = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_px(i));
%!                   load_case.loaded_nodes = mesh_group_i.nodes(:);
%!                   load_case.loads = zeros(numel(load_case.loaded_nodes), 3);
%!                   for j=1:3
%!                     load_case.loads(:, j) = F(j) * weight_press{1} / sum(weight_press{1});
%!                     load_case.loads(:, j + 3) = M(j) * weight_press{1} / sum(weight_press{1});
%!                   endfor
%!                   break
%!                 endif
%!               endfor
%!             endif
%!           else
%!             load_case.loads = zeros(0, 6);
%!             load_case.loaded_nodes = [];
%!           endif
%!           if (f_rigid_body_load || f_rigid_body_clamp)
%!             inum_elem_rbe3 = 0;
%!             inode_idx_rb = zeros(2, 1);
%!             ielem_idx_rb = zeros(2, 1);
%!           endif
%!           if (f_rigid_body_load)
%!             inode_idx_rb(1) = rows(mesh.nodes) + 1;
%!             ielem_idx_rb(1) = ++inum_elem_rbe3;
%!             opt_mbd_mesh.struct_nodes.type(inode_idx_rb(1)) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!             mesh.nodes(inode_idx_rb(1), 1:3) = [geometry.l + geometry.a, 0, 0];
%!             mesh.elements.rbe3(ielem_idx_rb(1)) = fem_pre_mesh_rbe3_from_surf(mesh, 3, inode_idx_rb(1), elem_type_surf{1});
%!           endif
%!           if (f_rigid_body_clamp)
%!             inode_idx_rb(2) = rows(mesh.nodes) + 1;
%!             ielem_idx_rb(2) = ++inum_elem_rbe3;
%!             opt_mbd_mesh.struct_nodes.type(inode_idx_rb(2)) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!             mesh.nodes(inode_idx_rb(2), 1:3) = [0, 0, 0];
%!             mesh.elements.rbe3(ielem_idx_rb(2)) = fem_pre_mesh_rbe3_from_surf(mesh, 2, inode_idx_rb(2), elem_type_surf{1});
%!           endif
%!           load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!           if (~f_rigid_body_clamp)
%!             for i=1:numel(elem_type_surf)
%!               if (~grp_idx_clamp(i))
%!                 continue;
%!               endif
%!               load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_clamp(i)).nodes, 1:3) = true;
%!             endfor
%!           endif
%!           mesh.nodes(:, 1:3) = mesh.nodes(:, 1:3) * R.';
%!           load_case.g = g - XPP;
%!           load_case.omega = W;
%!           load_case.omegadot = WP;
%!           grp_idx_p1 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!           grp_idx_p2 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!           grp_idx_p3 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!           switch (boundary_cond{idx_boundary_cond})
%!             case "symmetry"
%!               load_case_dof.locked_dof(:, :) = false;
%!               for i=1:numel(elem_type_surf)
%!                 if (grp_idx_clamp(i))
%!                   load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_clamp(i)).nodes, 1) = true;
%!                 endif
%!                 if (grp_idx_symmetry_xz(i))
%!                   load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_symmetry_xz(i)).nodes, 2) = true;
%!                 endif
%!                 if (grp_idx_symmetry_xy(i))
%!                   load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_symmetry_xy(i)).nodes, 3) = true;
%!                 endif
%!               endfor
%!             case "three point"
%!               load_case_dof.locked_dof(:, :) = false;
%!               load_case_dof.locked_dof(grp_idx_p1, 1:3) = true;
%!               load_case_dof.locked_dof(grp_idx_p2, 2:3) = true;
%!               load_case_dof.locked_dof(grp_idx_p3, 3) = true;
%!             otherwise
%!               error("unkown boundary condition");
%!           endswitch
%!           mesh.material_data = material;
%!           mesh.materials = struct();
%!           for i=1:numel(elem_type_solid)
%!             if (~isfield(mesh.elements, elem_type_solid{i}))
%!               continue;
%!             endif
%!             elem_mat = zeros(rows(getfield(mesh.elements, elem_type_solid{i})), 1, "int32");
%!             elem_mat(getfield(mesh.groups, elem_type_solid{i})(grp_idx_volume(i)).elements) = 1;
%!             mesh.materials = setfield(mesh.materials, elem_type_solid{i}, elem_mat);
%!           endfor
%!           opt_mbd_mesh.forces.time_function = "time";
%!           if (px || py || pz)
%!             opt_mbd_mesh.pressure_loads.time_function = opt_mbd_mesh.forces.time_function;
%!             load_case.pressure = struct();
%!             load_case.traction = struct();
%!             for i=1:numel(grp_idx_load_px)
%!               if (~isfield(mesh.elements, elem_type_surf{i}))
%!                 continue;
%!               endif
%!               if (grp_idx_load_px(i) > 0)
%!                 elem_px = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_px(i)).elements;
%!               else
%!                 elem_px = [];
%!               endif
%!               if (grp_idx_load_py(i) > 0)
%!                 elem_py = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_py(i)).elements;
%!               else
%!                 elem_py = [];
%!               endif
%!               if (grp_idx_load_pz(i) > 0)
%!                 elem_pz = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_pz(i)).elements;
%!               else
%!                 elem_pz = [];
%!               endif
%!               if (grp_idx_load_mx(i) > 0)
%!                 elem_mx = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_mx(i)).elements;
%!               else
%!                 elem_mx = [];
%!               endif
%!               if (grp_idx_load_my(i) > 0)
%!                 elem_my = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_my(i)).elements;
%!               else
%!                 elem_my = [];
%!               endif
%!               if (grp_idx_load_mz(i) > 0)
%!                 elem_mz = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_mz(i)).elements;
%!               else
%!                 elem_mz = [];
%!               endif
%!               elem_nodes = [getfield(mesh.elements, elem_type_surf{i})(elem_px, :);
%!                             getfield(mesh.elements, elem_type_surf{i})(elem_py, :);
%!                             getfield(mesh.elements, elem_type_surf{i})(elem_pz, :);
%!                             getfield(mesh.elements, elem_type_surf{i})(elem_mx, :);
%!                             getfield(mesh.elements, elem_type_surf{i})(elem_my, :);
%!                             getfield(mesh.elements, elem_type_surf{i})(elem_mz, :)];
%!               switch (load_type{idx_load_type})
%!               case "pressure"
%!                 elem_press = [repmat(px, numel(elem_px), columns(elem_nodes));
%!                               repmat(py, numel(elem_py), columns(elem_nodes));
%!                               repmat(pz, numel(elem_pz), columns(elem_nodes));
%!                               repmat(px, numel(elem_mx), columns(elem_nodes));
%!                               repmat(py, numel(elem_my), columns(elem_nodes));
%!                               repmat(pz, numel(elem_mz), columns(elem_nodes))];
%!                 load_case.pressure = setfield(load_case.pressure, ...
%!                                               elem_type_surf{i}, ...
%!                                               struct("elements", elem_nodes, ...
%!                                                      "p", elem_press));
%!               case "traction"
%!                 elem_trac = zeros(rows(elem_nodes), columns(elem_nodes), 3);
%!                 ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 1) = -px;
%!                 ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 2) = -py;
%!                 ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 3) = -pz;
%!                 ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 1) =  px;
%!                 ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 2) =  py;
%!                 ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 3) =  pz;
%!                 ioffset += numel(elem_mz);
%!                 load_case.traction = setfield(load_case.traction, ...
%!                                               elem_type_surf{i}, ...
%!                                               struct("elements", elem_nodes, ...
%!                                                      "f", elem_trac));
%!               endswitch
%!             endfor
%!           endif
%!           opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!           opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!           opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!           idx_joint = int32(0);
%!           if (f_rigid_body_load || f_rigid_body_clamp)
%!             unwind_protect
%!               [fd, msg] = fopen(elem_file, "at");
%!               if (fd == -1)
%!                 error("failed to open file \"%s\": %s", elem_file, msg);
%!               endif
%!               if (f_rigid_body_load)
%!                 ++idx_joint;
%!                 fprintf(fd, "force: %d, absolute, %d,\n", ++opt_mbd_mesh.forces.number, inode_idx_rb(1));
%!                 fprintf(fd, "\tposition, reference, node, 0., 0., %.16e,\n", geometry.zF);
%!                 fputs(fd, "\tcomponent");
%!                 for i=1:3
%!                   fprintf(fd, ", mult, %s, const, %.16e", opt_mbd_mesh.forces.time_function, F(i));
%!                 endfor
%!                 fputs(fd, ";\n\n");
%!                 fprintf(fd, "couple: %d, absolute, %d,\n", ++opt_mbd_mesh.forces.number, inode_idx_rb(1));
%!                 fprintf(fd, "\tposition, reference, node, 0., 0., %.16e,\n", geometry.zF);
%!                 fputs(fd, "\tcomponent");
%!                 for i=1:3
%!                   fprintf(fd, ", mult, %s, const, %.16e", opt_mbd_mesh.forces.time_function, M(i));
%!                 endfor
%!                 fputs(fd, ";\n\n");
%!               endif
%!               if (f_rigid_body_clamp)
%!                 ++idx_joint;
%!                 fprintf(fd, "joint: %d, clamp, %d, node, node;\n\n", ++idx_joint, inode_idx_rb(2));
%!               endif
%!             unwind_protect_cleanup
%!               if (fd ~= -1)
%!                 fclose(fd);
%!               endif
%!               fd = -1;
%!             end_unwind_protect
%!           endif
%!           unwind_protect
%!             [fd, msg] = fopen(set_file, "wt");
%!             if (fd == -1)
%!               error("failed to open file \"%s\": %s", set_file, msg);
%!             endif
%!             fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!             fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!             fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!             fprintf(fd, "set: integer number_of_forces = %d;\n", opt_mbd_mesh.forces.number);
%!             fprintf(fd, "set: integer number_of_joints = %d;\n", idx_joint);
%!             fprintf(fd, "set: integer number_of_pressure_loads = %d;\n", opt_mbd_mesh.pressure_loads.number);
%!             fprintf(fd, "set: integer number_of_traction_loads = %d;\n", opt_mbd_mesh.traction_loads.number);
%!             fprintf(fd, "set: real t1 = %.16e;\n", t1);
%!             fprintf(fd, "set: real t2 = %.16e;\n", t2);
%!             fprintf(fd, "set: real dt = %.16e;\n", dt);
%!             fprintf(fd, "set: real dto = %.16e;\n", dto);
%!             for i=1:3
%!               fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, g(i));
%!             endfor
%!             for i=1:3
%!               fprintf(fd, "set: real W%s = %.16e;\n", {"x","y","z"}{i}, W(i));
%!             endfor
%!             for i=1:3
%!               fprintf(fd, "set: real WP%s = %.16e;\n", {"x","y","z"}{i}, WP(i));
%!             endfor
%!             for i=1:3
%!               fprintf(fd, "set: real XPP%s = %.16e;\n", {"x","y","z"}{i}, XPP(i));
%!             endfor
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!             fd = -1;
%!           end_unwind_protect
%!           unwind_protect
%!             [fd, msg] = fopen(control_file, "wt");
%!             if (fd == -1)
%!               error("failed to open file \"%s\": %s", control_file, msg);
%!             endif
%!             switch (model)
%!               case "static"
%!                 fprintf(fd, "model: %s;\n", model);
%!               case "dynamic"
%!                 fprintf(fd, "# model: %s;\n", model);
%!             endswitch
%!             if (norm(W) || norm(XPP) || norm(WP))
%!               fputs(fd, " rigid body kinematics: drive,\n");
%!               fputs(fd, " angular velocity, component,\n");
%!               fputs(fd, " mult, time, const, Wx / t1,\n");
%!               fputs(fd, " mult, time, const, Wy / t1,\n");
%!               fputs(fd, " mult, time, const, Wz / t1,\n");
%!               fputs(fd, " acceleration, component,\n");
%!               fputs(fd, " mult, time, const, XPPx / t1,\n");
%!               fputs(fd, " mult, time, const, XPPy / t1,\n");
%!               fputs(fd, " mult, time, const, XPPz / t1,\n");
%!               fputs(fd, " angular acceleration, component,\n");
%!               fputs(fd, " mult, time, const, WPx / t1,\n");
%!               fputs(fd, " mult, time, const, WPy / t1,\n");
%!               fputs(fd, " mult, time, const, WPz / t1;\n");
%!             endif
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!             fd = -1;
%!           end_unwind_protect
%!           unwind_protect
%!             [fd, msg] = fopen(initial_value_file, "wt");
%!             if (fd == -1)
%!               error("failed to open file \"%s\": %s", initial_value_file, msg);
%!             endif
%!             fprintf(fd, "method: %s;\n", method);
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!             fd = -1;
%!           end_unwind_protect
%!           fd = -1;
%!           unwind_protect
%!             [fd, msg] = fopen(input_file, "wt");
%!             if (fd == -1)
%!               error("failed to open file \"%s\": %s", input_file, msg);
%!             endif
%!             fprintf(fd, "include: \"%s\";\n", set_file);
%!             fprintf(fd, "begin: data;\n");
%!             fprintf(fd, "        problem: initial value; # the default\n");
%!             fprintf(fd, "end: data;\n");
%!             fprintf(fd, "begin: initial value;\n");
%!             fprintf(fd, "        initial time: 0;\n");
%!             fprintf(fd, "        final time: t1 + t2;\n");
%!             fprintf(fd, "        time step: dt;\n");
%!             fprintf(fd, "        time step: dt;\n");
%!             fprintf(fd, "        max iterations: 50;\n");
%!             fprintf(fd, "        tolerance: 1e-7, test, minmax, 1e-7, test, minmax;\n");
%!             fprintf(fd, "        output: messages;\n");
%!             fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!             fprintf(fd, "        nonlinear solver: nox,\n");
%!             fprintf(fd, "                          modified, 30,\n");
%!             fprintf(fd, "                          keep jacobian matrix,\n");
%!             fprintf(fd, "                          use preconditioner as solver, no,\n");
%!             fprintf(fd, "                          linesearch method, backtrack,\n");
%!             fprintf(fd, "                          direction, newton,\n");
%!             fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!             fprintf(fd, "                          forcing term, type 2,\n");
%!             fprintf(fd, "                          linear solver tolerance, 1e-8,\n");
%!             fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!             fprintf(fd, "                          linear solver max iterations, 300,\n");
%!             fprintf(fd, "                          krylov subspace size, 300,\n");
%!             fprintf(fd, "                          minimum step, 1e-6,\n");
%!             fprintf(fd, "                          recovery step type, constant,\n");
%!             fprintf(fd, "                          recovery step, 1e-6,\n");
%!             fprintf(fd, "                          verbose, 3,\n");
%!             fprintf(fd, "                          print convergence info, no;\n");
%!             fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!             fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!             fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!             fprintf(fd, "        derivatives max iterations: 10;\n");
%!             fprintf(fd, "        threads: assembly, 4;\n");
%!             fprintf(fd, "        threads: solver, 4;\n");
%!             fprintf(fd, "        output: cpu time;\n");
%!             fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!             fprintf(fd, "end: initial value;\n");
%!             fprintf(fd, "begin: control data;\n");
%!             fprintf(fd, "       skip initial joint assembly;\n");
%!             fprintf(fd, "       output precision: 16;\n");
%!             fprintf(fd, "       output meter: closest next, 0., forever, dto;\n");
%!             fprintf(fd, "       include: \"%s\";\n", control_file);
%!             fprintf(fd, "       default output: all, solids, accelerations;\n");
%!             fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!             fprintf(fd, "       solids: number_of_solids;\n");
%!             fprintf(fd, "       genels: number_of_genels;\n");
%!             fprintf(fd, "       forces: number_of_forces;\n");
%!             fprintf(fd, "       pressure loads: number_of_pressure_loads;\n");
%!             fprintf(fd, "       traction loads: number_of_traction_loads;\n");
%!             fprintf(fd, "       joints: number_of_joints;\n");
%!             fprintf(fd, "       gravity;\n");
%!             fprintf(fd, "       use automatic differentiation;\n");
%!             fprintf(fd, "end: control data;\n");
%!             fprintf(fd, "include: \"%s\";\n", csl_file);
%!             fprintf(fd, "begin: nodes;\n");
%!             fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!             fprintf(fd, "end: nodes;\n");
%!             fprintf(fd, "begin: elements;\n");
%!             fprintf(fd, "       include: \"%s\";\n", elem_file);
%!             fprintf(fd, "       gravity: uniform, gx, gy, gz, string, \"Time / t1\";\n");
%!             fprintf(fd, "end: elements;\n");
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!             fd = -1;
%!           end_unwind_protect
%!           if (options.verbose)
%!             shell(sprintf("cat \"%s\" | nl", input_file));
%!             shell(sprintf("cat \"%s\" | nl", nodes_file));
%!             shell(sprintf("cat \"%s\" | nl", csl_file));
%!             shell(sprintf("cat \"%s\" | nl", elem_file));
%!           endif
%!           info = mbdyn_solver_run(input_file, opt_mbd);
%!           [mesh_sol, sol] = mbdyn_post_load_output_sol(output_file);
%!           [genel_id, genel_data] = mbdyn_post_load_output([output_file, ".gen"], 1, [], numel(sol.t), 1);
%!           tau_ref = [-px; -py; -pz; 0; 0; 0];
%!           tol = 1e-11;
%!           for i=1:numel(elem_type_solid)
%!             if (~isfield(mesh.elements, elem_type_solid{i}))
%!               continue;
%!             endif
%!             tau_res = getfield(sol.stress.tau, elem_type_solid{i});
%!             for j=1:size(tau_res, 1)
%!               for k=1:size(tau_res, 2)
%!                 assert(tau_res(j, k, :, end)(:), tau_ref, tol * norm(tau_ref));
%!               endfor
%!             endfor
%!           endfor
%!           tol_F = 1e-10;
%!           Fref = max(abs([geometry.w * geometry.h * px, geometry.l * geometry.h * py, geometry.l * geometry.w * pz]));
%!           for i=1:numel(genel_data)
%!             assert(all(all(abs(genel_data{i}) < tol_F * Fref)));
%!           endfor
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
