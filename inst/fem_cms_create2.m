## Copyright (C) 2019(-2025) Reinhard <octave-user@a1.net>
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
## @var{cms_opt}.modes.number @dots{} Number of normal modes to be considered
##
## @var{cms_opt}.element.name @dots{} MBDyn element name
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

  if (~isfield(cms_opt, "solver"))
    cms_opt.solver = fem_sol_select(true);
  endif

  for i=1:numel(cms_opt.nodes.interfaces)
    if (~isfield(cms_opt.nodes.interfaces(i), "enable_KTAU0") || isempty(cms_opt.nodes.interfaces(i).enable_KTAU0))
      cms_opt.nodes.interfaces(i).enable_KTAU0 = true(6, 1);
    endif
  endfor

  mat_type_stiffness = FEM_MAT_STIFFNESS_SYM_L;

  cms_opt.solver = fem_sol_select(true, cms_opt.solver);

  switch (cms_opt.solver)
    case {"umfpack", "lu", "mldivide"}
      mat_type_stiffness = FEM_MAT_STIFFNESS; ## enforce full matrix for unsymmetric solvers
  endswitch

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
                                           [mat_type_stiffness, ...
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

    load_case_tau0.tau0 = sol_tau.stress.tau;
    load_case_tau0.lambda = sol_tau.lambda;

    [~, KTAU0{i}] = fem_ass_matrix(mesh, ...
                                   dof_map_tau, ...
                                   [FEM_MAT_STIFFNESS, ...
                                    FEM_MAT_STIFFNESS_TAU0], ...
                                   load_case_tau0);
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
                                              mat_type_stiffness, ...
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
                                              mat_type_stiffness, ...
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

  switch (mat_type_stiffness)
    case {FEM_MAT_STIFFNESS_SYM, FEM_MAT_STIFFNESS_SYM_L}
      Ksym = fem_mat_sym(mat_ass_itf.K);
    otherwise
      Ksym = mat_ass_itf.K;
  endswitch

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
