## Copyright (C) 2020(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{load_case}, @var{bearing_surf}, @var{idx_modes}, @var{sol_eig}] = fem_ehd_pre_comp_mat_load_case2(@var{mesh}, @var{load_case}, @var{bearing_surf}, @var{options})
## Build pressure load cases required for fem_ehd_pre_comp_mat_unstruct.
##
## @var{mesh} @dots{} Return value from fem_pre_mesh_import or fem_pre_mesh_unstruct_create.
##
## @var{bearing_surf} @dots{} Struct array describing all bearing surfaces of this @var{mesh}.
##
## @var{bearing_surf}.group_idx @dots{} Index of a group of triangular elements in @var{mesh}.groups.tria6 used for applying pressure loads.
##
## @var{bearing_surf}.relative_tolerance @dots{} Optional field. If present, it will be checked if the nodes present in the @var{mesh} are covering the full bearing width.
##
## @var{bearing_surf}.absolute_tolerance @dots{} Optional field. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.X0 @dots{} Optional centre of the cylindrical bearing surface. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.R @dots{} Optional orientation of the cylindrical bearing surface. R(:, 3) will be the axis of the cylinder. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.w @dots{} Optional axial width of the cylindrical bearing surface. See also @var{bearing_surf}.relative_tolerance.
##
## @var{bearing_surf}.master_node_no @dots{} Master node number of an RBE3 element needed to apply rigid body displacements at the bearing surface if @var{options}.rigid_body_modes is equal to "flexible".
##
## @var{bearing_surf}.options.include_rigid_body_modes @dots{} This flag must be set to false if the bearing surface is related to the modal node of the body.
##
## @var{bearing_surf}.options.number_of_modes @dots{} Number of flexible modes used for each bearing surface respectively.
##
## @var{bearing_surf}.options.reference_pressure @dots{} Define the unit pressure applied to bearing surfaces.
##
## @var{bearing_surf}.options.mesh_size @dots{} Mesh size used for the hydrodynamic bearing element.
##
## @var{bearing_surf}.options.bearing_model @dots{} This value may be "EHD/FE" or "EHD/FD". In case of "EHD/FE" the number of output nodes will valid for a quadratic mesh.
##
## @var{options}.shift_A @dots{} If this value is not zero, an unconstrained eigenanalysis will be performed and the lowest six modes will be discarded.
##
## @var{options}.rigid_body_modes @dots{} This value may be "rigid" or "flexible". If rigid body modes are flexible, then RBE3 elements must be provided in the @var{mesh} and @var{bearing_surf}.master_node_no must be defined.
##
## @var{options}.rigid_body_modes_load_index @dots{} Consider all loads in @var{load_case}(@var{options}.rigid_body_modes_load_index) for rigid body modes.
##
## @var{options}.active_joint_idx_eig @dots{} Optional list of joints to be considered for eigenanalysis. By default, all joints are used.
##
## @var{options}.solver @dots{} Name of linear solver to use.
##
## @var{options}.number_of_threads @dots{} Number of threads to use for linear solver and assembly
##
## @var{options}.threshold_elem @dots{} Minimum number of elements for multithreaded assembly
##
## @var{options}.refine_max_iter @dots{} Maximum number of iterative refinement steps for the linear solver.
##
## @var{options}.verbose @dots{} Enable verbose output.
##
## @seealso{fem_ehd_pre_comp_mat_unstruct, fem_ehd_pre_comp_mat_load_case}
## @end deftypefn

function [mesh, load_case, bearing_surf, idx_modes, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, options)
  if (nargin < 3 || nargout > 5)
    print_usage();
  endif

  if (nargin < 4)
    options = struct();
  endif

  if (~isfield(options, "shift_A"))
    options.shift_A = 0;
  endif

  if (~isfield(options, "solver"))
    options.solver = fem_sol_select(true);
  endif

  if (~isfield(options, "number_of_threads"))
    options.number_of_threads = int32(1);
  endif

  if (~isfield(options, "threshold_elem"))
    options.threshold_elem = int32(10000);
  endif

  if (~isfield(options, "refine_max_iter"))
    options.refine_max_iter = int32(10);
  endif

  if (~isfield(options, "verbose"))
    options.verbose = int32(0);
  endif

  if (~isfield(options, "active_joint_idx_eig"))
    if (isfield(mesh.elements, "joints"))
      options.active_joint_idx_eig = 1:numel(mesh.elements.joints);
    else
      options.active_joint_idx_eig = [];
    endif
  endif

  if (~isfield(options, "rigid_body_modes"))
    options.rigid_body_modes = "flexible";
  endif

  if (~isfield(options, "rigid_body_modes_load_index"))
    options.rigid_body_modes_load_index = [];
  endif

  if (~isfield(options, "elem_type"))
    options.elem_type = "tria6";
  endif

  switch (options.rigid_body_modes)
    case "flexible"
      if (~isfield(bearing_surf, "master_node_no"))
        error("missing field bearing_surf.master_node_no needed for flexible rigid body modes");
      endif

      rigid_body_modes_flexible = true;
    case "rigid"
      rigid_body_modes_flexible = false;

      if (~isempty(options.rigid_body_modes_load_index))
        error("option rigid_body_modes_load_index can be used only with flexible interfaces");
      endif
    otherwise
      error("invalid option rigid_body_modes=\"%s\"", options.rigid_body_modes);
  endswitch

  options.solver = fem_sol_select(true, options.solver);
  options.symmetric = false;

  switch (options.solver)
    case {"umfpack", "lu", "mldivide"}
      mat_type_stiffness = FEM_MAT_STIFFNESS;
    otherwise
      mat_type_stiffness = FEM_MAT_STIFFNESS_SYM;
      options.symmetric = true;
  endswitch

  [load_case_press_dist, bearing_surf, idx_group] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf, options);

  load_case_press_dist(1).locked_dof = load_case.locked_dof;

  inum_elem_press = int32(0);

  for i=1:numel(bearing_surf)
    inum_elem_press += numel(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements);
    bearing_surf(i).nodes = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes;
  endfor

  elements = zeros(inum_elem_press, columns(getfield(mesh.elements, options.elem_type)), "int32");

  inum_elem_press = int32(0);
  ielem_idx = zeros(numel(bearing_surf), 2, "int32");

  for i=1:numel(bearing_surf)
    elem_idx = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements;
    idx_slot = inum_elem_press + (1:numel(elem_idx));
    elements(idx_slot, :) = getfield(mesh.elements, options.elem_type)(elem_idx, :);
    ielem_idx(i, :) = idx_slot([1, end]);
    inum_elem_press += numel(elem_idx);
  endfor

  load_case_press_mod = fem_pre_load_case_create_empty(numel(bearing_surf));

  for i=1:numel(bearing_surf)
    p = zeros(size(elements));
    p(ielem_idx(i, 1):ielem_idx(i, end), :) = 1;
    load_case_press_mod(i).pressure = struct(options.elem_type, struct("elements", elements, "p", p));
  endfor

  if (isfield(mesh.elements, "joints"))
    elem_joints = mesh.elements.joints;
  else
    elem_joints = struct("nodes", cell(1, 0), "C", cell(1, 0));
  endif

  mesh.elements.joints = elem_joints(options.active_joint_idx_eig);

  dof_map_press = fem_ass_dof_map(mesh, load_case(1));

  dof_map_press.parallel.threads_ass = options.number_of_threads;
  dof_map_press.parallel.threshold_elem = options.threshold_elem;

  [mat_ass_press.K, ...
   mat_ass_press.R] = fem_ass_matrix(mesh, ...
                                     dof_map_press, ...
                                     [mat_type_stiffness, ...
                                      FEM_VEC_LOAD_CONSISTENT], ...
                                     load_case_press_mod);


  mesh.elements.joints = elem_joints;

  diagA = zeros(rows(mesh.nodes), numel(bearing_surf));

  for i=1:columns(dof_map_press.ndof)
    dof_idx = dof_map_press.ndof(:, i);
    idx_act_dof = find(dof_idx > 0);
    diagA(idx_act_dof, :) += mat_ass_press.R(dof_idx(idx_act_dof), :).^2;
  endfor

  inum_modes_tot = numel(options.rigid_body_modes_load_index);
  inum_modes_press = int32(0);

  for i=1:numel(bearing_surf)
    inum_modes_tot += bearing_surf(i).options.number_of_modes;
    inum_modes_press += bearing_surf(i).options.number_of_modes;

    if (~isfield(bearing_surf(i).options, "include_rigid_body_modes"))
      bearing_surf(i).options.include_rigid_body_modes = true;
    endif

    if (bearing_surf(i).options.include_rigid_body_modes)
      inum_modes_tot += 6;
    endif
  endfor

  if (nargout >= 5)
    sol_eig.def = zeros(rows(mesh.nodes), columns(mesh.nodes), inum_modes_press);
    sol_eig.lambda = zeros(1, inum_modes_tot);
  endif

  sigma = 0;

  opt = struct();

  ir = int32([3, 2;
              1, 3;
              2, 1]);

  ic = int32([2, 3;
              3, 1;
              1, 2]);

  dc = [1, -1];

  ijoint_offset = numel(elem_joints);

  inum_joints = ijoint_offset;
  ijoint_idx = zeros(numel(bearing_surf), 2, "int32");

  for i=1:numel(bearing_surf)
    idx_joint = 1:numel(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes);
    ijoint_idx(i, :) = inum_joints + idx_joint([1, end]);
    inum_joints += numel(idx_joint);
  endfor

  empty_cell = cell(size(inum_joints));

  elem_joints_itf = struct("nodes", empty_cell, "C", empty_cell);

  elem_joints_itf(1:ijoint_offset) = elem_joints(1:ijoint_offset);

  for i=1:numel(bearing_surf)
    node_idx = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes;
    elem_joints_itf(ijoint_idx(i, 1) + (1:numel(node_idx)) - 1) = struct("nodes", mat2cell(node_idx, 1, ones(size(node_idx))), ...
                                                                         "C", repmat({[eye(3), zeros(3,3)]}, size(node_idx)));
  endfor

  load_case_displ = fem_pre_load_case_create_empty(inum_modes_tot);

  zero_U = struct("U", cellfun(@(C) zeros(rows(C), 1), {elem_joints_itf.C}, "UniformOutput", false));

  for i=1:numel(load_case_displ)
    load_case_displ(i).joints = zero_U;
  endfor

  inum_modes_tot = int32(0);
  inum_modes_press = int32(0);

  if (rigid_body_modes_flexible)
    mesh.elements.joints(ijoint_offset + (1:numel(bearing_surf))) = struct("nodes", {bearing_surf.master_node_no}, ...
                                                                           "C", repmat({eye(6)}, size(bearing_surf)));

    inum_modes_itf_flex = 6 * sum([[bearing_surf.options].include_rigid_body_modes]);

    load_case_itf_flex = fem_pre_load_case_create_empty(inum_modes_itf_flex);

    zero_U = struct("U", repmat({zeros(6, 1)}, 1, numel(bearing_surf)));

    for i=1:numel(options.rigid_body_modes_load_index)
      load_case(options.rigid_body_modes_load_index(i)).joints(ijoint_offset + (1:numel(bearing_surf))) = zero_U;
    endfor

    zero_U = struct("U", repmat({zeros(6, 1)}, size(mesh.elements.joints)));

    num_load_case_itf_flex = int32(0);
    idx_load_case_itf_flex = zeros(numel(bearing_surf), 2, "int32");

    for i=1:numel(bearing_surf)
      if (~bearing_surf(i).options.include_rigid_body_modes)
        continue;
      endif

      idx_load_case_itf_flex(i, :) = num_load_case_itf_flex + [1, 6];

      for j=1:6
        load_case_itf_flex(++num_load_case_itf_flex).joints = zero_U;
        load_case_itf_flex(num_load_case_itf_flex).joints(ijoint_offset + i).U(j) = 1;
      endfor
    endfor

    if (~isempty(options.rigid_body_modes_load_index))
      load_case_itf_flex = fem_pre_load_case_merge(load_case(options.rigid_body_modes_load_index), load_case_itf_flex);
    endif

    idx_load_case_itf_flex += numel(options.rigid_body_modes_load_index);

    clear zero_U;

    dof_map_itf_flex = fem_ass_dof_map(mesh, load_case(1));

    dof_map_itf_flex.parallel.threads_ass = options.number_of_threads;
    dof_map_itf_flex.threshold_elem = options.threshold_elem;

    [mat_ass_itf_flex.K, ...
     mat_ass_itf_flex.R] = fem_ass_matrix(mesh, ...
                                          dof_map_itf_flex, ...
                                          [mat_type_stiffness, ...
                                           FEM_VEC_LOAD_CONSISTENT], ...
                                          load_case_itf_flex);

    Kfact = fem_sol_factor(mat_ass_itf_flex.K, options);

    U = Kfact \ mat_ass_itf_flex.R;

    clear Kfact mat_ass_itf_flex load_case_itf_flex;

    Urb = cell(1, numel(bearing_surf));

    for i=1:numel(bearing_surf)
      inode_idx_bs = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes;
      idx_disp = (numel(options.rigid_body_modes_load_index) + 1):columns(U);
      Ubs = zeros(numel(inode_idx_bs) * 3, columns(U));

      for j=1:3
        idof_idx = dof_map_itf_flex.ndof(inode_idx_bs, j);
        idx_act_dof = find(idof_idx > 0);
        Ubs((idx_act_dof - 1) * 3 + j, :) = U(idof_idx(idx_act_dof), :);
      endfor

      if (bearing_surf(i).options.include_rigid_body_modes)
        Urb{i} = Ubs(:, idx_load_case_itf_flex(i, 1):idx_load_case_itf_flex(i, 2));
      endif

      for k=1:columns(Ubs)
        load_case_displ(inum_modes_tot + k).joints(ijoint_idx(i, 1) + (1:numel(inode_idx_bs)) - 1) = struct("U", mat2cell(Ubs(:, k), repmat(3, size(inode_idx_bs)), 1));
      endfor
    endfor

    for i=1:numel(options.rigid_body_modes_load_index)
      load_case_displ(inum_modes_tot + i).joints(1:ijoint_offset) = load_case(options.rigid_body_modes_load_index(i)).joints(1:ijoint_offset);
    endfor

    inum_modes_tot += columns(U);
  endif

  if (options.shift_A == 0)
    Kfact = fem_sol_factor(mat_ass_press.K, options);
  endif

  for i=1:numel(bearing_surf)
    Ap = zeros(dof_map_press.totdof, 1);

    for j=1:columns(dof_map_press.ndof)
      dof_idx = dof_map_press.ndof(:, j);
      idx_act_dof = find(dof_idx > 0);
      Ap(dof_idx(idx_act_dof)) = sqrt(diagA(idx_act_dof, i));
    endfor

    Ap = diag(Ap);

    if (options.shift_A == 0)
      gamma = 0;
    else
      gamma = options.shift_A * max(max(abs(mat_ass_press.K))) / max(abs(diag(Ap)));
      Kfact = fem_sol_factor(mat_ass_press.K + gamma * Ap, options);
    endif

    oper = cell(1, 2);
    oper{1} = @(x) Ap * x;
    oper{2} = @(x) Kfact \ x;

    num_modes = bearing_surf(i).options.number_of_modes;

    keep_rigid_body_normal_modes = options.shift_A == 0;

    if (~keep_rigid_body_normal_modes)
      num_modes += 6;
    endif

    if (num_modes)
      [Phi, kappa] = eig_sym(oper, columns(mat_ass_press.K), num_modes, sigma, opt);
    else
      Phi = zeros(rows(mat_ass_press.K), 0);
      kappa = [];
    endif

    clear oper Ap;

    if (options.shift_A ~= 0)
      clear Kfact;
    endif

    inode_idx_bs = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes;

    U = zeros(numel(inode_idx_bs) * 3, columns(Phi));

    for j=1:3
      dof_idx = dof_map_press.ndof(inode_idx_bs, j);
      idx_act_dof = find(dof_idx > 0);
      U((idx_act_dof - 1) * 3 + j, :) = Phi(dof_idx(idx_act_dof), :);
    endfor

    if (bearing_surf(i).options.include_rigid_body_modes && rigid_body_modes_flexible)
      Arb = Urb{i};
    else
      Arb = [repmat(eye(3), numel(inode_idx_bs), 1), zeros(3 * numel(inode_idx_bs), 3)];

      l = mesh.nodes(inode_idx_bs, 1:3) - bearing_surf(i).X0.';

      for j=1:3
        for k=1:2
          Arb(ir(j, k):3:end, ic(j, k) + 3) = -dc(k) * l(:, j);
        endfor
      endfor
    endif

    U -= Arb * (Arb \ U);

    if (~keep_rigid_body_normal_modes)
      U = U(:, 7:end);
    endif

    if (bearing_surf(i).options.include_rigid_body_modes && ~rigid_body_modes_flexible)
      U = [Arb, U];
    endif

    U *= diag(1 ./ max(abs(U), [], 1));

    for k=1:columns(U)
      load_case_displ(inum_modes_tot + k).joints(ijoint_idx(i, 1) + (1:numel(inode_idx_bs)) - 1) = struct("U", mat2cell(U(:, k), repmat(3, size(inode_idx_bs)), 1));
    endfor

    if (nargout >= 5)
      idx_mode_press = inum_modes_press + (1:columns(Phi));

      for j=1:columns(dof_map_press.ndof)
        idof_idx = dof_map_press.ndof(:, j);
        iactive_dof = find(idof_idx > 0);
        idof_idx = idof_idx(iactive_dof);

        sol_eig.def(iactive_dof, j, idx_mode_press) = Phi(idof_idx, :);
        clear idof_idx iactive_dof;
      endfor

      sol_eig.lambda(idx_mode_press) = diag(kappa) - gamma;
      inum_modes_press += columns(Phi);
    endif

    inum_modes_tot += columns(U);
    clear U;
  endfor

  clear mat_ass_press dof_map_press;

  mesh.elements.joints = elem_joints_itf;

  zero_U = struct("U", cellfun(@(C) zeros(rows(C), 1), {elem_joints_itf.C}, "UniformOutput", false));

  for i=1:numel(load_case_press_dist)
    load_case_press_dist(i).joints = zero_U;
  endfor

  load_case = fem_pre_load_case_merge(load_case_press_dist, load_case_displ);

  idx_modes = numel(load_case_press_dist) + (1:numel(load_case_displ));

  if (nargout >= 5)
    for i=1:size(sol_eig.def, 3)
      sol_eig.def(:, :, i) /= max(max(abs(sol_eig.def(:, :, i))));
    endfor
  endif
endfunction

