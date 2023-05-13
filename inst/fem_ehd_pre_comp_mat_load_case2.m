## Copyright (C) 2020(-2022) Reinhard <octave-user@a1.net>
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
    options.solver = "pastix";
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
      mat_type_stiffness = FEM_MAT_STIFFNESS_SYM_L;
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

%!test
%!  ## TEST 1
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    f_enable_constraint = [false, true];
%!    interfaces = {"flexible","rigid"};
%!    num_modes = int32([0, 10]);
%!    for l=1:numel(num_modes)
%!      for k=1:numel(f_enable_constraint)
%!        for j=1:numel(interfaces)
%!          clear bearing_surf cms_opt comp_mat dof_map grp_id_clamp grp_id_p1 grp_id_p2;
%!          clear load_case load_case_bearing mat_ass mat_ass_press mat_info mesh mesh_info mesh_size
%!          clear opt_modes sol_eig sol_eig_cms sol_stat;
%!          unwind_protect
%!            [fd, msg] = fopen([filename, ".geo"], "w");
%!            if (fd == -1)
%!              error("failed to open file \"%s.geo\"", filename);
%!            endif
%!            ri = 8e-3;
%!            ro = 10e-3;
%!            h = 12e-3;
%!            c = 2e-3;
%!            b = h - 2 * c;
%!            scale_def = 5e-3;
%!            mesh_size = 3e-3;
%!            fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!            fprintf(fd, "ri = %g;\n", ri);
%!            fprintf(fd, "ro = %g;\n", ro);
%!            fprintf(fd, "h = %g;\n", h);
%!            fprintf(fd, "c = %g;\n", c);
%!            fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!            fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!            fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!            fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!            fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!            fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!            fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!            fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!            fputs(fd, "Line(1) = {1,2};\n");
%!            fputs(fd, "Line(2) = {2,3};\n");
%!            fputs(fd, "Line(3) = {3,4};\n");
%!            fputs(fd, "Line(4) = {4,5};\n");
%!            fputs(fd, "Line(5) = {5,6};\n");
%!            fputs(fd, "Line(6) = {6,7};\n");
%!            fputs(fd, "Line(7) = {7,8};\n");
%!            fputs(fd, "Line(8) = {8,1};\n");
%!            fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!            fputs(fd, "Plane Surface(6) = {5};\n");
%!            fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%!            fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!            fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!            fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!            fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!            fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!            fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!          unwind_protect_cleanup
%!            if (fd ~= -1)
%!              fclose(fd);
%!              fd = -1;
%!            endif
%!          end_unwind_protect
%!          fprintf(stderr, "meshing ...\n");
%!          pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!          status = spawn_wait(pid);
%!          if (status ~= 0)
%!            error("gmsh failed with status %d", status);
%!          endif
%!          fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!          mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!          fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!          cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!          cms_opt.solver = "mldivide";
%!          switch (interfaces{j})
%!            case "flexible"
%!              cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!          endswitch
%!          grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%!          grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!          grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%!          bearing_surf(1).group_idx = grp_id_p1;
%!          bearing_surf(1).options.reference_pressure = 1e9;
%!          bearing_surf(1).options.mesh_size = 20e-3;
%!          bearing_surf(1).options.include_rigid_body_modes = false;
%!          bearing_surf(1).options.bearing_type = "shell";
%!          bearing_surf(1).options.matrix_type = "modal substruct total";
%!          bearing_surf(1).r = ri;
%!          bearing_surf(1).w = b;
%!          bearing_surf(1).X0 = [0; 0; b/2 + c];
%!          bearing_surf(1).R = eye(3);
%!          bearing_surf(1).relative_tolerance = 0;
%!          bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%!          bearing_surf(1).options.number_of_modes = 10;
%!          bearing_surf(2).group_idx = grp_id_p2;
%!          bearing_surf(2).options.reference_pressure = 1e9;
%!          bearing_surf(2).options.mesh_size = 20e-3;
%!          bearing_surf(2).options.include_rigid_body_modes = true;
%!          bearing_surf(2).options.bearing_type = "journal";
%!          bearing_surf(2).options.matrix_type = "modal substruct total";
%!          bearing_surf(2).r = ro;
%!          bearing_surf(2).w = b;
%!          bearing_surf(2).X0 = [0; 0; b/2 + c];
%!          bearing_surf(2).R = eye(3);
%!          bearing_surf(2).relative_tolerance = 0;
%!          bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%!          bearing_surf(2).options.number_of_modes = 12;
%!          switch (interfaces{j})
%!            case "flexible"
%!              bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!              bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!              for i=1:numel(bearing_surf)
%!                mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!              endfor
%!              for i=1:numel(bearing_surf)
%!                mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no);
%!              endfor
%!            otherwise
%!              mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!          endswitch
%!          cms_opt.inveriants = true;
%!          cms_opt.modes.number = num_modes(l);
%!          cms_opt.static_modes = false;
%!          cms_opt.modal_node_constraint = false;
%!          cms_opt.load_cases = "index";
%!          cms_opt.refine_max_iter = int32(10);
%!          load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!          switch (interfaces{j})
%!            case "flexible"
%!            otherwise
%!              load_case(1).locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%!          endswitch
%!          if (f_enable_constraint(k))
%!            load_case(1).locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%!          endif
%!          mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!          mesh.material_data.E = 210000e6;
%!          mesh.material_data.nu = 0.3;
%!          mesh.material_data.rho = 7850;
%!          if (f_enable_constraint(k))
%!            opt_modes.shift_A = 0;
%!          else
%!            opt_modes.shift_A = 1e-6;
%!          endif
%!          opt_modes.refine_max_iter = int32(10);
%!          opt_modes.verbose = int32(0);
%!          opt_modes.solver = cms_opt.solver;
%!          opt_modes.rigid_body_modes = interfaces{j};
%!          [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          dof_map = fem_ass_dof_map(mesh, load_case);
%!          [mat_ass.K, ...
%!           mat_ass.M, ...
%!           mat_ass.R, ...
%!           mat_info, ...
%!           mesh_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_STIFFNESS, ...
%!                                        FEM_MAT_MASS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case_bearing);
%!          sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!          [mesh, ...
%!           mat_ass, ...
%!           dof_map, ...
%!           sol_eig_cms, ...
%!           cms_optp] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_optp, bearing_surf);
%!        endfor
%!      endfor
%!    endfor
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 2
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    interfaces = {"rigid", "flexible"};
%!    tol_red = [5e-2, 2e-2];
%!    num_modes_cms = int32([0, 10]);
%!    num_modes = [50, 30];
%!    for k=1:numel(num_modes_cms)
%!      for j=1:numel(interfaces)
%!        clear Fred Ritf bearing_surf cms_opt comp_mat dof_map_comb dof_map_post err_red
%!        clear grp_id_p1 grp_id_p2 grp_idx_p1 grp_idx_p2 load_case load_case_bearing load_case_itf
%!        clear load_case_post mat_ass_post mat_ass_press mesh mesh_comb mesh_data mesh_post mesh_size
%!        clear opt_modes p1 p1red p2 p2red pid qred sol_comb sol_eig sol_eig_cms sol_post sol_red
%!        unwind_protect
%!          [fd, msg] = fopen([filename, ".geo"], "w");
%!          if (fd == -1)
%!            error("failed to open file \"%s.geo\"", filename);
%!          endif
%!          d = 14e-3;
%!          D = 19.5e-3;
%!          w = 5e-3;
%!          l = 47e-3;
%!          h = 5e-3;
%!          grp_id_p1 = 2;
%!          grp_id_p2 = 3;
%!          p1 = 1;
%!          p2 = 2;
%!          scale_def = 5e-3;
%!          mesh_size = 2e-3;
%!          fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!          fprintf(fd, "d = %g;\n", d);
%!          fprintf(fd, "D = %g;\n", D);
%!          fprintf(fd, "w = %g;\n", w);
%!          fprintf(fd, "l = %g;\n", l);
%!          fprintf(fd, "h = %g;\n", h);
%!          fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!          fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!          fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!          fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!          fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!          fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!          fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!          fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!          fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!          fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!          fputs(fd, "Line(11) = {13, 14};\n");
%!          fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!          fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!          fputs(fd, "Line(14) = {16, 11};\n");
%!          fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!          fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!          fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!          fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!          fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!          fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!          fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!          fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!          fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!          fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!        unwind_protect_cleanup
%!          if (fd ~= -1)
%!            fclose(fd);
%!            fd = -1;
%!          endif
%!        end_unwind_protect
%!        fprintf(stderr, "meshing ...\n");
%!        pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!        status = spawn_wait(pid);
%!        if (status ~= 0)
%!          error("gmsh failed with status %d", status);
%!        endif
%!        fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!        mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!        fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!        grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!        grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!        cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!        cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!        bearing_surf(1).group_idx = grp_idx_p1;
%!        bearing_surf(1).options.reference_pressure = 1e9;
%!        bearing_surf(1).options.mesh_size = 1e-3;
%!        bearing_surf(1).options.include_rigid_body_modes = true;
%!        bearing_surf(1).options.bearing_type = "shell";
%!        bearing_surf(1).options.matrix_type = "modal substruct total";
%!        bearing_surf(1).r = 0.5 * d;
%!        bearing_surf(1).w = w;
%!        bearing_surf(1).X0 = [l; 0; 0];
%!        bearing_surf(1).R = eye(3);
%!        bearing_surf(1).relative_tolerance = 0;
%!        bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!        bearing_surf(1).options.number_of_modes = num_modes(j);
%!        bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!        bearing_surf(2).group_idx = grp_idx_p2;
%!        bearing_surf(2).options.reference_pressure = 1e9;
%!        bearing_surf(2).options.mesh_size = 1e-3;
%!        bearing_surf(2).options.include_rigid_body_modes = false;
%!        bearing_surf(2).options.bearing_type = "shell";
%!        bearing_surf(2).options.matrix_type = "modal substruct total";
%!        bearing_surf(2).r = 0.5 * d;
%!        bearing_surf(2).w = w;
%!        bearing_surf(2).X0 = [0; 0; 0];
%!        bearing_surf(2).R = eye(3);
%!        bearing_surf(2).relative_tolerance = 0;
%!        bearing_surf(2).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!        bearing_surf(2).options.number_of_modes = num_modes(j);
%!        bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!        mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!        mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!        mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!        mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!        cms_opt.inveriants = true;
%!        cms_opt.modes.number = num_modes_cms(k);
%!        cms_opt.static_modes = false;
%!        cms_opt.modal_node_constraint = false;
%!        cms_opt.load_cases = "index";
%!        cms_opt.refine_max_iter = int32(10);
%!        load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!        mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!        mesh.material_data.E = 210000e6;
%!        mesh.material_data.nu = 0.3;
%!        mesh.material_data.rho = 7850;
%!        opt_modes.shift_A = 1e-6;
%!        opt_modes.refine_max_iter = int32(10);
%!        opt_modes.verbose = int32(0);
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        opt_modes.solver = "umfpack";
%!        cms_opt.solver = opt_modes.solver;
%!        [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!        [mesh, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!        assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!        assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!        comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%!        load_case_itf = fem_pre_load_case_create_empty(6);
%!        for i=1:6
%!          load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!          load_case_itf(i).loads = zeros(1, 6);
%!          load_case_itf(i).loads(i) = 1;
%!        endfor
%!        [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!        nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!        nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!        nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!        nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!        p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!        p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!        p1red2 = zeros((nx1 - 1) * nz1, 1);
%!        p2red2 = zeros((nx2 - 1) * nz2, 1);
%!        for i=1:nx1 - 1
%!          p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!        endfor
%!        for i=1:nx2 - 1
%!          p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!        endfor
%!        Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!                comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!        Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!        qred = mat_ass.Kred \ Fred;
%!        sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!        mesh_post = mesh;
%!        mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!        load_case_post = fem_pre_load_case_create_empty(10);
%!        load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!        for i=1:6
%!          load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!          load_case_post(i).loads = zeros(1, 6);
%!          load_case_post(i).loads(i) = 1;
%!        endfor
%!        x1 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!        y1 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!        x2 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!        y2 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!        Phi1 = atan2(y1, x1);
%!        Phi2 = atan2(y2, x2);
%!        load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!        load_case_post(7).pressure.tria6.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%!        load_case_post(8).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!        load_case_post(8).pressure.tria6.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!        load_case_post(9).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!        load_case_post(9).pressure.tria6.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6);
%!        load_case_post(10).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!        load_case_post(10).pressure.tria6.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!        mesh_post.elements.joints.C = eye(6);
%!        mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!        dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!        [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!        sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, opt_modes);
%!        mesh_data(1).mesh = mesh;
%!        mesh_data(1).dof_map = dof_map;
%!        mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!        mesh_data(2).mesh = mesh_post;
%!        mesh_data(2).dof_map = dof_map_post;
%!        [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!        sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!        sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!        sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!        for i=1:size(sol_comb.def, 3)
%!          sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!        endfor
%!        err_red = zeros(1, size(sol_post.def, 3));
%!        for i=1:size(sol_post.def, 3)
%!          err_red(i) = max(max(abs(sol_post.def(1:end - 2, :, i) - sol_red.def(1:end - 2, :, i)))) / max(max(abs(sol_post.def(1:end - 2, :, i))));
%!        endfor
%!        for i=1:size(sol_post.def, 3)
%!          fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!        endfor
%!        assert(all(err_red < tol_red(j)));
%!      endfor
%!    endfor
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 3
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!        error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      d1 = 8e-3;
%!      D1 = 12e-3;
%!      d2 = 14e-3;
%!      D2 = 19.5e-3;
%!      w = 5e-3;
%!      l = 47e-3;
%!      h = 5e-3;
%!      grp_id_p1 = 2;
%!      grp_id_p2 = 3;
%!      scale_def = 5e-3;
%!      mesh_size = 1.5e-3;
%!      fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "d1 = %g;\n", d1);
%!      fprintf(fd, "D1 = %g;\n", D1);
%!      fprintf(fd, "d2 = %g;\n", d2);
%!      fprintf(fd, "D2 = %g;\n", D2);
%!      fprintf(fd, "w = %g;\n", w);
%!      fprintf(fd, "l = %g;\n", l);
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(6)  = {         0.0,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(7)  = {         0.0,  0.5 * d2, -0.5 * w};\n");
%!      fputs(fd, "Point(8)  = {    0.5 * d2,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(9)  = {         0.0, -0.5 * d2, -0.5 * w};\n");
%!      fputs(fd, "Point(10) = {   -0.5 * d2,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(11) = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(12) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(13) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(14) = {Sqrt((D2/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(15) = {   -0.5 * D2,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(16) = {Sqrt((D2/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!      fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!      fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!      fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!      fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!      fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!      fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!      fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!      fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!      fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!      fputs(fd, "Line(11) = {13, 14};\n");
%!      fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!      fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!      fputs(fd, "Line(14) = {16, 11};\n");
%!      fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!      fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!      fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!      fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!      fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!      fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!      fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!      fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!      fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!      fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!    unwind_protect_cleanup
%!      if (fd ~= -1)
%!        fclose(fd);
%!        fd = -1;
%!      endif
%!    end_unwind_protect
%!    fprintf(stderr, "meshing ...\n");
%!    pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!    status = spawn_wait(pid);
%!    if (status ~= 0)
%!      error("gmsh failed with status %d", status);
%!    endif
%!    fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!    mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!    fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!    grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!    grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!    cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!    cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!    bearing_surf(1).group_idx = grp_idx_p1;
%!    bearing_surf(1).options.reference_pressure = 1e9;
%!    bearing_surf(1).options.mesh_size = 1e-3;
%!    bearing_surf(1).options.include_rigid_body_modes = true;
%!    bearing_surf(1).options.bearing_type = "shell";
%!    bearing_surf(1).options.matrix_type = "modal substruct total";
%!    bearing_surf(1).r = 0.5 * d1;
%!    bearing_surf(1).w = w;
%!    bearing_surf(1).X0 = [l; 0; 0];
%!    bearing_surf(1).R = eye(3);
%!    bearing_surf(1).relative_tolerance = 0;
%!    bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d1;
%!    bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!    bearing_surf(2).group_idx = grp_idx_p2;
%!    bearing_surf(2).options.reference_pressure = 1e9;
%!    bearing_surf(2).options.mesh_size = 1e-3;
%!    bearing_surf(2).options.include_rigid_body_modes = false;
%!    bearing_surf(2).options.bearing_type = "shell";
%!    bearing_surf(2).options.matrix_type = "modal substruct total";
%!    bearing_surf(2).r = 0.5 * d2;
%!    bearing_surf(2).w = w;
%!    bearing_surf(2).X0 = [0; 0; 0];
%!    bearing_surf(2).R = eye(3);
%!    bearing_surf(2).relative_tolerance = 0;
%!    bearing_surf(2).absolute_tolerance = sqrt(eps) * 0.5 * d2;
%!    bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!    mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!    mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!    mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!    cms_opt.inveriants = true;
%!    cms_opt.static_modes = false;
%!    cms_opt.modal_node_constraint = false;
%!    cms_opt.load_cases = "index";
%!    cms_opt.solver = "umfpack";
%!    cms_opt.refine_max_iter = int32(10);
%!    load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!    mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!    mesh.material_data.E = 210000e6;
%!    mesh.material_data.nu = 0.3;
%!    mesh.material_data.rho = 7850;
%!    opt_modes.shift_A = 1e-6;
%!    opt_modes.refine_max_iter = int32(10);
%!    opt_modes.verbose = int32(0);
%!    opt_modes.solver = cms_opt.solver;
%!    num_modes = [0, 30];
%!    interfaces = {"rigid", "flexible"};
%!    num_modes_cms = int32([0, 10]);
%!    k1 = 1;
%!    k2 = 0;
%!    num_cases = numel(num_modes) * numel(interfaces) * numel(num_modes_cms);
%!    err_red = zeros(7, num_cases);
%!    err_mod = zeros(1, num_cases);
%!    err_w = zeros(1, num_cases);
%!    icase = int32(0);
%!    for k=1:numel(num_modes_cms)
%!      cms_opt.modes.number = num_modes_cms(k);
%!      for j=1:numel(interfaces)
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        for l=1:numel(num_modes)
%!          for i=1:numel(bearing_surf)
%!            bearing_surf(i).options.number_of_modes = min(num_modes(l), floor(numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes) * 3 / 2 - 6));
%!          endfor
%!          [mesh_l, load_case_bearing, bearing_surf_l, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          [mesh_l, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh_l, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!          [lambda, idx_lambda] = sort(diag(lambda));
%!          qred = qred(:, idx_lambda);
%!          sol_eig_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          for i=1:size(sol_eig_red.def, 3)
%!            sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!          endfor
%!          sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh_l, mat_ass, dof_map, cms_opt, bearing_surf_l);
%!          load_case_itf = fem_pre_load_case_create_empty(6);
%!          for i=1:6
%!            load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!            load_case_itf(i).loads = zeros(1, 6);
%!            load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!            switch (i)
%!              case {4, 5}
%!                load_case_itf(i).loads(i) *= w;
%!              case 6
%!                load_case_itf(i).loads(i) *= d1;
%!            endswitch
%!          endfor
%!          [~, Ritf] = fem_ass_matrix(mesh_l, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!          nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!          nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!          nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!          nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!          X1 = mesh_l.nodes(mesh_l.groups.tria6(bearing_surf_l(1).group_idx).nodes, 1:3) - bearing_surf_l(1).X0.';
%!          Phi1 = atan2(X1(:, 2), X1(:, 1));
%!          p1 = k1 * sin(Phi1).^2 * bearing_surf_l(1).options.reference_pressure;
%!          Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!          z1g = comp_mat(1).bearing_surf.grid_z(:);
%!          p1red = zeros(numel(z1g), numel(Phi1g));
%!          for i=1:columns(p1red)
%!            p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!          endfor
%!          p1red = p1red(:);
%!          X2 = mesh_l.nodes(mesh_l.groups.tria6(bearing_surf_l(2).group_idx).nodes, 1:3) - bearing_surf_l(2).X0.';
%!          Phi2 = atan2(X2(:, 2), X2(:, 1));
%!          p2 = k2 * cos(Phi2) .* sin(pi * X2(:, 3) / (max(X2(:, 3)) - min(X2(:, 3)))) * bearing_surf_l(2).options.reference_pressure;
%!          Phi2g = comp_mat(2).bearing_surf.grid_x(:) / (0.5 * comp_mat(2).bearing_dimensions.bearing_diameter);
%!          z2g = comp_mat(2).bearing_surf.grid_z(:);
%!          p2red = zeros(numel(z2g), numel(Phi2g));
%!          for i=1:columns(p2red)
%!            p2red(:, i) = griddata([Phi2 - 2 * pi; Phi2; Phi2 + 2 * pi], repmat(X2(:, 3), 3, 1), repmat(p2, 3, 1), repmat(Phi2g(i), rows(z2g), 1), z2g);
%!          endfor
%!          p2red = p2red(:);
%!          Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf_l(1).options.reference_pressure + comp_mat(2).E(:, 1:end - nz2) * p2red(1:end -  nz2) / bearing_surf_l(2).options.reference_pressure;
%!          Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!          qred = mat_ass.Kred \ Fred;
%!          w1red = comp_mat(1).D * qred;
%!          w2red = comp_mat(2).D * qred;
%!          sol_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          mesh_post = mesh_l;
%!          mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!          load_case_post = fem_pre_load_case_create_empty(7);
%!          load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!          for i=1:6
%!            load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!            load_case_post(i).loads = load_case_itf(i).loads;
%!          endfor
%!          load_case_post(7).pressure.tria6.elements = [mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!                                                       mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :)];
%!          pn = zeros(rows(mesh_l.nodes), 1);
%!          pn(mesh_post.groups.tria6(grp_idx_p1).nodes) = p1;
%!          pn(mesh_post.groups.tria6(grp_idx_p2).nodes) = p2;
%!          load_case_post(7).pressure.tria6.p = [pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :));
%!                                                pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :))];
%!          mesh_post.elements.joints.C = eye(6);
%!          mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!          dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!          [mat_ass_post.M, ...
%!           mat_ass_post.K, ...
%!           mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                            dof_map_post, ...
%!                                            [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_post);
%!          sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!          node_idx1 = mesh_l.groups.tria6(bearing_surf(1).group_idx).nodes;
%!          w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!          for i=1:numel(node_idx1)
%!            ni = [mesh_l.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!            ni /= norm(ni);
%!            w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!          endfor
%!          w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!          for i=1:columns(w1postint)
%!            for m=1:numel(Phi1g)
%!              w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                             repmat(X1(:, 3), 3, 1), ...
%!                                                                             repmat(w1post(:, i), 3, 1), ...
%!                                                                             repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                             z1g);
%!            endfor
%!          endfor
%!          if (do_plot)
%!            for m=1:columns(w1red)
%!              for i=1:numel(z1g)
%!                figure("visible", "off");
%!                hold("on");
%!                plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;1");
%!                plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;0");
%!                ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!                title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!                xlabel("Phi [deg]");
%!                ylabel("w [um]");
%!              endfor
%!            endfor
%!          endif
%!          sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!          for i=1:size(sol_eig_post.def, 3)
%!            sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!          endfor
%!          num_modes_comp = min(10, floor(numel(sol_eig_red.f) / 3));
%!          err_mod(++icase) = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!          err_w(icase) = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!          mesh_data(1).mesh = mesh_l;
%!          mesh_data(1).dof_map = dof_map;
%!          mesh_data(2).mesh = mesh_post;
%!          mesh_data(2).dof_map = dof_map_post;
%!          mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!          [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!          sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!          for i=1:size(sol_post.def, 3)
%!            err_red(i, icase) = max(max(abs(sol_post.def(1:end - 2, 1:3, i) - sol_red.def(1:end - 2, 1:3, i)))) / max(max(abs(sol_post.def(1:end - 2, 1:3, i))));
%!          endfor
%!        endfor
%!        for l=1:numel(num_modes)
%!          fprintf(stderr, "%s interfaces using %d static pressure modes, %d dynamic modes:\n", interfaces{j}, num_modes(l), num_modes_cms(k));
%!          for i=1:size(sol_post.def, 3)
%!            fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i, icase));
%!          endfor
%!          fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod(icase));
%!          fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w(icase));
%!        endfor
%!      endfor
%!    endfor
%!    tol_red = 5e-2;
%!    tol_mod = 2e-2;
%!    tol_w = 1e-2;
%!    assert(all(err_red(:, end) < tol_red));
%!    assert(err_mod(end) < tol_mod);
%!    assert(err_w(end) < tol_w);
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 4
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!        error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      d1 = 8e-3;
%!      D1 = 12e-3;
%!      w = 5e-3;
%!      l = 47e-3;
%!      h = 5e-3;
%!      grp_id_p1 = 2;
%!      grp_id_clamp = 3;
%!      scale_def = 5e-3;
%!      mesh_size = 2.5e-3;
%!      fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "d1 = %g;\n", d1);
%!      fprintf(fd, "D1 = %g;\n", D1);
%!      fprintf(fd, "w = %g;\n", w);
%!      fprintf(fd, "l = %g;\n", l);
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(6)  = {           0,  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(7)  = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(8) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(9) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(10) = {          0, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!      fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!      fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!      fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!      fputs(fd, "Circle(5) = {7, 1, 8};\n");
%!      fputs(fd, "Circle(6) = {8, 1, 9};\n");
%!      fputs(fd, "Line(7) = {9, 10};\n");
%!      fputs(fd, "Line(8) = {10, 6};\n");
%!      fputs(fd, "Line(9) = {6, 7};\n");
%!      fputs(fd, "Curve Loop(10) = {5, 6, 7, 8, 9};\n");
%!      fputs(fd, "Curve Loop(11) = {1, 2, 3, 4};\n");
%!      fputs(fd, "Plane Surface(12) = {10, 11};\n");
%!      fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{12}; };\n");
%!      fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!      fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!      fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!      fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[7],tmp[8],tmp[9],tmp[10]};\n", grp_id_p1);
%!      fprintf(fd, "Physical Surface(\"clamp\", %d) = {tmp[5]};\n", grp_id_clamp);
%!    unwind_protect_cleanup
%!      if (fd ~= -1)
%!        fclose(fd);
%!        fd = -1;
%!      endif
%!    end_unwind_protect
%!    fprintf(stderr, "meshing ...\n");
%!    pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!    status = spawn_wait(pid);
%!    if (status ~= 0)
%!      error("gmsh failed with status %d", status);
%!    endif
%!    fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!    mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!    fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!    grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!    cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!    cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!    bearing_surf(1).group_idx = grp_idx_p1;
%!    bearing_surf(1).options.reference_pressure = 1e9;
%!    bearing_surf(1).options.mesh_size = 1e-3;
%!    bearing_surf(1).options.include_rigid_body_modes = true;
%!    bearing_surf(1).options.bearing_type = "shell";
%!    bearing_surf(1).options.matrix_type = "modal substruct total";
%!    bearing_surf(1).r = 0.5 * d1;
%!    bearing_surf(1).w = w;
%!    bearing_surf(1).X0 = [l; 0; 0];
%!    bearing_surf(1).R = eye(3);
%!    bearing_surf(1).relative_tolerance = 0;
%!    bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d1;
%!    bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!    mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!    mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_clamp, cms_opt.nodes.modal.number);
%!    mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!    mesh.elements.joints(1).nodes = cms_opt.nodes.modal.number;
%!    mesh.elements.joints(1).C = eye(6);
%!    cms_opt.inveriants = true;
%!    cms_opt.static_modes = false;
%!    cms_opt.modal_node_constraint = false;
%!    cms_opt.load_cases = "index";
%!    cms_opt.solver = "umfpack";
%!    cms_opt.refine_max_iter = int32(10);
%!    load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!    mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!    mesh.material_data.E = 210000e6;
%!    mesh.material_data.nu = 0.3;
%!    mesh.material_data.rho = 7850;
%!    opt_modes.shift_A = 0;
%!    opt_modes.refine_max_iter = int32(10);
%!    opt_modes.verbose = int32(0);
%!    opt_modes.solver = cms_opt.solver;
%!    opt_modes.active_joint_idx_eig = 1:numel(mesh.elements.joints);
%!    num_modes = [50];
%!    interfaces = {"rigid", "flexible"};
%!    num_modes_cms = int32([0, 50]);
%!    num_cases = numel(num_modes) * numel(interfaces) * numel(num_modes_cms);
%!    err_red = zeros(7, numel(num_modes));
%!    err_mod = zeros(1, numel(num_modes));
%!    err_w = zeros(1, numel(num_modes));
%!    k1 = 1;
%!    icase = int32(0);
%!    for k=1:numel(num_modes_cms)
%!      cms_opt.modes.number = num_modes_cms(k);
%!      for j=1:numel(interfaces)
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        for l=1:numel(num_modes)
%!          for i=1:numel(bearing_surf)
%!            bearing_surf(i).options.number_of_modes = min(num_modes(l), floor(numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes) * 3 / 2));
%!          endfor
%!          [mesh_l, load_case_bearing, bearing_surf_l, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          [mesh_l, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh_l, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!          [lambda, idx_lambda] = sort(diag(lambda));
%!          qred = qred(:, idx_lambda);
%!          sol_eig_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          for i=1:size(sol_eig_red.def, 3)
%!            sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!          endfor
%!          sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh_l, mat_ass, dof_map, cms_opt, bearing_surf_l);
%!          load_case_itf = fem_pre_load_case_create_empty(6);
%!          for i=1:6
%!            load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!            load_case_itf(i).loads = zeros(1, 6);
%!            load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!            switch (i)
%!              case {4, 5}
%!                load_case_itf(i).loads(i) *= w;
%!              case 6
%!                load_case_itf(i).loads(i) *= d1;
%!            endswitch
%!          endfor
%!          [Kitf, Ritf] = fem_ass_matrix(mesh_l, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!          nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!          nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!          X1 = mesh_l.nodes(mesh_l.groups.tria6(bearing_surf_l(1).group_idx).nodes, 1:3) - bearing_surf_l(1).X0.';
%!          Phi1 = atan2(X1(:, 2), X1(:, 1));
%!          p1 = k1 * sin(Phi1).^2 * bearing_surf_l(1).options.reference_pressure;
%!          Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!          z1g = comp_mat(1).bearing_surf.grid_z(:);
%!          p1red = zeros(numel(z1g), numel(Phi1g));
%!          for i=1:columns(p1red)
%!            p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!          endfor
%!          p1red = p1red(:);
%!          Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf_l(1).options.reference_pressure;
%!          Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!          qred = mat_ass.Kred \ Fred;
%!          w1red = comp_mat(1).D * qred;
%!          sol_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          mesh_post = mesh_l;
%!          mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!          load_case_post = fem_pre_load_case_create_empty(7);
%!          load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!          for i=1:6
%!            load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!            load_case_post(i).loads = load_case_itf(i).loads;
%!          endfor
%!          load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!          pn = zeros(rows(mesh_l.nodes), 1);
%!          pn(mesh_post.groups.tria6(grp_idx_p1).nodes) = p1;
%!          load_case_post(7).pressure.tria6.p = pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :));
%!          mesh_post.elements.joints.C = eye(6);
%!          mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!          dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!          [mat_ass_post.M, ...
%!           mat_ass_post.K, ...
%!           mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                            dof_map_post, ...
%!                                            [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_post);
%!          sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!          node_idx1 = mesh_l.groups.tria6(bearing_surf(1).group_idx).nodes;
%!          w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!          for i=1:numel(node_idx1)
%!            ni = [mesh_l.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!            ni /= norm(ni);
%!            w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!          endfor
%!          w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!          for i=1:columns(w1postint)
%!            for m=1:numel(Phi1g)
%!              w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                             repmat(X1(:, 3), 3, 1), ...
%!                                                                             repmat(w1post(:, i), 3, 1), ...
%!                                                                             repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                             z1g);
%!            endfor
%!          endfor
%!          if (do_plot)
%!            for m=1:columns(w1red)
%!              for i=1:numel(z1g)
%!                figure("visible", "off");
%!                hold("on");
%!                plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;1");
%!                plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;0");
%!                ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!                title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!                xlabel("Phi [deg]");
%!                ylabel("w [um]");
%!              endfor
%!            endfor
%!          endif
%!          sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!          for i=1:size(sol_eig_post.def, 3)
%!            sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!          endfor
%!          num_modes_comp = min(40, floor(numel(sol_eig_red.f) / 3));
%!          err_mod(++icase) = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!          err_w(icase) = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!          mesh_data(1).mesh = mesh_l;
%!          mesh_data(1).dof_map = dof_map;
%!          mesh_data(2).mesh = mesh_post;
%!          mesh_data(2).dof_map = dof_map_post;
%!          mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!          [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!          sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!          for i=1:size(sol_post.def, 3)
%!            err_red(i, icase) = max(max(abs(sol_post.def(1:end - 2, :, i) - sol_red.def(1:end - 2, :, i)))) / max(max(abs(sol_post.def(1:end - 2, :, i))));
%!          endfor
%!        endfor
%!        for l=1:numel(num_modes)
%!          fprintf(stderr, "%s interfaces using %d static pressure modes, %d dynamic modes:\n", interfaces{j}, num_modes(l), num_modes_cms(k));
%!          for i=1:size(sol_post.def, 3)
%!            fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i, icase));
%!          endfor
%!          fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod(icase));
%!          fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w(icase));
%!        endfor
%!      endfor
%!    endfor
%!    tol_red = 3e-2;
%!    tol_mod = 2e-2;
%!    tol_w = 2e-2;
%!    assert(all(err_red(:, end) < tol_red));
%!    assert(err_mod(end) < tol_mod);
%!    assert(err_w(end) < tol_w);
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 5
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!        error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      d1 = 8e-3;
%!      D1 = 12e-3;
%!      w = 5e-3;
%!      l = 47e-3;
%!      h = 5e-3;
%!      grp_id_p1 = 2;
%!      grp_id_clamp = 3;
%!      grp_id_itf2 = 4;
%!      scale_def = 5e-3;
%!      mesh_size = 2e-3;
%!      fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "d1 = %g;\n", d1);
%!      fprintf(fd, "D1 = %g;\n", D1);
%!      fprintf(fd, "w = %g;\n", w);
%!      fprintf(fd, "l = %g;\n", l);
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(6)  = {           0,  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(7)  = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(8) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(9) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(10) = {2 * h, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(11) = {2 * h, 4.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(12) = { h, 4.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(13) = { h, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(14) = {          0, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!      fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!      fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!      fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!      fputs(fd, "Circle(5) = {7, 1, 8};\n");
%!      fputs(fd, "Circle(6) = {8, 1, 9};\n");
%!      fputs(fd, "Line(7) = {9, 10};\n");
%!      fputs(fd, "Line(8) = {10, 11};\n");
%!      fputs(fd, "Line(9) = {11, 12};\n");
%!      fputs(fd, "Line(10) = {12, 13};\n");
%!      fputs(fd, "Line(11) = {13, 14};\n");
%!      fputs(fd, "Line(12) = {14, 6};\n");
%!      fputs(fd, "Line(13) = {6, 7};\n");
%!      fputs(fd, "Curve Loop(14) = {5, 6, 7, 8, 9, 10, 11, 12, 13};\n");
%!      fputs(fd, "Curve Loop(15) = {1, 2, 3, 4};\n");
%!      fputs(fd, "Plane Surface(16) = {14, 15};\n");
%!      fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{16}; };\n");
%!      fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!      fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!      fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!      fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[11],tmp[12],tmp[13],tmp[14]};\n", grp_id_p1);
%!      fprintf(fd, "Physical Surface(\"clamp\", %d) = {tmp[9]};\n", grp_id_clamp);
%!      fprintf(fd, "Physical Surface(\"itf2\", %d) = {tmp[6]};\n", grp_id_itf2);
%!    unwind_protect_cleanup
%!      if (fd ~= -1)
%!        fclose(fd);
%!        fd = -1;
%!      endif
%!    end_unwind_protect
%!    fprintf(stderr, "meshing ...\n");
%!    pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!    status = spawn_wait(pid);
%!    if (status ~= 0)
%!      error("gmsh failed with status %d", status);
%!    endif
%!    fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!    mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!    fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!    grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!    cms_opt.nodes.modal.number = rows(mesh.nodes) + 3;
%!    cms_opt.nodes.interfaces(1).number = rows(mesh.nodes) + 1;
%!    cms_opt.nodes.interfaces(2).number = rows(mesh.nodes) + 2;
%!    bearing_surf(1).group_idx = grp_idx_p1;
%!    bearing_surf(1).options.reference_pressure = 1e9;
%!    bearing_surf(1).options.mesh_size = 1e-3;
%!    bearing_surf(1).options.include_rigid_body_modes = true;
%!    bearing_surf(1).options.bearing_type = "shell";
%!    bearing_surf(1).options.matrix_type = "modal substruct total";
%!    bearing_surf(1).r = 0.5 * d1;
%!    bearing_surf(1).w = w;
%!    bearing_surf(1).X0 = [l; 0; 0];
%!    bearing_surf(1).R = eye(3);
%!    bearing_surf(1).relative_tolerance = 0;
%!    bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d1;
%!    bearing_surf(1).master_node_no = cms_opt.nodes.interfaces(1).number;
%!    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!    mesh.nodes(cms_opt.nodes.interfaces(1).number, 1:3) = bearing_surf(1).X0.';
%!    mesh.nodes(cms_opt.nodes.interfaces(2).number, 1:3) = [1.5 * h, 4.5 * h, 0];
%!    mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_clamp, cms_opt.nodes.modal.number);
%!    mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces(1).number);
%!    mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_itf2, cms_opt.nodes.interfaces(2).number);
%!    mesh.elements.joints(1).nodes = cms_opt.nodes.modal.number;
%!    mesh.elements.joints(1).C = eye(6);
%!    mesh.elements.joints(2).nodes = cms_opt.nodes.interfaces(2).number;
%!    mesh.elements.joints(2).C = eye(6);
%!    cms_opt.inveriants = true;
%!    cms_opt.static_modes = false;
%!    cms_opt.modal_node_constraint = false;
%!    cms_opt.load_cases = "index";
%!    cms_opt.solver = "umfpack";
%!    cms_opt.refine_max_iter = int32(10);
%!    cms_opt.verbose = int32(0);
%!    load_case = fem_pre_load_case_create_empty(6);
%!    load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!    zero_U = struct("U", cellfun(@(C) zeros(rows(C), 1), {mesh.elements.joints.C}, "UniformOutput", false));
%!    for i=1:6
%!      load_case(i).joints = zero_U;
%!      load_case(i).joints(2).U(i) = 1;
%!    endfor
%!    mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!    mesh.material_data.E = 210000e6;
%!    mesh.material_data.nu = 0.3;
%!    mesh.material_data.rho = 7850;
%!    opt_modes.shift_A = 0;
%!    opt_modes.refine_max_iter = cms_opt.refine_max_iter;
%!    opt_modes.verbose = cms_opt.verbose;
%!    opt_modes.solver = cms_opt.solver;
%!    opt_modes.active_joint_idx_eig = [1:2];
%!    opt_modes.rigid_body_modes_load_index = 1:6;
%!    num_modes = [10];
%!    interfaces = {"flexible"};
%!    num_modes_cms = int32([10]);
%!    num_cases = numel(num_modes) * numel(interfaces) * numel(num_modes_cms);
%!    err_red = zeros(7, num_cases);
%!    err_mod = zeros(1, num_cases);
%!    err_w = zeros(1, num_cases);
%!    k1 = 1;
%!    icase = int32(0);
%!    for k=1:numel(num_modes_cms)
%!      cms_opt.modes.number = num_modes_cms(k);
%!      for j=1:numel(interfaces)
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        for l=1:numel(num_modes)
%!          for i=1:numel(bearing_surf)
%!            bearing_surf(i).options.number_of_modes = min(num_modes(l), floor(numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes) * 3 / 2));
%!          endfor
%!          [mesh_l, load_case_bearing, bearing_surf_l, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          [mesh_l, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh_l, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!          [lambda, idx_lambda] = sort(diag(lambda));
%!          qred = qred(:, idx_lambda);
%!          sol_eig_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          for i=1:size(sol_eig_red.def, 3)
%!            sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!          endfor
%!          sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh_l, mat_ass, dof_map, cms_opt, bearing_surf_l);
%!          load_case_itf = fem_pre_load_case_create_empty(6);
%!          for i=1:6
%!            load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces(1).number;
%!            load_case_itf(i).loads = zeros(1, 6);
%!            load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!            switch (i)
%!              case {4, 5}
%!                load_case_itf(i).loads(i) *= w;
%!              case 6
%!                load_case_itf(i).loads(i) *= d1;
%!            endswitch
%!          endfor
%!          [Kitf, Ritf] = fem_ass_matrix(mesh_l, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!          nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!          nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!          X1 = mesh_l.nodes(mesh_l.groups.tria6(bearing_surf_l(1).group_idx).nodes, 1:3) - bearing_surf_l(1).X0.';
%!          Phi1 = atan2(X1(:, 2), X1(:, 1));
%!          p1 = k1 * sin(pi/4 + Phi1) * bearing_surf_l(1).options.reference_pressure;
%!          Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!          z1g = comp_mat(1).bearing_surf.grid_z(:);
%!          p1red = zeros(numel(z1g), numel(Phi1g));
%!          for i=1:columns(p1red)
%!            p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!          endfor
%!          p1red = p1red(:);
%!          Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf_l(1).options.reference_pressure;
%!          Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!          qred = mat_ass.Kred \ Fred;
%!          w1red = comp_mat(1).D * qred;
%!          sol_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          mesh_post = mesh_l;
%!          mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!          load_case_post = fem_pre_load_case_create_empty(7);
%!          load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!          for i=1:6
%!            load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!            load_case_post(i).loads = load_case_itf(i).loads;
%!          endfor
%!          load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!          pn = zeros(rows(mesh_l.nodes), 1);
%!          pn(mesh_post.groups.tria6(grp_idx_p1).nodes) = p1;
%!          load_case_post(7).pressure.tria6.p = pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :));
%!          mesh_post.elements.joints.C = eye(6);
%!          mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!          dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!          [mat_ass_post.M, ...
%!           mat_ass_post.K, ...
%!           mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                            dof_map_post, ...
%!                                            [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_post);
%!          sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!          node_idx1 = mesh_l.groups.tria6(bearing_surf(1).group_idx).nodes;
%!          w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!          for i=1:numel(node_idx1)
%!            ni = [mesh_l.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!            ni /= norm(ni);
%!            w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!          endfor
%!          w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!          for i=1:columns(w1postint)
%!            for m=1:numel(Phi1g)
%!              w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                             repmat(X1(:, 3), 3, 1), ...
%!                                                                             repmat(w1post(:, i), 3, 1), ...
%!                                                                             repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                             z1g);
%!            endfor
%!          endfor
%!          if (do_plot)
%!            for m=1:columns(w1red)
%!              for i=1:numel(z1g)
%!                figure("visible", "off");
%!                hold("on");
%!                plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;1");
%!                plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;0");
%!                ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!                title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!                xlabel("Phi [deg]");
%!                ylabel("w [um]");
%!              endfor
%!            endfor
%!          endif
%!          sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!          for i=1:size(sol_eig_post.def, 3)
%!            sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!          endfor
%!          num_modes_comp = min(40, floor(numel(sol_eig_red.f) / 3));
%!          err_mod(++icase) = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!          err_w(icase) = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!          mesh_data(1).mesh = mesh_l;
%!          mesh_data(1).dof_map = dof_map;
%!          mesh_data(2).mesh = mesh_post;
%!          mesh_data(2).dof_map = dof_map_post;
%!          mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!          [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!          sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!          for i=1:size(sol_post.def, 3)
%!            err_red(i, icase) = norm(sol_post.def(:, :, i) - sol_red.def(:, :, i)) / norm(sol_post.def(:, :, i));
%!          endfor
%!        endfor
%!        for l=1:numel(num_modes)
%!          fprintf(stderr, "%s interfaces using %d static pressure modes, %d dynamic modes:\n", interfaces{j}, num_modes(l), num_modes_cms(k));
%!          for i=1:size(sol_post.def, 3)
%!            fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i, icase));
%!          endfor
%!          fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod(icase));
%!          fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w(icase));
%!        endfor
%!      endfor
%!    endfor
%!    tol_red = 1e-2;
%!    tol_mod = 1e-2;
%!    tol_w = 1e-2;
%!    assert(all(err_red(:, end) < tol_red));
%!    assert(err_mod(end) < tol_mod);
%!    assert(err_w(end) < tol_w);
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 6
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    interfaces = {"rigid", "flexible"};
%!    tol_red = [5e-2, 2e-2];
%!    num_modes_cms = int32([0, 10]);
%!    num_modes = [50, 30];
%!    for k=1:numel(num_modes_cms)
%!      for j=1:numel(interfaces)
%!        clear Fred Ritf bearing_surf cms_opt comp_mat dof_map_comb dof_map_post err_red
%!        clear grp_id_p1 grp_id_p2 grp_idx_p1 grp_idx_p2 load_case load_case_bearing load_case_itf
%!        clear load_case_post mat_ass_post mat_ass_press mesh mesh_comb mesh_data mesh_post mesh_size
%!        clear opt_modes p1 p1red p2 p2red pid qred sol_comb sol_eig sol_eig_cms sol_post sol_red
%!        unwind_protect
%!          [fd, msg] = fopen([filename, ".geo"], "w");
%!          if (fd == -1)
%!            error("failed to open file \"%s.geo\"", filename);
%!          endif
%!          d = 14e-3;
%!          D = 19.5e-3;
%!          w = 5e-3;
%!          l = 47e-3;
%!          h = 5e-3;
%!          grp_id_p1 = 2;
%!          grp_id_p2 = 3;
%!          p1 = 1;
%!          p2 = 2;
%!          scale_def = 5e-3;
%!          mesh_size = 8e-3;
%!          fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!          fprintf(fd, "d = %g;\n", d);
%!          fprintf(fd, "D = %g;\n", D);
%!          fprintf(fd, "w = %g;\n", w);
%!          fprintf(fd, "l = %g;\n", l);
%!          fprintf(fd, "h = %g;\n", h);
%!          fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!          fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!          fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!          fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!          fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!          fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!          fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!          fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!          fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!          fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!          fputs(fd, "Line(11) = {13, 14};\n");
%!          fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!          fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!          fputs(fd, "Line(14) = {16, 11};\n");
%!          fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!          fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!          fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!          fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!          fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!          fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!          fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!          fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!          fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!          fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!        unwind_protect_cleanup
%!          if (fd ~= -1)
%!            fclose(fd);
%!            fd = -1;
%!          endif
%!        end_unwind_protect
%!        fprintf(stderr, "meshing ...\n");
%!        pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!        status = spawn_wait(pid);
%!        if (status ~= 0)
%!          error("gmsh failed with status %d", status);
%!        endif
%!        fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!        mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!        fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!        grp_idx_p1 = find([[mesh.groups.tria10].id] == grp_id_p1);
%!        grp_idx_p2 = find([[mesh.groups.tria10].id] == grp_id_p2);
%!        cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!        cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!        bearing_surf(1).group_idx = grp_idx_p1;
%!        bearing_surf(1).options.reference_pressure = 1e9;
%!        bearing_surf(1).options.mesh_size = 1e-3;
%!        bearing_surf(1).options.include_rigid_body_modes = true;
%!        bearing_surf(1).options.bearing_type = "shell";
%!        bearing_surf(1).options.matrix_type = "modal substruct total";
%!        bearing_surf(1).r = 0.5 * d;
%!        bearing_surf(1).w = w;
%!        bearing_surf(1).X0 = [l; 0; 0];
%!        bearing_surf(1).R = eye(3);
%!        bearing_surf(1).relative_tolerance = 0;
%!        bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!        bearing_surf(1).options.number_of_modes = num_modes(j);
%!        bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!        bearing_surf(2).group_idx = grp_idx_p2;
%!        bearing_surf(2).options.reference_pressure = 1e9;
%!        bearing_surf(2).options.mesh_size = 1e-3;
%!        bearing_surf(2).options.include_rigid_body_modes = false;
%!        bearing_surf(2).options.bearing_type = "shell";
%!        bearing_surf(2).options.matrix_type = "modal substruct total";
%!        bearing_surf(2).r = 0.5 * d;
%!        bearing_surf(2).w = w;
%!        bearing_surf(2).X0 = [0; 0; 0];
%!        bearing_surf(2).R = eye(3);
%!        bearing_surf(2).relative_tolerance = 0;
%!        bearing_surf(2).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!        bearing_surf(2).options.number_of_modes = num_modes(j);
%!        bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!        mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!        mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!        mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number, "tria10");
%!        mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number, "tria10");
%!        cms_opt.inveriants = true;
%!        cms_opt.modes.number = num_modes_cms(k);
%!        cms_opt.static_modes = false;
%!        cms_opt.modal_node_constraint = false;
%!        cms_opt.load_cases = "index";
%!        cms_opt.refine_max_iter = int32(10);
%!        load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!        mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!        mesh.material_data.E = 210000e6;
%!        mesh.material_data.nu = 0.3;
%!        mesh.material_data.rho = 7850;
%!        opt_modes.shift_A = 1e-6;
%!        opt_modes.refine_max_iter = int32(10);
%!        opt_modes.verbose = int32(0);
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        opt_modes.solver = "umfpack";
%!        opt_modes.elem_type = "tria10";
%!        cms_opt.solver = opt_modes.solver;
%!        [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!        [mesh, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!        assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!        assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!        comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%!        load_case_itf = fem_pre_load_case_create_empty(6);
%!        for i=1:6
%!          load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!          load_case_itf(i).loads = zeros(1, 6);
%!          load_case_itf(i).loads(i) = 1;
%!        endfor
%!        [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!        nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!        nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!        nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!        nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!        p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!        p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!        p1red2 = zeros((nx1 - 1) * nz1, 1);
%!        p2red2 = zeros((nx2 - 1) * nz2, 1);
%!        for i=1:nx1 - 1
%!          p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!        endfor
%!        for i=1:nx2 - 1
%!          p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!        endfor
%!        Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!                comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!        Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!        qred = mat_ass.Kred \ Fred;
%!        sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!        mesh_post = mesh;
%!        mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!        load_case_post = fem_pre_load_case_create_empty(10);
%!        load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!        for i=1:6
%!          load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!          load_case_post(i).loads = zeros(1, 6);
%!          load_case_post(i).loads(i) = 1;
%!        endfor
%!        x1 = mesh.nodes(:, 1)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!        y1 = mesh.nodes(:, 2)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!        x2 = mesh.nodes(:, 1)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!        y2 = mesh.nodes(:, 2)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!        Phi1 = atan2(y1, x1);
%!        Phi2 = atan2(y2, x2);
%!        load_case_post(7).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!        load_case_post(7).pressure.tria10.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria10(grp_idx_p1).elements), 10);
%!        load_case_post(8).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!        load_case_post(8).pressure.tria10.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!        load_case_post(9).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p2).elements, :);
%!        load_case_post(9).pressure.tria10.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria10(grp_idx_p2).elements), 10);
%!        load_case_post(10).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p2).elements, :);
%!        load_case_post(10).pressure.tria10.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!        mesh_post.elements.joints.C = eye(6);
%!        mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!        dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!        [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!        sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, opt_modes);
%!        mesh_data(1).mesh = mesh;
%!        mesh_data(1).dof_map = dof_map;
%!        mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!        mesh_data(2).mesh = mesh_post;
%!        mesh_data(2).dof_map = dof_map_post;
%!        [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!        sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!        sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!        sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!        for i=1:size(sol_comb.def, 3)
%!          sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!        endfor
%!        err_red = zeros(1, size(sol_post.def, 3));
%!        for i=1:size(sol_post.def, 3)
%!          err_red(i) = max(max(abs(sol_post.def(1:end - 2, :, i) - sol_red.def(1:end - 2, :, i)))) / max(max(abs(sol_post.def(1:end - 2, :, i))));
%!        endfor
%!        for i=1:size(sol_post.def, 3)
%!          fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!        endfor
%!        assert(all(err_red < tol_red(j)));
%!      endfor
%!    endfor
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 7
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    f_enable_constraint = [false, true];
%!    interfaces = {"flexible","rigid"};
%!    num_modes = int32([0, 10]);
%!    for l=1:numel(num_modes)
%!      for k=1:numel(f_enable_constraint)
%!        for j=1:numel(interfaces)
%!          clear bearing_surf cms_opt comp_mat dof_map grp_id_clamp grp_id_p1 grp_id_p2;
%!          clear load_case load_case_bearing mat_ass mat_ass_press mat_info mesh mesh_info mesh_size
%!          clear opt_modes sol_eig sol_eig_cms sol_stat;
%!          unwind_protect
%!            [fd, msg] = fopen([filename, ".geo"], "w");
%!            if (fd == -1)
%!              error("failed to open file \"%s.geo\"", filename);
%!            endif
%!            ri = 8e-3;
%!            ro = 10e-3;
%!            h = 12e-3;
%!            c = 2e-3;
%!            b = h - 2 * c;
%!            scale_def = 5e-3;
%!            mesh_size = 6e-3;
%!            fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!            fprintf(fd, "ri = %g;\n", ri);
%!            fprintf(fd, "ro = %g;\n", ro);
%!            fprintf(fd, "h = %g;\n", h);
%!            fprintf(fd, "c = %g;\n", c);
%!            fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!            fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!            fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!            fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!            fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!            fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!            fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!            fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!            fputs(fd, "Line(1) = {1,2};\n");
%!            fputs(fd, "Line(2) = {2,3};\n");
%!            fputs(fd, "Line(3) = {3,4};\n");
%!            fputs(fd, "Line(4) = {4,5};\n");
%!            fputs(fd, "Line(5) = {5,6};\n");
%!            fputs(fd, "Line(6) = {6,7};\n");
%!            fputs(fd, "Line(7) = {7,8};\n");
%!            fputs(fd, "Line(8) = {8,1};\n");
%!            fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!            fputs(fd, "Plane Surface(6) = {5};\n");
%!            fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%!            fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!            fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!            fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!            fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!            fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!            fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!          unwind_protect_cleanup
%!            if (fd ~= -1)
%!              fclose(fd);
%!              fd = -1;
%!            endif
%!          end_unwind_protect
%!          fprintf(stderr, "meshing ...\n");
%!          pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!          status = spawn_wait(pid);
%!          if (status ~= 0)
%!            error("gmsh failed with status %d", status);
%!          endif
%!          fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!          mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!          fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!          cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!          cms_opt.solver = "mldivide";
%!          switch (interfaces{j})
%!            case "flexible"
%!              cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!          endswitch
%!          grp_id_clamp = find([[mesh.groups.tria10].id] == 1);
%!          grp_id_p1 = find([[mesh.groups.tria10].id] == 3);
%!          grp_id_p2 = find([[mesh.groups.tria10].id] == 2);
%!          bearing_surf(1).group_idx = grp_id_p1;
%!          bearing_surf(1).options.reference_pressure = 1e9;
%!          bearing_surf(1).options.mesh_size = 20e-3;
%!          bearing_surf(1).options.include_rigid_body_modes = false;
%!          bearing_surf(1).options.bearing_type = "shell";
%!          bearing_surf(1).options.matrix_type = "modal substruct total";
%!          bearing_surf(1).r = ri;
%!          bearing_surf(1).w = b;
%!          bearing_surf(1).X0 = [0; 0; b/2 + c];
%!          bearing_surf(1).R = eye(3);
%!          bearing_surf(1).relative_tolerance = 0;
%!          bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%!          bearing_surf(1).options.number_of_modes = 10;
%!          bearing_surf(2).group_idx = grp_id_p2;
%!          bearing_surf(2).options.reference_pressure = 1e9;
%!          bearing_surf(2).options.mesh_size = 20e-3;
%!          bearing_surf(2).options.include_rigid_body_modes = true;
%!          bearing_surf(2).options.bearing_type = "journal";
%!          bearing_surf(2).options.matrix_type = "modal substruct total";
%!          bearing_surf(2).r = ro;
%!          bearing_surf(2).w = b;
%!          bearing_surf(2).X0 = [0; 0; b/2 + c];
%!          bearing_surf(2).R = eye(3);
%!          bearing_surf(2).relative_tolerance = 0;
%!          bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%!          bearing_surf(2).options.number_of_modes = 12;
%!          switch (interfaces{j})
%!            case "flexible"
%!              bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!              bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!              for i=1:numel(bearing_surf)
%!                mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!              endfor
%!              for i=1:numel(bearing_surf)
%!                mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no, "tria10");
%!              endfor
%!            otherwise
%!              mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!          endswitch
%!          cms_opt.inveriants = true;
%!          cms_opt.modes.number = num_modes(l);
%!          cms_opt.static_modes = false;
%!          cms_opt.modal_node_constraint = false;
%!          cms_opt.load_cases = "index";
%!          cms_opt.refine_max_iter = int32(10);
%!          load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!          switch (interfaces{j})
%!            case "flexible"
%!            otherwise
%!              load_case(1).locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%!          endswitch
%!          if (f_enable_constraint(k))
%!            load_case(1).locked_dof(mesh.groups.tria10(grp_id_clamp).nodes, :) = true;
%!          endif
%!          mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!          mesh.material_data.E = 210000e6;
%!          mesh.material_data.nu = 0.3;
%!          mesh.material_data.rho = 7850;
%!          if (f_enable_constraint(k))
%!            opt_modes.shift_A = 0;
%!          else
%!            opt_modes.shift_A = 1e-6;
%!          endif
%!          opt_modes.refine_max_iter = int32(10);
%!          opt_modes.verbose = int32(0);
%!          opt_modes.solver = cms_opt.solver;
%!          opt_modes.rigid_body_modes = interfaces{j};
%!          opt_modes.elem_type = "tria10";
%!          [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          dof_map = fem_ass_dof_map(mesh, load_case);
%!          [mat_ass.K, ...
%!           mat_ass.M, ...
%!           mat_ass.R, ...
%!           mat_info, ...
%!           mesh_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_STIFFNESS, ...
%!                                        FEM_MAT_MASS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case_bearing);
%!          sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!          [mesh, ...
%!           mat_ass, ...
%!           dof_map, ...
%!           sol_eig_cms, ...
%!           cms_optp] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_optp, bearing_surf);
%!        endfor
%!      endfor
%!    endfor
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 8
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!        error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      d1 = 8e-3;
%!      D1 = 12e-3;
%!      d2 = 14e-3;
%!      D2 = 19.5e-3;
%!      w = 5e-3;
%!      l = 47e-3;
%!      h = 5e-3;
%!      grp_id_p1 = 2;
%!      grp_id_p2 = 3;
%!      scale_def = 5e-3;
%!      mesh_size = 2.5e-3;
%!      fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "d1 = %g;\n", d1);
%!      fprintf(fd, "D1 = %g;\n", D1);
%!      fprintf(fd, "d2 = %g;\n", d2);
%!      fprintf(fd, "D2 = %g;\n", D2);
%!      fprintf(fd, "w = %g;\n", w);
%!      fprintf(fd, "l = %g;\n", l);
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(6)  = {         0.0,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(7)  = {         0.0,  0.5 * d2, -0.5 * w};\n");
%!      fputs(fd, "Point(8)  = {    0.5 * d2,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(9)  = {         0.0, -0.5 * d2, -0.5 * w};\n");
%!      fputs(fd, "Point(10) = {   -0.5 * d2,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(11) = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(12) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(13) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(14) = {Sqrt((D2/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(15) = {   -0.5 * D2,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(16) = {Sqrt((D2/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!      fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!      fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!      fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!      fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!      fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!      fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!      fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!      fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!      fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!      fputs(fd, "Line(11) = {13, 14};\n");
%!      fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!      fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!      fputs(fd, "Line(14) = {16, 11};\n");
%!      fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!      fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!      fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!      fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!      fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!      fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!      fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!      fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!      fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!      fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!    unwind_protect_cleanup
%!      if (fd ~= -1)
%!        fclose(fd);
%!        fd = -1;
%!      endif
%!    end_unwind_protect
%!    fprintf(stderr, "meshing ...\n");
%!    pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!    status = spawn_wait(pid);
%!    if (status ~= 0)
%!      error("gmsh failed with status %d", status);
%!    endif
%!    fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!    mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!    fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!    grp_idx_p1 = find([[mesh.groups.tria10].id] == grp_id_p1);
%!    grp_idx_p2 = find([[mesh.groups.tria10].id] == grp_id_p2);
%!    cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!    cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!    bearing_surf(1).group_idx = grp_idx_p1;
%!    bearing_surf(1).options.reference_pressure = 1e9;
%!    bearing_surf(1).options.mesh_size = 1e-3;
%!    bearing_surf(1).options.include_rigid_body_modes = true;
%!    bearing_surf(1).options.bearing_type = "shell";
%!    bearing_surf(1).options.matrix_type = "modal substruct total";
%!    bearing_surf(1).r = 0.5 * d1;
%!    bearing_surf(1).w = w;
%!    bearing_surf(1).X0 = [l; 0; 0];
%!    bearing_surf(1).R = eye(3);
%!    bearing_surf(1).relative_tolerance = 0;
%!    bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d1;
%!    bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!    bearing_surf(2).group_idx = grp_idx_p2;
%!    bearing_surf(2).options.reference_pressure = 1e9;
%!    bearing_surf(2).options.mesh_size = 1e-3;
%!    bearing_surf(2).options.include_rigid_body_modes = false;
%!    bearing_surf(2).options.bearing_type = "shell";
%!    bearing_surf(2).options.matrix_type = "modal substruct total";
%!    bearing_surf(2).r = 0.5 * d2;
%!    bearing_surf(2).w = w;
%!    bearing_surf(2).X0 = [0; 0; 0];
%!    bearing_surf(2).R = eye(3);
%!    bearing_surf(2).relative_tolerance = 0;
%!    bearing_surf(2).absolute_tolerance = sqrt(eps) * 0.5 * d2;
%!    bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!    mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!    mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number, "tria10");
%!    mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number, "tria10");
%!    cms_opt.inveriants = true;
%!    cms_opt.static_modes = false;
%!    cms_opt.modal_node_constraint = false;
%!    cms_opt.load_cases = "index";
%!    cms_opt.solver = "umfpack";
%!    cms_opt.refine_max_iter = int32(10);
%!    load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!    mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!    mesh.material_data.E = 210000e6;
%!    mesh.material_data.nu = 0.3;
%!    mesh.material_data.rho = 7850;
%!    opt_modes.shift_A = 1e-6;
%!    opt_modes.refine_max_iter = int32(10);
%!    opt_modes.verbose = int32(0);
%!    opt_modes.solver = cms_opt.solver;
%!    opt_modes.elem_type = "tria10";
%!    num_modes = [0, 30];
%!    interfaces = {"rigid", "flexible"};
%!    num_modes_cms = int32([0, 10]);
%!    k1 = 1;
%!    k2 = 0;
%!    num_cases = numel(num_modes) * numel(interfaces) * numel(num_modes_cms);
%!    err_red = zeros(7, num_cases);
%!    err_mod = zeros(1, num_cases);
%!    err_w = zeros(1, num_cases);
%!    icase = int32(0);
%!    for k=1:numel(num_modes_cms)
%!      cms_opt.modes.number = num_modes_cms(k);
%!      for j=1:numel(interfaces)
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        for l=1:numel(num_modes)
%!          for i=1:numel(bearing_surf)
%!            bearing_surf(i).options.number_of_modes = min(num_modes(l), floor(numel(mesh.groups.tria10(bearing_surf(i).group_idx).nodes) * 3 / 2 - 6));
%!          endfor
%!          [mesh_l, load_case_bearing, bearing_surf_l, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          [mesh_l, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh_l, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!          [lambda, idx_lambda] = sort(diag(lambda));
%!          qred = qred(:, idx_lambda);
%!          sol_eig_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          for i=1:size(sol_eig_red.def, 3)
%!            sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!          endfor
%!          sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh_l, mat_ass, dof_map, cms_opt, bearing_surf_l);
%!          load_case_itf = fem_pre_load_case_create_empty(6);
%!          for i=1:6
%!            load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!            load_case_itf(i).loads = zeros(1, 6);
%!            load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!            switch (i)
%!              case {4, 5}
%!                load_case_itf(i).loads(i) *= w;
%!              case 6
%!                load_case_itf(i).loads(i) *= d1;
%!            endswitch
%!          endfor
%!          [~, Ritf] = fem_ass_matrix(mesh_l, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!          nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!          nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!          nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!          nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!          X1 = mesh_l.nodes(mesh_l.groups.tria10(bearing_surf_l(1).group_idx).nodes, 1:3) - bearing_surf_l(1).X0.';
%!          Phi1 = atan2(X1(:, 2), X1(:, 1));
%!          p1 = k1 * sin(Phi1).^2 * bearing_surf_l(1).options.reference_pressure;
%!          Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!          z1g = comp_mat(1).bearing_surf.grid_z(:);
%!          p1red = zeros(numel(z1g), numel(Phi1g));
%!          for i=1:columns(p1red)
%!            p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!          endfor
%!          p1red = p1red(:);
%!          X2 = mesh_l.nodes(mesh_l.groups.tria10(bearing_surf_l(2).group_idx).nodes, 1:3) - bearing_surf_l(2).X0.';
%!          Phi2 = atan2(X2(:, 2), X2(:, 1));
%!          p2 = k2 * cos(Phi2) .* sin(pi * X2(:, 3) / (max(X2(:, 3)) - min(X2(:, 3)))) * bearing_surf_l(2).options.reference_pressure;
%!          Phi2g = comp_mat(2).bearing_surf.grid_x(:) / (0.5 * comp_mat(2).bearing_dimensions.bearing_diameter);
%!          z2g = comp_mat(2).bearing_surf.grid_z(:);
%!          p2red = zeros(numel(z2g), numel(Phi2g));
%!          for i=1:columns(p2red)
%!            p2red(:, i) = griddata([Phi2 - 2 * pi; Phi2; Phi2 + 2 * pi], repmat(X2(:, 3), 3, 1), repmat(p2, 3, 1), repmat(Phi2g(i), rows(z2g), 1), z2g);
%!          endfor
%!          p2red = p2red(:);
%!          Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf_l(1).options.reference_pressure + comp_mat(2).E(:, 1:end - nz2) * p2red(1:end -  nz2) / bearing_surf_l(2).options.reference_pressure;
%!          Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!          qred = mat_ass.Kred \ Fred;
%!          w1red = comp_mat(1).D * qred;
%!          w2red = comp_mat(2).D * qred;
%!          sol_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          mesh_post = mesh_l;
%!          mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!          load_case_post = fem_pre_load_case_create_empty(7);
%!          load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!          for i=1:6
%!            load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!            load_case_post(i).loads = load_case_itf(i).loads;
%!          endfor
%!          load_case_post(7).pressure.tria10.elements = [mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!                                                       mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p2).elements, :)];
%!          pn = zeros(rows(mesh_l.nodes), 1);
%!          pn(mesh_post.groups.tria10(grp_idx_p1).nodes) = p1;
%!          pn(mesh_post.groups.tria10(grp_idx_p2).nodes) = p2;
%!          load_case_post(7).pressure.tria10.p = [pn(mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :));
%!                                                pn(mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p2).elements, :))];
%!          mesh_post.elements.joints.C = eye(6);
%!          mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!          dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!          [mat_ass_post.M, ...
%!           mat_ass_post.K, ...
%!           mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                            dof_map_post, ...
%!                                            [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_post);
%!          sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!          node_idx1 = mesh_l.groups.tria10(bearing_surf(1).group_idx).nodes;
%!          w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!          for i=1:numel(node_idx1)
%!            ni = [mesh_l.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!            ni /= norm(ni);
%!            w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!          endfor
%!          w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!          for i=1:columns(w1postint)
%!            for m=1:numel(Phi1g)
%!              w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                             repmat(X1(:, 3), 3, 1), ...
%!                                                                             repmat(w1post(:, i), 3, 1), ...
%!                                                                             repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                             z1g);
%!            endfor
%!          endfor
%!          if (do_plot)
%!            for m=1:columns(w1red)
%!              for i=1:numel(z1g)
%!                figure("visible", "off");
%!                hold("on");
%!                plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;1");
%!                plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;0");
%!                ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!                title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!                xlabel("Phi [deg]");
%!                ylabel("w [um]");
%!              endfor
%!            endfor
%!          endif
%!          sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!          for i=1:size(sol_eig_post.def, 3)
%!            sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!          endfor
%!          num_modes_comp = min(10, floor(numel(sol_eig_red.f) / 3));
%!          err_mod(++icase) = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!          err_w(icase) = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!          mesh_data(1).mesh = mesh_l;
%!          mesh_data(1).dof_map = dof_map;
%!          mesh_data(2).mesh = mesh_post;
%!          mesh_data(2).dof_map = dof_map_post;
%!          mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!          [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!          sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!          for i=1:size(sol_post.def, 3)
%!            err_red(i, icase) = max(max(abs(sol_post.def(1:end - 2, 1:3, i) - sol_red.def(1:end - 2, 1:3, i)))) / max(max(abs(sol_post.def(1:end - 2, 1:3, i))));
%!          endfor
%!        endfor
%!        for l=1:numel(num_modes)
%!          fprintf(stderr, "%s interfaces using %d static pressure modes, %d dynamic modes:\n", interfaces{j}, num_modes(l), num_modes_cms(k));
%!          for i=1:size(sol_post.def, 3)
%!            fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i, icase));
%!          endfor
%!          fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod(icase));
%!          fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w(icase));
%!        endfor
%!      endfor
%!    endfor
%!    tol_red = 5e-2;
%!    tol_mod = 2e-2;
%!    tol_w = 1e-2;
%!    assert(all(err_red(:, end) < tol_red));
%!    assert(err_mod(end) < tol_mod);
%!    assert(err_w(end) < tol_w);
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 9
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!        error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      d1 = 8e-3;
%!      D1 = 12e-3;
%!      w = 5e-3;
%!      l = 47e-3;
%!      h = 5e-3;
%!      grp_id_p1 = 2;
%!      grp_id_clamp = 3;
%!      scale_def = 5e-3;
%!      mesh_size = 5e-3;
%!      fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "d1 = %g;\n", d1);
%!      fprintf(fd, "D1 = %g;\n", D1);
%!      fprintf(fd, "w = %g;\n", w);
%!      fprintf(fd, "l = %g;\n", l);
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(6)  = {           0,  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(7)  = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(8) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(9) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(10) = {          0, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!      fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!      fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!      fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!      fputs(fd, "Circle(5) = {7, 1, 8};\n");
%!      fputs(fd, "Circle(6) = {8, 1, 9};\n");
%!      fputs(fd, "Line(7) = {9, 10};\n");
%!      fputs(fd, "Line(8) = {10, 6};\n");
%!      fputs(fd, "Line(9) = {6, 7};\n");
%!      fputs(fd, "Curve Loop(10) = {5, 6, 7, 8, 9};\n");
%!      fputs(fd, "Curve Loop(11) = {1, 2, 3, 4};\n");
%!      fputs(fd, "Plane Surface(12) = {10, 11};\n");
%!      fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{12}; };\n");
%!      fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!      fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!      fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!      fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[7],tmp[8],tmp[9],tmp[10]};\n", grp_id_p1);
%!      fprintf(fd, "Physical Surface(\"clamp\", %d) = {tmp[5]};\n", grp_id_clamp);
%!    unwind_protect_cleanup
%!      if (fd ~= -1)
%!        fclose(fd);
%!        fd = -1;
%!      endif
%!    end_unwind_protect
%!    fprintf(stderr, "meshing ...\n");
%!    pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!    status = spawn_wait(pid);
%!    if (status ~= 0)
%!      error("gmsh failed with status %d", status);
%!    endif
%!    fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!    mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!    fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!    grp_idx_p1 = find([[mesh.groups.tria10].id] == grp_id_p1);
%!    cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!    cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!    bearing_surf(1).group_idx = grp_idx_p1;
%!    bearing_surf(1).options.reference_pressure = 1e9;
%!    bearing_surf(1).options.mesh_size = 1e-3;
%!    bearing_surf(1).options.include_rigid_body_modes = true;
%!    bearing_surf(1).options.bearing_type = "shell";
%!    bearing_surf(1).options.matrix_type = "modal substruct total";
%!    bearing_surf(1).r = 0.5 * d1;
%!    bearing_surf(1).w = w;
%!    bearing_surf(1).X0 = [l; 0; 0];
%!    bearing_surf(1).R = eye(3);
%!    bearing_surf(1).relative_tolerance = 0;
%!    bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d1;
%!    bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!    mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!    mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_clamp, cms_opt.nodes.modal.number, "tria10");
%!    mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number, "tria10");
%!    mesh.elements.joints(1).nodes = cms_opt.nodes.modal.number;
%!    mesh.elements.joints(1).C = eye(6);
%!    cms_opt.inveriants = true;
%!    cms_opt.static_modes = false;
%!    cms_opt.modal_node_constraint = false;
%!    cms_opt.load_cases = "index";
%!    cms_opt.solver = "umfpack";
%!    cms_opt.refine_max_iter = int32(10);
%!    load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!    mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!    mesh.material_data.E = 210000e6;
%!    mesh.material_data.nu = 0.3;
%!    mesh.material_data.rho = 7850;
%!    opt_modes.shift_A = 0;
%!    opt_modes.refine_max_iter = int32(10);
%!    opt_modes.verbose = int32(0);
%!    opt_modes.solver = cms_opt.solver;
%!    opt_modes.active_joint_idx_eig = 1:numel(mesh.elements.joints);
%!    opt_modes.elem_type = "tria10";
%!    num_modes = [70];
%!    interfaces = {"rigid", "flexible"};
%!    num_modes_cms = int32([0, 50]);
%!    num_cases = numel(num_modes) * numel(interfaces) * numel(num_modes_cms);
%!    err_red = zeros(7, numel(num_modes));
%!    err_mod = zeros(1, numel(num_modes));
%!    err_w = zeros(1, numel(num_modes));
%!    k1 = 1;
%!    icase = int32(0);
%!    for k=1:numel(num_modes_cms)
%!      cms_opt.modes.number = num_modes_cms(k);
%!      for j=1:numel(interfaces)
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        for l=1:numel(num_modes)
%!          for i=1:numel(bearing_surf)
%!            bearing_surf(i).options.number_of_modes = min(num_modes(l), floor(numel(mesh.groups.tria10(bearing_surf(i).group_idx).nodes) * 3 / 2));
%!          endfor
%!          [mesh_l, load_case_bearing, bearing_surf_l, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          [mesh_l, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh_l, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!          [lambda, idx_lambda] = sort(diag(lambda));
%!          qred = qred(:, idx_lambda);
%!          sol_eig_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          for i=1:size(sol_eig_red.def, 3)
%!            sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!          endfor
%!          sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh_l, mat_ass, dof_map, cms_opt, bearing_surf_l);
%!          load_case_itf = fem_pre_load_case_create_empty(6);
%!          for i=1:6
%!            load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!            load_case_itf(i).loads = zeros(1, 6);
%!            load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!            switch (i)
%!              case {4, 5}
%!                load_case_itf(i).loads(i) *= w;
%!              case 6
%!                load_case_itf(i).loads(i) *= d1;
%!            endswitch
%!          endfor
%!          [Kitf, Ritf] = fem_ass_matrix(mesh_l, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!          nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!          nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!          X1 = mesh_l.nodes(mesh_l.groups.tria10(bearing_surf_l(1).group_idx).nodes, 1:3) - bearing_surf_l(1).X0.';
%!          Phi1 = atan2(X1(:, 2), X1(:, 1));
%!          p1 = k1 * sin(Phi1).^2 * bearing_surf_l(1).options.reference_pressure;
%!          Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!          z1g = comp_mat(1).bearing_surf.grid_z(:);
%!          p1red = zeros(numel(z1g), numel(Phi1g));
%!          for i=1:columns(p1red)
%!            p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!          endfor
%!          p1red = p1red(:);
%!          Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf_l(1).options.reference_pressure;
%!          Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!          qred = mat_ass.Kred \ Fred;
%!          w1red = comp_mat(1).D * qred;
%!          sol_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          mesh_post = mesh_l;
%!          mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!          load_case_post = fem_pre_load_case_create_empty(7);
%!          load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!          for i=1:6
%!            load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!            load_case_post(i).loads = load_case_itf(i).loads;
%!          endfor
%!          load_case_post(7).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!          pn = zeros(rows(mesh_l.nodes), 1);
%!          pn(mesh_post.groups.tria10(grp_idx_p1).nodes) = p1;
%!          load_case_post(7).pressure.tria10.p = pn(mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :));
%!          mesh_post.elements.joints.C = eye(6);
%!          mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!          dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!          [mat_ass_post.M, ...
%!           mat_ass_post.K, ...
%!           mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                            dof_map_post, ...
%!                                            [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_post);
%!          sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!          node_idx1 = mesh_l.groups.tria10(bearing_surf(1).group_idx).nodes;
%!          w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!          for i=1:numel(node_idx1)
%!            ni = [mesh_l.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!            ni /= norm(ni);
%!            w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!          endfor
%!          w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!          for i=1:columns(w1postint)
%!            for m=1:numel(Phi1g)
%!              w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                             repmat(X1(:, 3), 3, 1), ...
%!                                                                             repmat(w1post(:, i), 3, 1), ...
%!                                                                             repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                             z1g);
%!            endfor
%!          endfor
%!          if (do_plot)
%!            for m=1:columns(w1red)
%!              for i=1:numel(z1g)
%!                figure("visible", "off");
%!                hold("on");
%!                plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;1");
%!                plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;0");
%!                ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!                title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!                xlabel("Phi [deg]");
%!                ylabel("w [um]");
%!              endfor
%!            endfor
%!          endif
%!          sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!          for i=1:size(sol_eig_post.def, 3)
%!            sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!          endfor
%!          num_modes_comp = min(40, floor(numel(sol_eig_red.f) / 3));
%!          err_mod(++icase) = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!          err_w(icase) = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!          mesh_data(1).mesh = mesh_l;
%!          mesh_data(1).dof_map = dof_map;
%!          mesh_data(2).mesh = mesh_post;
%!          mesh_data(2).dof_map = dof_map_post;
%!          mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!          [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!          sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!          for i=1:size(sol_post.def, 3)
%!            err_red(i, icase) = max(max(abs(sol_post.def(1:end - 2, :, i) - sol_red.def(1:end - 2, :, i)))) / max(max(abs(sol_post.def(1:end - 2, :, i))));
%!          endfor
%!        endfor
%!        for l=1:numel(num_modes)
%!          fprintf(stderr, "%s interfaces using %d static pressure modes, %d dynamic modes:\n", interfaces{j}, num_modes(l), num_modes_cms(k));
%!          for i=1:size(sol_post.def, 3)
%!            fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i, icase));
%!          endfor
%!          fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod(icase));
%!          fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w(icase));
%!        endfor
%!      endfor
%!    endfor
%!    tol_red = 3e-2;
%!    tol_mod = 2e-2;
%!    tol_w = 2e-2;
%!    assert(all(err_red(:, end) < tol_red));
%!    assert(err_mod(end) < tol_mod);
%!    assert(err_w(end) < tol_w);
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!test
%!  ## TEST 10
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    unwind_protect
%!      [fd, msg] = fopen([filename, ".geo"], "w");
%!      if (fd == -1)
%!        error("failed to open file \"%s.geo\"", filename);
%!      endif
%!      d1 = 8e-3;
%!      D1 = 12e-3;
%!      w = 5e-3;
%!      l = 47e-3;
%!      h = 5e-3;
%!      grp_id_p1 = 2;
%!      grp_id_clamp = 3;
%!      grp_id_itf2 = 4;
%!      scale_def = 5e-3;
%!      mesh_size = 2e-3;
%!      fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!      fprintf(fd, "d1 = %g;\n", d1);
%!      fprintf(fd, "D1 = %g;\n", D1);
%!      fprintf(fd, "w = %g;\n", w);
%!      fprintf(fd, "l = %g;\n", l);
%!      fprintf(fd, "h = %g;\n", h);
%!      fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!      fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(6)  = {           0,  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(7)  = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(8) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!      fputs(fd, "Point(9) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(10) = {2 * h, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(11) = {2 * h, 4.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(12) = { h, 4.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(13) = { h, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Point(14) = {          0, 0.5 * h, -0.5 * w};\n");
%!      fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!      fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!      fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!      fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!      fputs(fd, "Circle(5) = {7, 1, 8};\n");
%!      fputs(fd, "Circle(6) = {8, 1, 9};\n");
%!      fputs(fd, "Line(7) = {9, 10};\n");
%!      fputs(fd, "Line(8) = {10, 11};\n");
%!      fputs(fd, "Line(9) = {11, 12};\n");
%!      fputs(fd, "Line(10) = {12, 13};\n");
%!      fputs(fd, "Line(11) = {13, 14};\n");
%!      fputs(fd, "Line(12) = {14, 6};\n");
%!      fputs(fd, "Line(13) = {6, 7};\n");
%!      fputs(fd, "Curve Loop(14) = {5, 6, 7, 8, 9, 10, 11, 12, 13};\n");
%!      fputs(fd, "Curve Loop(15) = {1, 2, 3, 4};\n");
%!      fputs(fd, "Plane Surface(16) = {14, 15};\n");
%!      fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{16}; };\n");
%!      fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!      fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!      fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!      fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[11],tmp[12],tmp[13],tmp[14]};\n", grp_id_p1);
%!      fprintf(fd, "Physical Surface(\"clamp\", %d) = {tmp[9]};\n", grp_id_clamp);
%!      fprintf(fd, "Physical Surface(\"itf2\", %d) = {tmp[6]};\n", grp_id_itf2);
%!    unwind_protect_cleanup
%!      if (fd ~= -1)
%!        fclose(fd);
%!        fd = -1;
%!      endif
%!    end_unwind_protect
%!    fprintf(stderr, "meshing ...\n");
%!    pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!    status = spawn_wait(pid);
%!    if (status ~= 0)
%!      error("gmsh failed with status %d", status);
%!    endif
%!    fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!    mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!    fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!    grp_idx_p1 = find([[mesh.groups.tria10].id] == grp_id_p1);
%!    cms_opt.nodes.modal.number = rows(mesh.nodes) + 3;
%!    cms_opt.nodes.interfaces(1).number = rows(mesh.nodes) + 1;
%!    cms_opt.nodes.interfaces(2).number = rows(mesh.nodes) + 2;
%!    bearing_surf(1).group_idx = grp_idx_p1;
%!    bearing_surf(1).options.reference_pressure = 1e9;
%!    bearing_surf(1).options.mesh_size = 1e-3;
%!    bearing_surf(1).options.include_rigid_body_modes = true;
%!    bearing_surf(1).options.bearing_type = "shell";
%!    bearing_surf(1).options.matrix_type = "modal substruct total";
%!    bearing_surf(1).r = 0.5 * d1;
%!    bearing_surf(1).w = w;
%!    bearing_surf(1).X0 = [l; 0; 0];
%!    bearing_surf(1).R = eye(3);
%!    bearing_surf(1).relative_tolerance = 0;
%!    bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d1;
%!    bearing_surf(1).master_node_no = cms_opt.nodes.interfaces(1).number;
%!    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!    mesh.nodes(cms_opt.nodes.interfaces(1).number, 1:3) = bearing_surf(1).X0.';
%!    mesh.nodes(cms_opt.nodes.interfaces(2).number, 1:3) = [1.5 * h, 4.5 * h, 0];
%!    mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_clamp, cms_opt.nodes.modal.number, "tria10");
%!    mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces(1).number, "tria10");
%!    mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_itf2, cms_opt.nodes.interfaces(2).number, "tria10");
%!    mesh.elements.joints(1).nodes = cms_opt.nodes.modal.number;
%!    mesh.elements.joints(1).C = eye(6);
%!    mesh.elements.joints(2).nodes = cms_opt.nodes.interfaces(2).number;
%!    mesh.elements.joints(2).C = eye(6);
%!    cms_opt.inveriants = true;
%!    cms_opt.static_modes = false;
%!    cms_opt.modal_node_constraint = false;
%!    cms_opt.load_cases = "index";
%!    cms_opt.solver = "umfpack";
%!    cms_opt.refine_max_iter = int32(10);
%!    cms_opt.verbose = int32(0);
%!    load_case = fem_pre_load_case_create_empty(6);
%!    load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!    zero_U = struct("U", cellfun(@(C) zeros(rows(C), 1), {mesh.elements.joints.C}, "UniformOutput", false));
%!    for i=1:6
%!      load_case(i).joints = zero_U;
%!      load_case(i).joints(2).U(i) = 1;
%!    endfor
%!    mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!    mesh.material_data.E = 210000e6;
%!    mesh.material_data.nu = 0.3;
%!    mesh.material_data.rho = 7850;
%!    opt_modes.shift_A = 0;
%!    opt_modes.refine_max_iter = cms_opt.refine_max_iter;
%!    opt_modes.verbose = cms_opt.verbose;
%!    opt_modes.solver = cms_opt.solver;
%!    opt_modes.active_joint_idx_eig = [1:2];
%!    opt_modes.rigid_body_modes_load_index = 1:6;
%!    opt_modes.elem_type = "tria10";
%!    num_modes = [10];
%!    interfaces = {"flexible"};
%!    num_modes_cms = int32([10]);
%!    num_cases = numel(num_modes) * numel(interfaces) * numel(num_modes_cms);
%!    err_red = zeros(7, num_cases);
%!    err_mod = zeros(1, num_cases);
%!    err_w = zeros(1, num_cases);
%!    k1 = 1;
%!    icase = int32(0);
%!    for k=1:numel(num_modes_cms)
%!      cms_opt.modes.number = num_modes_cms(k);
%!      for j=1:numel(interfaces)
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        for l=1:numel(num_modes)
%!          for i=1:numel(bearing_surf)
%!            bearing_surf(i).options.number_of_modes = min(num_modes(l), floor(numel(mesh.groups.tria10(bearing_surf(i).group_idx).nodes) * 3 / 2));
%!          endfor
%!          [mesh_l, load_case_bearing, bearing_surf_l, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!          [mesh_l, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh_l, load_case_bearing, cms_opt);
%!          assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!          [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!          [lambda, idx_lambda] = sort(diag(lambda));
%!          qred = qred(:, idx_lambda);
%!          sol_eig_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          for i=1:size(sol_eig_red.def, 3)
%!            sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!          endfor
%!          sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!          comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh_l, mat_ass, dof_map, cms_opt, bearing_surf_l);
%!          load_case_itf = fem_pre_load_case_create_empty(6);
%!          for i=1:6
%!            load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces(1).number;
%!            load_case_itf(i).loads = zeros(1, 6);
%!            load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!            switch (i)
%!              case {4, 5}
%!                load_case_itf(i).loads(i) *= w;
%!              case 6
%!                load_case_itf(i).loads(i) *= d1;
%!            endswitch
%!          endfor
%!          [Kitf, Ritf] = fem_ass_matrix(mesh_l, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!          nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!          nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!          X1 = mesh_l.nodes(mesh_l.groups.tria10(bearing_surf_l(1).group_idx).nodes, 1:3) - bearing_surf_l(1).X0.';
%!          Phi1 = atan2(X1(:, 2), X1(:, 1));
%!          p1 = k1 * sin(pi/4 + Phi1) * bearing_surf_l(1).options.reference_pressure;
%!          Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!          z1g = comp_mat(1).bearing_surf.grid_z(:);
%!          p1red = zeros(numel(z1g), numel(Phi1g));
%!          for i=1:columns(p1red)
%!            p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!          endfor
%!          p1red = p1red(:);
%!          Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf_l(1).options.reference_pressure;
%!          Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!          qred = mat_ass.Kred \ Fred;
%!          w1red = comp_mat(1).D * qred;
%!          sol_red.def = fem_post_cms_expand_body(mesh_l, dof_map, mat_ass, qred);
%!          mesh_post = mesh_l;
%!          mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!          load_case_post = fem_pre_load_case_create_empty(7);
%!          load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!          for i=1:6
%!            load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!            load_case_post(i).loads = load_case_itf(i).loads;
%!          endfor
%!          load_case_post(7).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!          pn = zeros(rows(mesh_l.nodes), 1);
%!          pn(mesh_post.groups.tria10(grp_idx_p1).nodes) = p1;
%!          load_case_post(7).pressure.tria10.p = pn(mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :));
%!          mesh_post.elements.joints.C = eye(6);
%!          mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!          dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!          [mat_ass_post.M, ...
%!           mat_ass_post.K, ...
%!           mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                            dof_map_post, ...
%!                                            [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_post);
%!          sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!          node_idx1 = mesh_l.groups.tria10(bearing_surf(1).group_idx).nodes;
%!          w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!          for i=1:numel(node_idx1)
%!            ni = [mesh_l.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!            ni /= norm(ni);
%!            w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!          endfor
%!          w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!          for i=1:columns(w1postint)
%!            for m=1:numel(Phi1g)
%!              w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                             repmat(X1(:, 3), 3, 1), ...
%!                                                                             repmat(w1post(:, i), 3, 1), ...
%!                                                                             repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                             z1g);
%!            endfor
%!          endfor
%!          if (do_plot)
%!            for m=1:columns(w1red)
%!              for i=1:numel(z1g)
%!                figure("visible", "off");
%!                hold("on");
%!                plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;1");
%!                plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;0");
%!                ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!                title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!                xlabel("Phi [deg]");
%!                ylabel("w [um]");
%!              endfor
%!            endfor
%!          endif
%!          sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!          for i=1:size(sol_eig_post.def, 3)
%!            sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!          endfor
%!          num_modes_comp = min(40, floor(numel(sol_eig_red.f) / 3));
%!          err_mod(++icase) = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!          err_w(icase) = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!          mesh_data(1).mesh = mesh_l;
%!          mesh_data(1).dof_map = dof_map;
%!          mesh_data(2).mesh = mesh_post;
%!          mesh_data(2).dof_map = dof_map_post;
%!          mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!          [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!          sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!          sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!          for i=1:size(sol_post.def, 3)
%!            err_red(i, icase) = norm(sol_post.def(:, :, i) - sol_red.def(:, :, i)) / norm(sol_post.def(:, :, i));
%!          endfor
%!        endfor
%!        for l=1:numel(num_modes)
%!          fprintf(stderr, "%s interfaces using %d static pressure modes, %d dynamic modes:\n", interfaces{j}, num_modes(l), num_modes_cms(k));
%!          for i=1:size(sol_post.def, 3)
%!            fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i, icase));
%!          endfor
%!          fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod(icase));
%!          fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w(icase));
%!        endfor
%!      endfor
%!    endfor
%!    tol_red = 1e-2;
%!    tol_mod = 1e-2;
%!    tol_w = 1e-2;
%!    assert(all(err_red(:, end) < tol_red));
%!    assert(err_mod(end) < tol_mod);
%!    assert(err_w(end) < tol_w);
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect

%!demo
%!  ## DEMO1
%!  do_plot = false;
%!  if (do_plot)
%!    close all;
%!  endif
%!  fd = -1;
%!  filename = "";
%!  unwind_protect
%!    filename = tempname();
%!    if (ispc())
%!      filename(filename == "\\") = "/";
%!    endif
%!    interfaces = {"flexible"};
%!    tol_red = [2.5e-2];
%!    num_modes_cms = int32([10]);
%!    num_modes = [15];
%!    for k=1:numel(num_modes_cms)
%!      for j=1:numel(interfaces)
%!        clear Fred Ritf bearing_surf cms_opt comp_mat dof_map_comb dof_map_post err_red
%!        clear grp_id_p1 grp_id_p2 grp_idx_p1 grp_idx_p2 load_case load_case_bearing load_case_itf
%!        clear load_case_post mat_ass_post mat_ass_press mesh mesh_comb mesh_data mesh_post mesh_size
%!        clear opt_modes p1 p1red p2 p2red pid qred sol_comb sol_eig sol_eig_cms sol_post sol_red
%!        unwind_protect
%!          [fd, msg] = fopen([filename, ".geo"], "w");
%!          if (fd == -1)
%!            error("failed to open file \"%s.geo\"", filename);
%!          endif
%!          d = 14e-3;
%!          D = 19.5e-3;
%!          w = 5e-3;
%!          l = 47e-3;
%!          h = 5e-3;
%!          grp_id_p1 = 2;
%!          grp_id_p2 = 3;
%!          p1 = 1;
%!          p2 = 2;
%!          scale_def = 5e-3;
%!          mesh_size = 1.5e-3;
%!          fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!          fprintf(fd, "d = %g;\n", d);
%!          fprintf(fd, "D = %g;\n", D);
%!          fprintf(fd, "w = %g;\n", w);
%!          fprintf(fd, "l = %g;\n", l);
%!          fprintf(fd, "h = %g;\n", h);
%!          fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!          fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!          fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!          fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!          fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!          fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!          fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!          fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!          fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!          fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!          fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!          fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!          fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!          fputs(fd, "Line(11) = {13, 14};\n");
%!          fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!          fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!          fputs(fd, "Line(14) = {16, 11};\n");
%!          fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!          fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!          fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!          fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!          fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!          fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!          fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!          fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!          fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!          fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!        unwind_protect_cleanup
%!          if (fd ~= -1)
%!            fclose(fd);
%!            fd = -1;
%!          endif
%!        end_unwind_protect
%!        fprintf(stderr, "meshing ...\n");
%!        pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!        status = spawn_wait(pid);
%!        if (status ~= 0)
%!          error("gmsh failed with status %d", status);
%!        endif
%!        fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!        mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!        fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!        grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!        grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!        cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!        cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!        cms_opt.floating_frame = false;
%!        cms_opt.algorithm = "diag-shift-invert";
%!        bearing_surf(1).group_idx = grp_idx_p1;
%!        bearing_surf(1).options.reference_pressure = 1e9;
%!        bearing_surf(1).options.mesh_size = 1e-3;
%!        bearing_surf(1).options.include_rigid_body_modes = true;
%!        bearing_surf(1).options.bearing_type = "shell";
%!        bearing_surf(1).options.matrix_type = "modal substruct total";
%!        bearing_surf(1).r = 0.5 * d;
%!        bearing_surf(1).w = w;
%!        bearing_surf(1).X0 = [l; 0; 0];
%!        bearing_surf(1).R = eye(3);
%!        bearing_surf(1).relative_tolerance = 0;
%!        bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!        bearing_surf(1).options.number_of_modes = num_modes(j);
%!        bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!        bearing_surf(2).group_idx = grp_idx_p2;
%!        bearing_surf(2).options.reference_pressure = 1e9;
%!        bearing_surf(2).options.mesh_size = 1e-3;
%!        bearing_surf(2).options.include_rigid_body_modes = false;
%!        bearing_surf(2).options.bearing_type = "shell";
%!        bearing_surf(2).options.matrix_type = "modal substruct total";
%!        bearing_surf(2).r = 0.5 * d;
%!        bearing_surf(2).w = w;
%!        bearing_surf(2).X0 = [0; 0; 0];
%!        bearing_surf(2).R = eye(3);
%!        bearing_surf(2).relative_tolerance = 0;
%!        bearing_surf(2).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!        bearing_surf(2).options.number_of_modes = num_modes;
%!        bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!        mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!        mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!        mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!        mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!        cms_opt.inveriants = true;
%!        cms_opt.modes.number = num_modes_cms(k);
%!        cms_opt.static_modes = false;
%!        cms_opt.modal_node_constraint = false;
%!        cms_opt.load_cases = "index";
%!        cms_opt.refine_max_iter = int32(10);
%!        load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!        mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!        mesh.material_data.E = 210000e6;
%!        mesh.material_data.nu = 0.3;
%!        mesh.material_data.rho = 7850;
%!        opt_modes.shift_A = 1e-6;
%!        opt_modes.refine_max_iter = int32(10);
%!        opt_modes.verbose = int32(0);
%!        opt_modes.rigid_body_modes = interfaces{j};
%!        opt_modes.solver = "pastix";
%!        opt_modes.number_of_threads = int32(4);
%!        cms_opt.solver = opt_modes.solver;
%!        cms_opt.number_of_threads = opt_modes.number_of_threads;
%!        [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!        [mesh, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!        assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!        assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!        [qred, lambda_red] = eig(mat_ass.Kred, mat_ass.Mred);
%!        [lambda_red, idx_lambda_red] = sort(diag(lambda_red));
%!        qred = qred(:, idx_lambda_red);
%!        sol_red_modal.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!        sol_red_modal.f = sqrt(lambda_red) / (2 * pi);
%!        for i=1:size(sol_red_modal.def, 3)
%!          sol_red_modal.def(:, :, i) *= 10e-3 / max(max(abs(sol_red_modal.def(:, 1:3, i))));
%!        endfor
%!        comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%!        load_case_itf = fem_pre_load_case_create_empty(6);
%!        for i=1:6
%!          load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!          load_case_itf(i).loads = zeros(1, 6);
%!          load_case_itf(i).loads(i) = 1;
%!        endfor
%!        [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!        nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!        nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!        nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!        nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!        p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!        p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!        p1red2 = zeros((nx1 - 1) * nz1, 1);
%!        p2red2 = zeros((nx2 - 1) * nz2, 1);
%!        for i=1:nx1 - 1
%!          p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!        endfor
%!        for i=1:nx2 - 1
%!          p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!        endfor
%!        Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!                comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!        Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!        qred = mat_ass.Kred \ Fred;
%!        sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!        mesh_post = mesh;
%!        mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!        load_case_post = fem_pre_load_case_create_empty(10);
%!        load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!        for i=1:6
%!          load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!          load_case_post(i).loads = zeros(1, 6);
%!          load_case_post(i).loads(i) = 1;
%!        endfor
%!        x1 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!        y1 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!        x2 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!        y2 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!        Phi1 = atan2(y1, x1);
%!        Phi2 = atan2(y2, x2);
%!        load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!        load_case_post(7).pressure.tria6.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%!        load_case_post(8).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!        load_case_post(8).pressure.tria6.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!        load_case_post(9).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!        load_case_post(9).pressure.tria6.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6);
%!        load_case_post(10).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!        load_case_post(10).pressure.tria6.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!        mesh_post.elements.joints.C = eye(6);
%!        mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!        dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!        dof_map_post.parallel.threads_ass = cms_opt.number_of_threads;
%!        [mat_ass_post.M, ...
%!         mat_ass_post.K, ...
%!         mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                          dof_map_post, ...
%!                                          [FEM_MAT_MASS, ...
%!                                           FEM_MAT_STIFFNESS, ...
%!                                           FEM_VEC_LOAD_CONSISTENT], ...
%!                                          load_case_post);
%!        sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, opt_modes);
%!        sol_post_modal = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, cms_opt.modes.number + 6, 0, sqrt(eps), "shift-invert", opt_modes);
%!        for i=1:size(sol_post_modal.def, 3)
%!          sol_post_modal.def(:, :, i) *= 10e-3 / max(max(abs(sol_post_modal.def(:, 1:3, i))));
%!          if (norm(sol_post_modal.def(:, :, i) + sol_red_modal.def(:, :, i)) < norm(sol_post_modal.def(:, :, i) - sol_red_modal.def(:, :, i)))
%!            sol_post_modal.def(:, :, i) *= -1;
%!          endif
%!        endfor
%!        mesh_data(1).mesh = mesh;
%!        mesh_data(1).dof_map = dof_map;
%!        mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!        mesh_data(2).mesh = mesh_post;
%!        mesh_data(2).dof_map = dof_map_post;
%!        [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!        sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!        sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!        sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!        for i=1:size(sol_comb.def, 3)
%!          sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!        endfor
%!        sol_comb_modal.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post_modal.def, 3));
%!        sol_comb_modal.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post_modal.def)), :, :) = sol_post_modal.def;
%!        sol_comb_modal.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red_modal.def)), :, :) = sol_red_modal.def(:, :, 1:size(sol_post_modal.def, 3));
%!        for i=1:size(sol_comb.def, 3)
%!          sol_comb_modal.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!        endfor
%!        err_red = zeros(1, size(sol_post.def, 3));
%!        for i=1:size(sol_post.def, 3)
%!          err_red(i) = max(max(abs(sol_post.def(1:end - 2, :, i) - sol_red.def(1:end - 2, :, i)))) / max(max(abs(sol_post.def(1:end - 2, :, i))));
%!        endfor
%!        for i=1:size(sol_post.def, 3)
%!          fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!        endfor
%!        err_mod_freq = (sol_red_modal.f(1:numel(sol_post_modal.f)) ./ sol_post_modal.f - 1);
%!        MAC = zeros(1, numel(sol_post_modal.f));
%!        for i=1:numel(sol_post_modal.f)
%!          ve = sol_post_modal.def(:, :, i)(:);
%!          vo = sol_red_modal.def(:, :, i)(:);
%!          MAC(i) = (ve.' * vo)^2 / ((ve.' * ve) * (vo.' * vo));
%!        endfor
%!        for i=1:numel(sol_post_modal.f)
%!          fprintf(stderr, ...
%!                  "mode %d: reduced: %.1fHz full: %.1fHz difference freq %.1f%% MAC %.4f\n", ...
%!                  i, ...
%!                  sol_red_modal.f(i), ...
%!                  sol_post_modal.f(i), ...
%!                  100 * (sol_red_modal.f(i) / sol_post_modal.f(i) - 1), ...
%!                  MAC(i));
%!        endfor
%!        assert(all(err_red < tol_red(j)));
%!      endfor
%!    endfor
%!  unwind_protect_cleanup
%!    if (numel(filename))
%!      fn = dir([filename, "*"]);
%!      for i=1:numel(fn)
%!        status = unlink(fullfile(fn(i).folder, fn(i).name));
%!        if (status ~= 0)
%!          warning("failed to remove file \"%s\"", fn(i).name);
%!        endif
%!      endfor
%!    endif
%!  end_unwind_protect
