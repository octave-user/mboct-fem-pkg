function [mesh, load_case, bearing_surf, idx_modes, sol_eig, mat_ass_press, dof_map_press] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, options)
  if (nargin < 4)
    options = struct();
  endif

  if (~isfield(options, "shift_A"))
    options.shift_A = 0;
  endif

  if (~isfield(options, "solver"))
    options.solver = struct();
  endif

  [load_case_press_dist, bearing_surf, idx_group] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);

  load_case_press_dist(1).locked_dof = load_case.locked_dof;
  
  inum_elem_press = int32(0);
  
  for i=1:numel(bearing_surf)
    inum_elem_press += numel(mesh.groups.tria6(bearing_surf(i).group_idx).elements);
    bearing_surf(i).nodes = mesh.groups.tria6(bearing_surf(i).group_idx).nodes;
  endfor
  
  elements = zeros(inum_elem_press, 6, "int32");

  inum_elem_press = int32(0);
  ielem_idx = zeros(numel(bearing_surf), 2, "int32");
  
  for i=1:numel(bearing_surf)
    elem_idx = mesh.groups.tria6(bearing_surf(i).group_idx).elements;
    idx_slot = inum_elem_press + (1:numel(elem_idx));
    elements(idx_slot, :) = mesh.elements.tria6(elem_idx, :);
    ielem_idx(i, :) = idx_slot([1, end]);
    inum_elem_press += numel(elem_idx);
  endfor
  
  load_case_press_mod = fem_pre_load_case_create_empty(numel(bearing_surf));
  
  for i=1:numel(bearing_surf)
    load_case_press_mod(i).pressure.tria6.elements = elements;
    load_case_press_mod(i).pressure.tria6.p = zeros(size(elements));
    load_case_press_mod(i).pressure.tria6.p(ielem_idx(i, 1):ielem_idx(i, end), :) = 1;
  endfor
  
  dof_map_press = fem_ass_dof_map(mesh, load_case);

  [mat_ass_press.K, ...
   mat_ass_press.R] = fem_ass_matrix(mesh, ...
				     dof_map_press, ...
				     [FEM_MAT_STIFFNESS, ...
				      FEM_VEC_LOAD_CONSISTENT], ...
				     load_case_press_mod);

  
  diagA = zeros(rows(mesh.nodes), numel(bearing_surf));

  for i=1:columns(dof_map_press.ndof)
    dof_idx = dof_map_press.ndof(:, i);
    idx_act_dof = find(dof_idx > 0);
    diagA(idx_act_dof, :) += mat_ass_press.R(dof_idx(idx_act_dof), :).^2;
  endfor
  
  inum_modes_tot = int32(0);
  inum_modes_press = int32(0);
  
  for i=1:numel(bearing_surf)
    inum_modes_tot += bearing_surf(i).number_of_modes;
    inum_modes_press += bearing_surf(i).number_of_modes;

    if (~isfield(bearing_surf(i).options, "include_rigid_body_modes"))
      bearing_surf(i).options.include_rigid_body_modes = true;
    endif

    if (bearing_surf(i).options.include_rigid_body_modes)
      inum_modes_tot += 6;
    endif
    
    if (options.shift_A)
      inum_modes_tot -= 6;
    endif
  endfor

  if (nargout >= 5)
    sol_eig.def = zeros(rows(mesh.nodes), columns(mesh.nodes), inum_modes_press);
    sol_eig.lambda = zeros(1, inum_modes_tot);
  endif
  
  if (options.shift_A == 0)
    Kfact = fem_sol_factor(mat_ass_press.K, options.solver);
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

  if (isfield(load_case, "joints"))
    ijoint_offset = numel(load_case.joints);
  else
    ijoint_offset = 0;
  endif
  
  inum_joints = ijoint_offset;
  ijoint_idx = zeros(numel(bearing_surf), 2, "int32");
  
  for i=1:numel(bearing_surf)
    idx_joint = 1:numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes);
    ijoint_idx(i, :) = inum_joints + idx_joint([1, end]);
    inum_joints += numel(idx_joint);
  endfor

  empty_cell = cell(1, inum_joints);

  elem_joints = struct("nodes", empty_cell, "C", empty_cell);

  for i=1:ijoint_offset
    elem_joints(i) = mesh.joints(i);
  endfor
  
  for i=1:numel(bearing_surf)
    node_idx = mesh.groups.tria6(bearing_surf(i).group_idx).nodes;
    for j=1:numel(node_idx)
      elem_joints(ijoint_idx(i, 1) + j - 1).nodes = node_idx(j);
      elem_joints(ijoint_idx(i, 1) + j - 1).C = [eye(3), zeros(3, 3)];
    endfor
  endfor

  load_case_displ = fem_pre_load_case_create_empty(inum_modes_tot);

  for i=1:numel(load_case_displ)
    load_case_displ(i).joints = struct("U", repmat({zeros(3, 1)}, 1, inum_joints));
  endfor
  
  inum_modes_tot = int32(0);
  inum_modes_press = int32(0);
  
  for i=1:numel(bearing_surf)
    Ap = zeros(dof_map_press.totdof, 1);
    
    for j=1:columns(dof_map_press.ndof)
      dof_idx = dof_map_press.ndof(:, j);
      idx_act_dof = find(dof_idx > 0);      
      Ap(dof_idx(idx_act_dof)) = sqrt(diagA(idx_act_dof, i));
    endfor

    Ap = diag(Ap);

    if (options.shift_A == 0)
      Kp = mat_ass_press.K;
      gamma = 0;
    else
      gamma = options.shift_A * max(max(abs(mat_ass_press.K))) / max(abs(diag(Ap)));
      Kp = mat_ass_press.K + gamma * Ap;
      Kfact = fem_sol_factor(Kp, options.solver);
    endif
      
    oper = cell(1, 2);  
    oper{1} = @(x) Ap * x;
    oper{2} = @(x) Kfact \ x;
    
    [Phi, kappa] = eig_sym(oper, columns(Kp), bearing_surf(i).number_of_modes, sigma, opt);

    inode_idx_bs = mesh.groups.tria6(bearing_surf(i).group_idx).nodes;

    U = zeros(numel(inode_idx_bs) * 3, columns(Phi));

    for j=1:3
      dof_idx = dof_map_press.ndof(inode_idx_bs, j);
      idx_act_dof = find(dof_idx > 0);
      U((idx_act_dof - 1) * 3 + j, :) = Phi(dof_idx(idx_act_dof), :);
    endfor

    Arb = [repmat(eye(3), numel(inode_idx_bs), 1), zeros(3 * numel(inode_idx_bs), 3)];

    l = mesh.nodes(inode_idx_bs, 1:3) - bearing_surf(i).X0.';

    for j=1:3
      for k=1:2
	Arb(ir(j, k):3:end, ic(j, k) + 3) = -dc(k) * l(:, j);
      endfor
    endfor

    q = (Arb.' * Arb) \ (Arb.' * U);
    U -= Arb * q;

    if (options.shift_A)
      U = U(:, 7:end);
    endif
    
    if (bearing_surf(i).options.include_rigid_body_modes)
      U = [Arb, U];
    endif

    U *= diag(1 ./ max(abs(U), [], 1));
    
    for k=1:columns(U)
      for j=1:numel(inode_idx_bs)
	load_case_displ(inum_modes_tot + k).joints(ijoint_idx(i, 1) + j - 1).U = U((j - 1) * 3 + (1:3), k);
      endfor
    endfor
    
    if (nargout >= 5)
      idx_mode_press = inum_modes_press + (1:columns(Phi));
      
      for j=1:columns(dof_map_press.ndof)
	idof_idx = dof_map_press.ndof(:, j);
	iactive_dof = find(idof_idx > 0);
	idof_idx = idof_idx(iactive_dof);
	
	sol_eig.def(iactive_dof, j, idx_mode_press) = Phi(idof_idx, :);
      endfor

      sol_eig.lambda(idx_mode_press) = diag(kappa) - gamma;
      inum_modes_press += columns(Phi);
    endif
    
    inum_modes_tot += columns(U);
  endfor

  mesh.elements.joints = elem_joints;
  
  load_case = fem_pre_load_case_merge(load_case_press_dist, load_case_displ);

  idx_modes = numel(load_case_press_dist) + (1:numel(load_case_displ));

  if (nargout >= 5)
    for i=1:size(sol_eig.def, 3)
      sol_eig.def(:, :, i) /= max(max(abs(sol_eig.def(:, :, i))));
    endfor
  endif
endfunction

%!demo
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! scale_def = 5e-3;
%! mesh_size = 1.25e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!  if (fd ~= -1)
%!    fclose(fd);
%!  endif
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");

%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  error("gmsh failed with status %d", status);
%! endif

%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! bearing_surf(1).group_idx = grp_id_p1;
%! bearing_surf(1).options.reference_pressure = 1e9;
%! bearing_surf(1).options.mesh_size = 20e-3;
%! bearing_surf(1).options.include_rigid_body_modes = false;
%! bearing_surf(1).options.bearing_type = "shell";
%! bearing_surf(1).options.matrix_type = "modal substruct total";
%! bearing_surf(1).r = ri;
%! bearing_surf(1).w = b;
%! bearing_surf(1).X0 = [0; 0; b/2 + c];
%! bearing_surf(1).R = eye(3);
%! bearing_surf(1).relative_tolerance = 0;
%! bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! bearing_surf(1).number_of_modes = 10;
%! bearing_surf(2).group_idx = grp_id_p2;
%! bearing_surf(2).options.reference_pressure = 1e9;
%! bearing_surf(2).options.mesh_size = 20e-3;
%! bearing_surf(2).options.include_rigid_body_modes = true;
%! bearing_surf(2).options.bearing_type = "journal";
%! bearing_surf(2).options.matrix_type = "modal substruct total";
%! bearing_surf(2).r = ro;
%! bearing_surf(2).w = b;
%! bearing_surf(2).X0 = [0; 0; b/2 + c];
%! bearing_surf(2).R = eye(3);
%! bearing_surf(2).relative_tolerance = 0;
%! bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! bearing_surf(2).number_of_modes = 12;
%! cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%! mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(1).X0.';
%! cms_opt.inveriants = true;
%! cms_opt.modes.number = 10;
%! cms_opt.static_modes = false;
%! cms_opt.load_cases = "index";
%! cms_opt.refine_max_iter = int32(10);
%! f_enable_constraint = false;
%! load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%! load_case(1).locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%! if (f_enable_constraint)
%!   load_case(1).locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! endif
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! if (~f_enable_constraint)
%!   opt_modes.shift_A = 1e-10;
%! endif
%! opt_modes.solver.refine_max_iter = int32(10);
%! opt_modes.solver.verbose = PASTIX_API_VERBOSE_NOT;
%! #debug_on_error(true);
%! #dbstop("fem_ehd_pre_comp_mat_unstruct", 52);
%! #dbstop("fem_cms_create",45);
%! [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.R, ...
%!  mat_info, ...
%!  mesh_info] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_MAT_MASS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case_bearing);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! [mesh, ...
%!  mat_ass, ...
%!  dof_map, ...
%!  sol_eig_cms, ...
%!  cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%! comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
