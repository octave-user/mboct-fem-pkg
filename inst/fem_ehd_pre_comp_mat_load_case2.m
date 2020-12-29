function [mesh, load_case, bearing_surf, idx_modes, sol_eig, mat_ass_press, dof_map_press] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, options)
  if (nargin < 4)
    options = struct();
  endif

  if (~isfield(options, "shift_A"))
    options.shift_A = 0;
  endif

  if (~isfield(options, "solver"))
    options.solver = "pastix";
  endif

  if (~isfield(options, "refine_max_iter"))
    options.refine_max_iter = int32(10);
  endif

  if (~isfield(options, "verbose"))
    options.verbose = int32(0);
  endif

  if (~isfield(options, "rigid_body_modes"))
    options.rigid_body_modes = "rigid";
  endif

  switch (options.rigid_body_modes)
    case "flexible"
      if (~isfield(bearing_surf, "master_node_no"))
	error("missing field bearing_surf.master_node_no needed for flexible rigid body modes");
      endif

      rigid_body_modes_flexible = true;
    case "rigid"
      rigid_body_modes_flexible = false;
    otherwise
      error("invalid option rigid_body_modes=\"%s\"", options.rigid_body_modes);
  endswitch

  options.solver = fem_sol_select(true, options.solver);

  switch (options.solver)
    case {"umfpack", "lu", "mldivide"}
      mat_type_stiffness = FEM_MAT_STIFFNESS;
    otherwise
      mat_type_stiffness = FEM_MAT_STIFFNESS_SYM_L;
  endswitch

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
				     [mat_type_stiffness, ...
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

  if (options.shift_A == 0)
    Kfact = fem_sol_factor(mat_ass_press.K, options);
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

  if (isfield(mesh.elements, "joints"))
    ijoint_offset = numel(mesh.elements.joints);
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
    elem_joints(i) = mesh.elements.joints(i);
  endfor

  for i=1:numel(bearing_surf)
    node_idx = mesh.groups.tria6(bearing_surf(i).group_idx).nodes;
    for j=1:numel(node_idx)
      elem_joints(ijoint_idx(i, 1) + j - 1).nodes = node_idx(j);
      elem_joints(ijoint_idx(i, 1) + j - 1).C = [eye(3), zeros(3, 3)];
    endfor
  endfor

  load_case_displ = fem_pre_load_case_create_empty(inum_modes_tot);

  for j=1:numel(elem_joints)
    U = zeros(rows(elem_joints(j).C), 1);

    for i=1:numel(load_case_displ)
      load_case_displ(i).joints(j).U = U;
    endfor
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
      gamma = 0;
    else
      gamma = options.shift_A * max(max(abs(mat_ass_press.K))) / max(abs(diag(Ap)));
      Kfact = fem_sol_factor(mat_ass_press.K + gamma * Ap, options);
    endif

    oper = cell(1, 2);
    oper{1} = @(x) Ap * x;
    oper{2} = @(x) Kfact \ x;

    num_modes = bearing_surf(i).options.number_of_modes;

    if (options.shift_A ~= 0)
      num_modes += 6;
    endif

    [Phi, kappa] = eig_sym(oper, columns(mat_ass_press.K), num_modes, sigma, opt);

    clear oper Ap;

    if (options.shift_A ~= 0)
      clear Kfact;
    endif

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

    if (bearing_surf(i).options.include_rigid_body_modes && ~rigid_body_modes_flexible)
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
	clear idof_idx iactive_dof;
      endfor

      sol_eig.lambda(idx_mode_press) = diag(kappa) - gamma;
      inum_modes_press += columns(Phi);
    endif

    inum_modes_tot += columns(U);
    clear U;
  endfor

  if (nargout < 6)
    clear mat_ass_press;
  endif

  if (rigid_body_modes_flexible)
    inum_modes_itf_flex = int32(0);

    for i=1:numel(bearing_surf)
      mesh.elements.joints(ijoint_offset + i).nodes = bearing_surf(i).master_node_no;
      mesh.elements.joints(ijoint_offset + i).C = eye(6);

      if (bearing_surf(i).options.include_rigid_body_modes)
	inum_modes_itf_flex += 6;
      endif
    endfor

    load_case_itf_flex = fem_pre_load_case_create_empty(inum_modes_itf_flex);

    zero_U = struct("U", repmat({zeros(6, 1)}, 1, numel(mesh.elements.joints)));

    idx_load_case_itf_flex = int32(0);

    for i=1:numel(bearing_surf)
      if (~bearing_surf(i).options.include_rigid_body_modes)
	continue;
      endif

      for j=1:6
	load_case_itf_flex(++idx_load_case_itf_flex).joints = zero_U;
	load_case_itf_flex(idx_load_case_itf_flex).joints(ijoint_offset + i).U(j) = 1;
      endfor
    endfor

    clear zero_U;

    dof_map_itf_flex = fem_ass_dof_map(mesh, load_case);

    [mat_ass_itf_flex.K, ...
     mat_ass_itf_flex.R] = fem_ass_matrix(mesh, ...
					  dof_map_itf_flex, ...
					  [mat_type_stiffness, ...
					   FEM_VEC_LOAD_CONSISTENT], ...
					  load_case_itf_flex);

    Kfact = fem_sol_factor(mat_ass_itf_flex.K, options);

    U = Kfact \ mat_ass_itf_flex.R;

    clear Kfact mat_ass_itf_flex load_case_itf_flex;

    for i=1:numel(bearing_surf)
      inode_idx_bs = mesh.groups.tria6(bearing_surf(i).group_idx).nodes;
      Ubs = zeros(3, numel(inode_idx_bs), columns(U));

      for j=1:3
	idof_idx = dof_map_itf_flex.ndof(inode_idx_bs, j);
	iactive_dof = find(idof_idx > 0);
	Ubs(j, iactive_dof, :) = U(idof_idx(iactive_dof), :);
      endfor

      for k=1:size(Ubs, 3)
	for j=1:numel(inode_idx_bs)
	  load_case_displ(inum_modes_tot + k).joints(ijoint_idx(i, 1) + j - 1).U = Ubs(:, j, k);
	endfor
      endfor
    endfor
  endif

  mesh.elements.joints = elem_joints;

  for j=1:numel(elem_joints)
    U = zeros(rows(elem_joints(j).C), 1);

    for i=1:numel(load_case_press_dist)
      load_case_press_dist(i).joints(j).U = U;
    endfor
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
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   f_enable_constraint = [false, true];
%!   interfaces = {"rigid", "flexible"};
%!   num_modes = int32([0, 10]);
%!   for l=1:numel(num_modes)
%!     for k=1:numel(f_enable_constraint)
%!       for j=1:numel(interfaces)
%! 	clear bearing_surf cms_opt comp_mat dof_map dof_map_press grp_id_clamp grp_id_p1 grp_id_p2;
%! 	clear load_case load_case_bearing mat_ass mat_ass_press mat_info mesh mesh_info mesh_size
%! 	clear opt_modes sol_eig sol_eig_cms sol_stat;
%! 	unwind_protect
%! 	  [fd, msg] = fopen([filename, ".geo"], "wt");
%! 	  if (fd == -1)
%! 	    error("failed to open file \"%s.geo\"", filename);
%! 	  endif
%! 	  ri = 8e-3;
%! 	  ro = 10e-3;
%! 	  h = 12e-3;
%! 	  c = 2e-3;
%! 	  b = h - 2 * c;
%! 	  scale_def = 5e-3;
%! 	  mesh_size = 3e-3;
%! 	  fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! 	  fprintf(fd, "ri = %g;\n", ri);
%! 	  fprintf(fd, "ro = %g;\n", ro);
%! 	  fprintf(fd, "h = %g;\n", h);
%! 	  fprintf(fd, "c = %g;\n", c);
%! 	  fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! 	  fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! 	  fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! 	  fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! 	  fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! 	  fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! 	  fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! 	  fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! 	  fputs(fd, "Line(1) = {1,2};\n");
%! 	  fputs(fd, "Line(2) = {2,3};\n");
%! 	  fputs(fd, "Line(3) = {3,4};\n");
%! 	  fputs(fd, "Line(4) = {4,5};\n");
%! 	  fputs(fd, "Line(5) = {5,6};\n");
%! 	  fputs(fd, "Line(6) = {6,7};\n");
%! 	  fputs(fd, "Line(7) = {7,8};\n");
%! 	  fputs(fd, "Line(8) = {8,1};\n");
%! 	  fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! 	  fputs(fd, "Plane Surface(6) = {5};\n");
%! 	  fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%! 	  fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%! 	  fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! 	  fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! 	  fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! 	  fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! 	unwind_protect_cleanup
%! 	  if (fd ~= -1)
%! 	    fclose(fd);
%! 	    fd = -1;
%! 	  endif
%! 	end_unwind_protect
%! 	fprintf(stderr, "meshing ...\n");
%! 	pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! 	status = spawn_wait(pid);
%! 	if (status ~= 0)
%! 	  error("gmsh failed with status %d", status);
%! 	endif
%! 	fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! 	mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! 	fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! 	cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%! 	switch (interfaces{j})
%! 	  case "flexible"
%! 	    cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%! 	endswitch
%! 	grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! 	grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! 	grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! 	bearing_surf(1).group_idx = grp_id_p1;
%! 	bearing_surf(1).options.reference_pressure = 1e9;
%! 	bearing_surf(1).options.mesh_size = 20e-3;
%! 	bearing_surf(1).options.include_rigid_body_modes = false;
%! 	bearing_surf(1).options.bearing_type = "shell";
%! 	bearing_surf(1).options.matrix_type = "modal substruct total";
%! 	bearing_surf(1).r = ri;
%! 	bearing_surf(1).w = b;
%! 	bearing_surf(1).X0 = [0; 0; b/2 + c];
%! 	bearing_surf(1).R = eye(3);
%! 	bearing_surf(1).relative_tolerance = 0;
%! 	bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! 	bearing_surf(1).options.number_of_modes = 10;
%! 	bearing_surf(2).group_idx = grp_id_p2;
%! 	bearing_surf(2).options.reference_pressure = 1e9;
%! 	bearing_surf(2).options.mesh_size = 20e-3;
%! 	bearing_surf(2).options.include_rigid_body_modes = true;
%! 	bearing_surf(2).options.bearing_type = "journal";
%! 	bearing_surf(2).options.matrix_type = "modal substruct total";
%! 	bearing_surf(2).r = ro;
%! 	bearing_surf(2).w = b;
%! 	bearing_surf(2).X0 = [0; 0; b/2 + c];
%! 	bearing_surf(2).R = eye(3);
%! 	bearing_surf(2).relative_tolerance = 0;
%! 	bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! 	bearing_surf(2).options.number_of_modes = 12;
%! 	switch (interfaces{j})
%! 	  case "flexible"
%! 	    bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%! 	    bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%! 	    for i=1:numel(bearing_surf)
%!               mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%! 	    endfor
%! 	    for i=1:numel(bearing_surf)
%!               mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no);
%! 	    endfor
%! 	  otherwise
%! 	    mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%! 	endswitch
%! 	cms_opt.inveriants = true;
%! 	cms_opt.modes.number = num_modes(l);
%! 	cms_opt.static_modes = false;
%! 	cms_opt.load_cases = "index";
%! 	cms_opt.refine_max_iter = int32(10);
%! 	load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%! 	switch (interfaces{j})
%! 	  case "flexible"
%! 	  otherwise
%! 	    load_case(1).locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%! 	endswitch
%! 	if (f_enable_constraint(k))
%! 	  load_case(1).locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! 	endif
%! 	mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! 	mesh.material_data.E = 210000e6;
%! 	mesh.material_data.nu = 0.3;
%! 	mesh.material_data.rho = 7850;
%! 	if (f_enable_constraint(k))
%! 	  opt_modes.shift_A = 0;
%! 	else
%! 	  opt_modes.shift_A = 1e-6;
%! 	endif
%! 	opt_modes.refine_max_iter = int32(10);
%! 	opt_modes.verbose = int32(0);
%! 	opt_modes.rigid_body_modes = interfaces{j};
%! 	[mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig, mat_ass_press, dof_map_press] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%! 	dof_map = fem_ass_dof_map(mesh, load_case);
%! 	[mat_ass.K, ...
%! 	 mat_ass.M, ...
%! 	 mat_ass.R, ...
%! 	 mat_info, ...
%! 	 mesh_info] = fem_ass_matrix(mesh, ...
%! 				     dof_map, ...
%! 				     [FEM_MAT_STIFFNESS, ...
%! 				      FEM_MAT_MASS, ...
%! 				      FEM_VEC_LOAD_CONSISTENT], ...
%! 				     load_case_bearing);
%! 	sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! 	[mesh, ...
%! 	 mat_ass, ...
%! 	 dof_map, ...
%! 	 sol_eig_cms, ...
%! 	 cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%! 	comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%!       endfor
%!     endfor
%!   endfor
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

%!test
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   interfaces = {"rigid", "flexible"};
%!   num_modes_cms = int32([0, 10]);
%!   for k=1:numel(num_modes_cms)
%!     for j=1:numel(interfaces)
%!       clear Fred Ritf bearing_surf cms_opt comp_mat dof_map_comb dof_map_post dof_map_press err_red
%!       clear grp_id_p1 grp_id_p2 grp_idx_p1 grp_idx_p2 load_case load_case_bearing load_case_itf
%!       clear load_case_post mat_ass_post mat_ass_press mesh mesh_comb mesh_data mesh_post mesh_size
%!       clear opt_modes p1 p1red p2 p2red pid qred sol_comb sol_eig sol_eig_cms sol_post sol_red tol_red
%!       unwind_protect
%! 	[fd, msg] = fopen([filename, ".geo"], "wt");
%! 	if (fd == -1)
%! 	  error("failed to open file \"%s.geo\"", filename);
%! 	endif
%! 	d = 14e-3;
%! 	D = 19.5e-3;
%! 	w = 5e-3;
%! 	l = 47e-3;
%! 	h = 5e-3;
%! 	grp_id_p1 = 2;
%! 	grp_id_p2 = 3;
%! 	p1 = 0;
%! 	p2 = 1;
%! 	scale_def = 5e-3;
%! 	mesh_size = 7e-3;
%! 	num_modes = 100;
%! 	fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%! 	fprintf(fd, "d = %g;\n", d);
%! 	fprintf(fd, "D = %g;\n", D);
%! 	fprintf(fd, "w = %g;\n", w);
%! 	fprintf(fd, "l = %g;\n", l);
%! 	fprintf(fd, "h = %g;\n", h);
%! 	fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%! 	fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%! 	fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%! 	fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%! 	fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%! 	fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%! 	fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%! 	fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%! 	fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%! 	fputs(fd, "Circle(1) = {2, 1, 3};\n");
%! 	fputs(fd, "Circle(2) = {3, 1, 4};\n");
%! 	fputs(fd, "Circle(3) = {4, 1, 5};\n");
%! 	fputs(fd, "Circle(4) = {5, 1, 2};\n");
%! 	fputs(fd, "Circle(5) = {7, 6, 8};\n");
%! 	fputs(fd, "Circle(6) = {8, 6, 9};\n");
%! 	fputs(fd, "Circle(7) = {9, 6, 10};\n");
%! 	fputs(fd, "Circle(8) = {10, 6, 7};\n");
%! 	fputs(fd, "Circle(9) = {11, 1, 12};\n");
%! 	fputs(fd, "Circle(10) = {12, 1, 13};\n");
%! 	fputs(fd, "Line(11) = {13, 14};\n");
%! 	fputs(fd, "Circle(12) = {14, 6, 15};\n");
%! 	fputs(fd, "Circle(13) = {15, 6, 16};\n");
%! 	fputs(fd, "Line(14) = {16, 11};\n");
%! 	fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%! 	fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%! 	fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%! 	fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%! 	fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%! 	fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%! 	fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%! 	fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%! 	fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!       unwind_protect_cleanup
%! 	if (fd ~= -1)
%! 	  fclose(fd);
%! 	  fd = -1;
%! 	endif
%!       end_unwind_protect
%!       fprintf(stderr, "meshing ...\n");
%!       pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%! 	   error("gmsh failed with status %d", status);
%!       endif
%!       fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!       mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!       fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!       grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!       grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!       cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!       cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!       bearing_surf(1).group_idx = grp_idx_p1;
%!       bearing_surf(1).options.reference_pressure = 1e9;
%!       bearing_surf(1).options.mesh_size = 1e-3;
%!       bearing_surf(1).options.include_rigid_body_modes = true;
%!       bearing_surf(1).options.bearing_type = "shell";
%!       bearing_surf(1).options.matrix_type = "modal substruct total";
%!       bearing_surf(1).r = 0.5 * d;
%!       bearing_surf(1).w = w;
%!       bearing_surf(1).X0 = [l; 0; 0];
%!       bearing_surf(1).R = eye(3);
%!       bearing_surf(1).relative_tolerance = 0;
%!       bearing_surf(1).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!       bearing_surf(1).options.number_of_modes = num_modes;
%!       bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!       bearing_surf(2).group_idx = grp_idx_p2;
%!       bearing_surf(2).options.reference_pressure = 1e9;
%!       bearing_surf(2).options.mesh_size = 1e-3;
%!       bearing_surf(2).options.include_rigid_body_modes = false;
%!       bearing_surf(2).options.bearing_type = "shell";
%!       bearing_surf(2).options.matrix_type = "modal substruct total";
%!       bearing_surf(2).r = 0.5 * d;
%!       bearing_surf(2).w = w;
%!       bearing_surf(2).X0 = [0; 0; 0];
%!       bearing_surf(2).R = eye(3);
%!       bearing_surf(2).relative_tolerance = 0;
%!       bearing_surf(2).absolute_tolerance = sqrt(eps) * 0.5 * d;
%!       bearing_surf(2).options.number_of_modes = num_modes;
%!       bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!       mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!       mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!       mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!       mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!       cms_opt.inveriants = true;
%!       cms_opt.modes.number = num_modes_cms(k);
%!       cms_opt.static_modes = false;
%!       cms_opt.modal_node_constraint = false;
%!       cms_opt.load_cases = "index";
%!       cms_opt.refine_max_iter = int32(10);
%!       load_case(1).locked_dof = false(rows(mesh.nodes), 6);
%!       mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!       mesh.material_data.E = 210000e6;
%!       mesh.material_data.nu = 0.3;
%!       mesh.material_data.rho = 7850;
%!       opt_modes.shift_A = 1e-6;
%!       opt_modes.refine_max_iter = int32(10);
%!       opt_modes.verbose = int32(0);
%!       opt_modes.rigid_body_modes = interfaces{j};
%!       [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig, mat_ass_press, dof_map_press] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!       [mesh, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!       comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%!       load_case_itf = fem_pre_load_case_create_empty(6);
%!       for i=1:6
%! 	   load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%! 	   load_case_itf(i).loads = zeros(1, 6);
%! 	   load_case_itf(i).loads(i) = 1;
%!       endfor
%!       Ritf = fem_ass_matrix(mesh, dof_map, FEM_VEC_LOAD_CONSISTENT, load_case_itf);
%!       nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!       nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!       nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!       nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!       p1red = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!       p2red = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!       Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red / bearing_surf(1).options.reference_pressure + comp_mat(2).E(:, 1:end - nz2) * p2red / bearing_surf(2).options.reference_pressure;
%!       Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!       qred = mat_ass.Kred \ Fred;
%!       sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!       mesh_post = mesh;
%!       mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!       load_case_post = fem_pre_load_case_create_empty(7);
%!       load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!       for i=1:6
%! 	   load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%! 	   load_case_post(i).loads = zeros(1, 6);
%! 	   load_case_post(i).loads(i) = 1;
%!       endfor
%!       load_case_post(7).pressure.tria6.elements = [mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%! 						   mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :)];
%!       load_case_post(7).pressure.tria6.p = [repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%! 					    repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6)];
%!       mesh_post.elements.joints.C = eye(6);
%!       mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!       dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!       [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!       sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post);
%!       mesh_data(1).mesh = mesh;
%!       mesh_data(1).dof_map = dof_map;
%!       mesh_data(2).mesh = mesh_post;
%!       mesh_data(2).dof_map = dof_map_post;
%!       [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!       sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!       sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!       sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!       err_red = zeros(1, size(sol_post.def, 3));
%!       for i=1:size(sol_post.def, 3)
%! 	   err_red(i) = max(max(abs(sol_post.def(1:end - 2, :, i) - sol_red.def(1:end - 2, :, i)))) / max(max(abs(sol_post.def(1:end - 2, :, i))));
%!       endfor
%!       for i=1:size(sol_post.def, 3)
%! 	   fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!       endfor
%!       tol_red = 3e-2;
%!       assert(all(err_red < tol_red));
%!     endfor
%!   endfor
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
