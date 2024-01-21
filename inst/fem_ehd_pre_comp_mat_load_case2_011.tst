## fem_ehd_pre_comp_mat_load_case2.m:11
%!test
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
%!        opt_modes.number_of_threads = mbdyn_solver_num_threads_default();
%!        cms_opt.solver = opt_modes.solver;
%!        cms_opt.number_of_threads = opt_modes.number_of_threads;
%!        [mesh, load_case_bearing, bearing_surf, cms_opt.load_cases_index, sol_eig] = fem_ehd_pre_comp_mat_load_case2(mesh, load_case, bearing_surf, opt_modes);
%!        [mesh, mat_ass, dof_map, sol_eig_cms, cms_opt] = fem_cms_create(mesh, load_case_bearing, cms_opt);
%!        assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!        assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
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
%!        sol_post_modal = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, cms_opt.modes.number + 6, 0, sqrt(eps), "shift-invert", opt_modes.solver);
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
%!        assert_simple(all(err_red < tol_red(j)));
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
