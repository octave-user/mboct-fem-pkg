## fem_ehd_pre_comp_mat_load_case2.m:04
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
%!          assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
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
%!    assert_simple(all(err_red(:, end) < tol_red));
%!    assert_simple(err_mod(end) < tol_mod);
%!    assert_simple(err_w(end) < tol_w);
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
