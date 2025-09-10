## fem_ehd_pre_comp_mat_load_case2.m:07
%!test
%! try
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
%!    num_modes = int32([0, 5]);
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
%!            ro = 18e-3;
%!            h = 28e-3;
%!            c = 10e-3;
%!            b = h - 2 * c;
%!            scale_def = 5e-3;
%!            mesh_size = 100e-3;
%!            fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!            fprintf(fd, "ri = %g;\n", ri);
%!            fprintf(fd, "ro = %g;\n", ro);
%!            fprintf(fd, "h = %g;\n", h);
%!            fprintf(fd, "c = %g;\n", c);
%!            fprintf(fd, "m = %g;\n", mesh_size);
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
%!            fputs(fd, "MeshSize{PointsOf{Volume{tmp[1]};}} = m;\n");
%!            fputs(fd, "Mesh.ElementOrder = 3;\n");
%!            fputs(fd, "Mesh.OptimizeThreshold = 0.99;\n");
%!          unwind_protect_cleanup
%!            if (fd ~= -1)
%!              fclose(fd);
%!              fd = -1;
%!            endif
%!          end_unwind_protect
%!          fprintf(stderr, "meshing ...\n");
%!          pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!          status = spawn_wait(pid);
%!          if (status ~= 0)
%!            error("gmsh failed with status %d", status);
%!          endif
%!          fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!          mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!          fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!          cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!          cms_opt.solver = "pardiso";
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
%!          assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!          assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
