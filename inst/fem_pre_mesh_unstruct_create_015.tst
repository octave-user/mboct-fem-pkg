## fem_pre_mesh_unstruct_create.m: 15
%!test
%! try
%! ## TEST 15
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h = 0.25e-3;
%! geo.Do = 4e-3;
%! geo.Di = 2e-3;
%! geo.L = 4e-3;
%! filename = "";
%! f_run_post_proc = false;
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     geo_file = [filename, ".geo"];
%!     fd = fopen(geo_file, "w");
%!     if (fd == -1)
%!       error("failed to open file %s", geo_file);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * Di, h};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * Do, h};\n");
%!     fputs(fd, "Point(3) = { 0.5 * L, 0.0, 0.5 * Do, h};\n");
%!     fputs(fd, "Point(4) = { 0.5 * L, 0.0, 0.5 * Di, h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{6}; Layers{Ceil(Do * Pi / 4 / h)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",0) = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\",1) = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\",2) = {6};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {v[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   opt.mesh.elem_type = {"tria6", "quad9", "penta18", "tet10h"};
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   elem_types = fieldnames(mesh.elements);
%!   mesh.materials = struct();
%!   for i=1:numel(elem_types)
%!     mesh.materials = setfield(mesh.materials, elem_types{i}, ones(rows(getfield(mesh.elements, elem_types{i})), 1, "int32"));
%!   endfor
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   node_id_disp = [];
%!   if (isfield(mesh.groups, "tria6"))
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==1)).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==2)).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==3)).nodes, 1) = true;
%!   node_id_disp = unique([node_id_disp, mesh.groups.tria6(find([mesh.groups.tria6.id]==4)).nodes]);
%!   endif
%!   if (isfield(mesh.groups, "quad9"))
%!   load_case.locked_dof(mesh.groups.quad9(find([mesh.groups.quad9.id]==3)).nodes, 1) = true;
%!   node_id_disp = unique([node_id_disp, mesh.groups.quad9(find([mesh.groups.quad9.id]==4)).nodes]);
%!   endif
%!   empty_cell = cell(1, numel(node_id_disp));
%!   mesh.elements.joints = struct("C", empty_cell, "nodes", empty_cell);
%!   load_case.joints = struct("U", empty_cell);
%!   for i=1:numel(node_id_disp)
%!     mesh.elements.joints(i).C = [1, 0, 0, 0, 0, 0];
%!     mesh.elements.joints(i).nodes = node_id_disp(i);
%!     load_case.joints(i).U = geo.L / E;
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mtot, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_SCA_TOT_MASS], ...
%!                                load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   m_a = (geo.Do^2 - geo.Di^2) * pi / 4 * geo.L / 4 * rho;
%!   assert_simple(mat_ass.mtot, m_a, eps^0.4 * m_a);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   if (f_run_post_proc)
%!     fem_post_sol_external(mesh, sol_stat, opt_post);
%!     figure("visible", "off");
%!     [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!     imshow(img, alpha);
%!     title("deformed mesh/van Mises stress");
%!     figure_list();
%!   endif
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
