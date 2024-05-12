## fem_pre_mesh_unstruct_create.m:06
%!test
%! try
%! ## TEST6
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.2e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0.01e-3; ## FIXME: di=0 not supported by Gmsh 4.5.5
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.62;
%! B = 3.5;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.d^2 * pi / 4);
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
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * L, 0.0, 0.5 * d, h1};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Line(4) = {4,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,1};\n");
%!     fputs(fd, "Line Loop(7) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(8) = {7};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{8}; Layers{Ceil(d * Pi / 4 / h0)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{8,v[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {8};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[6]};\n");
%!   unwind_protect_cleanup
%!     fclose(fd);
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.quad8(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.quad8(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.quad8(3).nodes, 1) = true;
%!   load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(4).elements, :);
%!   load_case.pressure.quad8.p = repmat(p, rows(load_case.pressure.quad8.elements), 8);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.iso20)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
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
%!   assert_simple(Kt, Kt_a, 0.15 * Kt_a);
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
