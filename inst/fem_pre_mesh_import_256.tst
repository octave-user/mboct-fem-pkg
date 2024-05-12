## fem_pre_mesh_import.m:256
%!test
%! try
%! ## TEST 256: dynamic patch test of tet20
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     a = 10;
%!     b = 10;
%!     c = 10;
%!     mesh_size = 10;
%!     N = 20;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a = %g;\n", a);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "c = %g;\n", c);
%!     fprintf(fd, "h = %g;\n", mesh_size);
%!     fputs(fd, "Point(1) = {0,0,0,h};\n");
%!     fputs(fd, "Point(2) = {a,0,0,h};\n");
%!     fputs(fd, "Point(3) = {a,b,0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0,h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; };\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!     fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=1;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   E = 210000;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850e-9;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   dof_map.parallel.threshold_elem = int32(1000);
%!   [mat_ass.M, ...
%!    mat_ass.D, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mtot, ...
%!    mat_ass.J, ...
%!    mat_ass.S] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_DAMPING, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_SCA_TOT_MASS, ...
%!                                 FEM_MAT_INERTIA_J, ...
%!                                 FEM_VEC_INERTIA_M1], ...
%!                                load_case);
%!   opt_sol.solver = "pastix";
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.verbose = int32(0);
%!   opt_sol.pre_scaling = true;
%!   opt_sol.refine_max_iter = int32(250);
%!   opt_sol.rho = -(max(max(abs(mat_ass.K))) / max(max(abs(mat_ass.M))))^0.5;
%!   opt_sol.tolerance = sqrt(eps);
%!   opt_sol.algorithm = "shift-invert";
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, opt_sol);
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   mref = a * b * c * mesh.material_data.rho;
%!   lcg = [a; b; c] / 2;
%!   J1ref = mref * (b^2 + c^2) / 12;
%!   J2ref = mref * (a^2 + c^2) / 12;
%!   J3ref = mref * (a^2 + b^2) / 12;
%!   Jref = diag([J1ref, J2ref, J3ref]) - skew(lcg) * skew(lcg) * mref;
%!   Sref = mref * lcg;
%!   assert_simple(mat_ass.mtot, mref, eps^0.8 * mref);
%!   assert_simple(mat_ass.J, Jref, eps^0.8 * norm(Jref));
%!   assert_simple(mat_ass.S, Sref, eps^0.8 * norm(Sref));
%!   tolf = eps^0.3;
%!   tolt = eps^0.4;
%!   assert_simple(max(max(max(max(abs(sol_eig.stress.tau.tet20(:, :, :, 1:6)))))) < tolt * max(max(max(max(abs(sol_eig.stress.tau.tet20(:, :, :, 7:end)))))));
%!   assert_simple(all(sol_eig.f(1:6) < tolf * max(sol_eig.f(7:10))));
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
