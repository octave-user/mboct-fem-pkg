## fem_pre_mesh_import.m:453
%!test
%! try
%! ## TEST 453: dynamic patch test of penta18
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
%!     a = 8e-3;
%!     b = 15e-3;
%!     c = 12e-3;
%!     mesh_size = 7e-3;
%!     scale_def = 1e-3;
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
%!     fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; Layers{Ceil(c / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!     fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(1).nodes, 1:3) = false;
%!   mesh.materials.penta18 = ones(rows(mesh.elements.penta18), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.dm, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS, ...
%!                                         FEM_MAT_MASS, ...
%!                                         FEM_MAT_MASS_LUMPED, ...
%!                                         FEM_SCA_TOT_MASS], ...
%!                                       load_case);
%!   opt_modal.rho = -(max(max(abs(mat_ass.K))) / max(max(abs(mat_ass.M))))^0.5;
%!   opt_modal.solver = "pardiso";
%!   opt_modal.refine_max_iter = int32(250);
%!   opt_modal.pre_scaling = true;
%!   opt_modal.verbose = false;
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 10, opt_modal);
%!   sol_eig.def *= scale_def / max(max(max(abs(sol_eig.def))));
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig);
%!   mat_ass.M = mat_ass.Mdiag;
%!   sol_eig2 = fem_sol_modal(mesh, dof_map, mat_ass, 10, opt_modal);
%!   sol_eig2.def *= scale_def / max(max(max(abs(sol_eig2.def))));
%!   [sol_eig2.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig2);
%!   dmref = a * b * c * mesh.material_data.rho;
%!   assert_simple(mat_ass.dm, dmref, eps * dmref);
%!   assert_simple(sum(diag(mat_ass.Mdiag))/3, dmref, eps^0.9 * dmref);
%!   tolt = eps^0.5;
%!   assert_simple(max(max(max(max(abs(sol_eig.stress.tau.penta18(:, :, :, 1:6)))))) < tolt * max(max(max(max(abs(sol_eig.stress.tau.penta18(:, :, :, 7:end)))))));
%!   assert_simple(max(max(max(max(abs(sol_eig2.stress.tau.penta18(:, :, :, 1:6)))))) < tolt * max(max(max(max(abs(sol_eig2.stress.tau.penta18(:, :, :, 7:end)))))));
%!   tolf = eps^0.4;
%!   assert_simple(all(sol_eig.f(1:6) < tolf * max(sol_eig.f(7:10))));
%!   assert_simple(all(sol_eig2.f(1:6) < tolf * max(sol_eig2.f(7:10))));
%!   assert_simple(sum(sum(mat_ass.M)) / 3, mat_ass.dm, 1e-10 * mat_ass.dm);
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
