## fem_pre_mesh_import.m:94
%!test
%! try
%! ## TEST 94: thermal strain tet10
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
%!     dT = 400;
%!     mesh_size = 7e-3;
%!     scale_def = 10e-3;
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
%!   opt_mesh.elem_type = {"tet10", "tria6"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.tria6].id] == 4);
%!   grp_id_clamp_y = find([[mesh.groups.tria6].id] == 1);
%!   grp_id_clamp_z = find([[mesh.groups.tria6].id] == 5);
%!   if (~isempty(grp_id_clamp_x))
%!     load_case.locked_dof(mesh.groups.tria6(grp_id_clamp_x).nodes, 1) = true;
%!   endif
%!   if (~isempty(grp_id_clamp_y))
%!     load_case.locked_dof(mesh.groups.tria6(grp_id_clamp_y).nodes, 2) = true;
%!   endif
%!   if (~isempty(grp_id_clamp_z))
%!     load_case.locked_dof(mesh.groups.tria6(grp_id_clamp_z).nodes, 3) = true;
%!   endif
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.gamma = 1.26e-5;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   load_case.dTheta = repmat(dT, rows(mesh.nodes), 1);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress, ...
%!    sol_stat.strain] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH, ...
%!                                       FEM_VEC_STRAIN_TOTAL], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   load_case2.epsilon0 = sol_stat.strain.epsilon;
%!   [mat_ass2.K, ...
%!    mat_ass2.R] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS, ...
%!                                  FEM_VEC_LOAD_CONSISTENT], ...
%!                                 load_case2);
%!   [sol_stat2] = fem_sol_static(mesh, dof_map, mat_ass2);
%!   [sol_stat2.stress, ...
%!    sol_stat2.strain] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_VEC_STRESS_CAUCH, ...
%!                                        FEM_VEC_STRAIN_TOTAL], ...
%!                                       load_case2, ...
%!                                       sol_stat2);
%!   tol = eps^0.7;
%!   epsilon_th = dT * mesh.material_data.gamma;
%!   for i=1:3
%!     assert_simple(sol_stat.def(:, i), mesh.nodes(:, i) * dT * mesh.material_data.gamma, tol * max(abs([a, b, c] * epsilon_th)));
%!     for j=1:columns(sol_stat.strain.epsilon.tet10)
%!       assert_simple(sol_stat.strain.epsilon.tet10(:, j, i), repmat(dT * mesh.material_data.gamma, rows(sol_stat.strain.epsilon.tet10), 1), tol * abs(epsilon_th));
%!       assert_simple(sol_stat.strain.epsilon.tet10(:, j, i + 3), zeros(rows(sol_stat.strain.epsilon.tet10), 1), tol * abs(epsilon_th));
%!     endfor
%!   endfor
%!   for i=1:3
%!     assert_simple(sol_stat2.def(:, i), mesh.nodes(:, i) * dT * mesh.material_data.gamma, tol * max(abs([a, b, c] * epsilon_th)));
%!     for j=1:columns(sol_stat2.strain.epsilon.tet10)
%!       assert_simple(sol_stat2.strain.epsilon.tet10(:, j, i), repmat(dT * mesh.material_data.gamma, rows(sol_stat2.strain.epsilon.tet10), 1), tol * abs(epsilon_th));
%!       assert_simple(sol_stat2.strain.epsilon.tet10(:, j, i + 3), zeros(rows(sol_stat2.strain.epsilon.tet10), 1), tol * abs(epsilon_th));
%!     endfor
%!   endfor
%!   assert_simple(sol_stat2.def, sol_stat.def, tol * max(max(max(abs(sol_stat.def)))));
%!   assert_simple(sol_stat2.strain.epsilon.tet10, sol_stat.strain.epsilon.tet10, tol * max(max(max(abs(sol_stat.strain.epsilon.tet10)))));
%!   assert_simple(sol_stat2.strain.epsilonm.tet10, sol_stat.strain.epsilonm.tet10, tol * max(max(max(abs(sol_stat.strain.epsilonm.tet10)))));
%!   assert_simple(max(max(max(abs(sol_stat.stress.tau.tet10)))) < tol * abs(E * epsilon_th));
%!   assert_simple(max(max(max(abs(sol_stat2.stress.tau.tet10)))) < tol * abs(E * epsilon_th));
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
