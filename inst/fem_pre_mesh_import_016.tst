## fem_pre_mesh_import.m:16
%!test
%! try
%! ## TEST 16
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
%!     ## K.J.Bathe page 328 4.20a
%!     mesh_size = 1.5;
%!     p1 = 0.006;
%!     E = 55;
%!     nu = 0.3;
%!     rho = 1000e-12;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%!     fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%!     fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%!     fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%!     fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%!     fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%!     fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%!     fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%!     fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%!     fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%!     fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%!     fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,11,4};\n");
%!     fputs(fd, "Line(4) = {4,5};\n");
%!     fputs(fd, "Line(5) = {5,6};\n");
%!     fputs(fd, "Line(6) = {6,7};\n");
%!     fputs(fd, "Circle(7) = {7,12,8};\n");
%!     fputs(fd, "Line(8) = {8,9};\n");
%!     fputs(fd, "Line(9) = {9,10};\n");
%!     fputs(fd, "Line(10) = {10,1};\n");
%!     fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%!     fputs(fd, "Plane Surface(14) = {11};\n");
%!     fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; };\n", mesh_size);
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%!     fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%!     fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%!     fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%!     fputs(fd, "Mesh.HighOrderIterMax = 200;\n");
%!     fputs(fd, "Mesh.HighOrderPassMax = 100;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold = 0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%!   grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!   elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%!   elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%!   load_case.pressure.tria6.elements = elno_p1;
%!   load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_VEC_LOAD_LUMPED], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   grp_id_displacement = find([[mesh.groups.tria6].id] == 4);
%!   elem_id_displacement = mesh.groups.tria6(grp_id_displacement).elements;
%!   elno_id_displacement = mesh.elements.tria6(elem_id_displacement, :);
%!   delta = mean(sol_stat.def(elno_id_displacement, 3));
%!   grp_id_stress = find([[mesh.groups.tria6].id] == 5);
%!   elem_id_stress = mesh.groups.tria6(grp_id_stress).elements;
%!   elno_id_stress = mesh.elements.tria6(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(mesh.elements.tet10 == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(sol_stat.stress.taum.tet10(ridx(j), cidx(j), :), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   sigma1_max = 0;
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max = max(sigma1_max, max(eig(TAU)));
%!   endfor
%!   fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   ## K.J.Bathe page 329 4.20b
%!   sigma1_max_ref = 0.6056;
%!   delta_ref = -1.669;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   assert_simple(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%!   assert_simple(delta, delta_ref, 0.04 * abs(delta_ref));
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
