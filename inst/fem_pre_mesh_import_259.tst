## fem_pre_mesh_import.m:259
%!test
%! try
%! ## TEST 259: static patch test of tet20
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
%!     px = 25.79e6;
%!     py = 7.83e6;
%!     pz = 1.3758e6;
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
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opts.elem_type = {"tet20", "tria10"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opts));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.tria10].id] == 4);
%!   grp_id_clamp_y = find([[mesh.groups.tria10].id] == 1);
%!   grp_id_clamp_z = find([[mesh.groups.tria10].id] == 5);
%!   load_case.locked_dof(mesh.groups.tria10(grp_id_clamp_x).nodes, 1) = true;
%!   load_case.locked_dof(mesh.groups.tria10(grp_id_clamp_y).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.tria10(grp_id_clamp_z).nodes, 3) = true;
%!   grp_id_px = find([[mesh.groups.tria10].id] == 2);
%!   grp_id_py = find([[mesh.groups.tria10].id] == 3);
%!   grp_id_pz = find([[mesh.groups.tria10].id] == 6);
%!   elem_id_px = mesh.groups.tria10(grp_id_px).elements;
%!   elem_id_py = mesh.groups.tria10(grp_id_py).elements;
%!   elem_id_pz = mesh.groups.tria10(grp_id_pz).elements;
%!   elno_px = mesh.elements.tria10(elem_id_px, :);
%!   elno_py = mesh.elements.tria10(elem_id_py, :);
%!   elno_pz = mesh.elements.tria10(elem_id_pz, :);
%!   load_case.pressure.tria10.elements = [elno_px; elno_py; elno_pz];
%!   load_case.pressure.tria10.p = [repmat(px, rows(elno_px), columns(elno_px));
%!                                 repmat(py, rows(elno_py), columns(elno_py));
%!                                 repmat(pz, rows(elno_pz), columns(elno_pz))];
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                 load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   if (do_plot)
%!     figure("visible", "off");
%!     fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("static deflection - consistent pressure load");
%!     figure_list();
%!   endif
%!   tol = eps^0.7;
%!   assert_simple(max(max(max(abs(sol_stat.stress.tau.tet20(:, :, 1) / -px - 1)))) < tol);
%!   assert_simple(max(max(max(abs(sol_stat.stress.tau.tet20(:, :, 2) / -py - 1)))) < tol);
%!   assert_simple(max(max(max(abs(sol_stat.stress.tau.tet20(:, :, 3) / -pz - 1)))) < tol);
%!   assert_simple(max(max(max(abs(sol_stat.stress.tau.tet20(:, :, 4:6) / max(abs([px,py,pz])))))) < tol);
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
