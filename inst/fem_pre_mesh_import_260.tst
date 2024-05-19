## fem_pre_mesh_import.m:260
%!test
%! try
%! ## TEST260
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
%!     ri = 8e-3;
%!     ro = 10e-3;
%!     h = 12e-3;
%!     c = 2e-3;
%!     b = h - 2 * c;
%!     p1 = 25.79e6;
%!     p2 = 7.83e6;
%!     p3 = 1.3758e6;
%!     scale_def = 5e-3;
%!     mesh_size = 3e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "ri = %g;\n", ri);
%!     fprintf(fd, "ro = %g;\n", ro);
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "c = %g;\n", c);
%!     fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!     fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!     fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!     fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!     fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!     fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!     fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!     fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,5};\n");
%!     fputs(fd, "Line(5) = {5,6};\n");
%!     fputs(fd, "Line(6) = {6,7};\n");
%!     fputs(fd, "Line(7) = {7,8};\n");
%!     fputs(fd, "Line(8) = {8,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!     fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria10].id] == 1);
%!   load_case.locked_dof(mesh.groups.tria10(grp_id_clamp).nodes, :) = true;
%!   grp_id_p1 = find([[mesh.groups.tria10].id] == 3);
%!   grp_id_p2 = find([[mesh.groups.tria10].id] == 2);
%!   grp_id_p3 = find([[mesh.groups.tria10].id] == 4);
%!   elem_id_p1 = mesh.groups.tria10(grp_id_p1).elements;
%!   elem_id_p2 = mesh.groups.tria10(grp_id_p2).elements;
%!   elem_id_p3 = mesh.groups.tria10(grp_id_p3).elements;
%!   elno_p1 = mesh.elements.tria10(elem_id_p1, :);
%!   elno_p2 = mesh.elements.tria10(elem_id_p2, :);
%!   elno_p3 = mesh.elements.tria10(elem_id_p3, :);
%!   load_case.pressure.tria10.elements = [elno_p1; elno_p2; elno_p3];
%!   load_case.pressure.tria10.p = [repmat(p1, rows(elno_p1), columns(elno_p1));
%!                                 repmat(p2, rows(elno_p2), columns(elno_p2));
%!                                 repmat(p3, rows(elno_p3), columns(elno_p3))];
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   if (do_plot)
%!     figure("visible", "off");
%!     fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("static deflection - consistent pressure load");
%!   endif
%!   X1 = mesh.nodes(unique(elno_p1), 1:3).';
%!   X2 = mesh.nodes(unique(elno_p2), 1:3).';
%!   X3 = mesh.nodes(unique(elno_p3), 1:3).';
%!   dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%!   dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%!   dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%!   F1_con = full(mat_ass.R(dof1)).';
%!   F2_con = full(mat_ass.R(dof2)).';
%!   F3_con = full(mat_ass.R(dof3)).';
%!   M1_con = cross(X1, F1_con);
%!   M2_con = cross(X2, F2_con);
%!   M3_con = cross(X3, F3_con);
%!   Ftot1_con = sum(F1_con, 2);
%!   Ftot2_con = sum(F2_con, 2);
%!   Ftot3_con = sum(F3_con, 2);
%!   Mtot1_con = sum(M1_con, 2);
%!   Mtot2_con = sum(M2_con, 2);
%!   Mtot3_con = sum(M3_con, 2);
%!   ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);
%!   F1_an = [ri * b * p1;
%!            ri * b * p1;
%!            0];
%!   M1_an = [-ri * b * p1 * (c + b/2);
%!            ri * b * p1 * (c + b/2);
%!            0];
%!   F2_an = [-ro * b * p2;
%!            -ro * b * p2;
%!            0];
%!   M2_an = [ ro * b * p2 * (c + b/2);
%!             -ro * b * p2 * (c + b/2);
%!             0];
%!   F3_an = [0;
%!            0;
%!            -p3 * (ro^2 - ri^2) * pi / 4];
%!   M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            0];
%!   assert_simple(Ftot1_con, F1_an, eps^0.9 * norm(F1_an));
%!   assert_simple(Ftot2_con, F2_an, eps^0.9 * norm(F2_an));
%!   assert_simple(Mtot1_con, M1_an, eps^0.9 * norm(M1_an));
%!   assert_simple(Mtot2_con, M2_an, eps^0.9 * norm(M2_an));
%!   assert_simple(Ftot3_con, F3_an, eps^0.3 * norm(F3_an));
%!   assert_simple(Mtot3_con, M3_an, eps^0.3 * norm(M3_an));
%!   if (do_plot)
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
