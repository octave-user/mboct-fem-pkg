## fem_pre_mesh_import.m:381
%!test
%! try
%! ## TEST 381
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   animate = false;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     E = 200000e6;
%!     nu = 0.3;
%!     rho = 7850;
%!     R = 100e-3;
%!     h = 10e-3;
%!     t = 2e-3;
%!     N = 10;
%!     scale_def = 0.5 * R;
%!     Phi = 0 * pi / 180;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "R = %g;\n", R);
%!     fprintf(fd, "t = %g;\n", t);
%!     fprintf(fd, "Phi = %g;\n", Phi);
%!     fputs(fd, "Point(1) = {  0.0, 0.0,    0.0, h};\n");
%!     fputs(fd, "Point(2) = {  R * Sin(Phi), 0.0,     -R * Cos(Phi), h};\n");
%!     fputs(fd, "Point(3) = {    R, 0.0,    0.0, h};\n");
%!     fputs(fd, "Point(4) = {  R * Sin(Phi), 0.0,      R * Cos(Phi), h};\n");
%!     fputs(fd, "Point(5) = { (R + t) * Sin(Phi), 0.0,  (R + t) * Cos(Phi), h};\n");
%!     fputs(fd, "Point(6) = { R + t, 0.0,    0.0, h};\n");
%!     fputs(fd, "Point(7) = {  (R + t) * Sin(Phi), 0.0, -(R + t) * Cos(Phi), h};\n");
%!     fputs(fd, "Circle(1) = {2,1,3};\n");
%!     fputs(fd, "Circle(2) = {3,1,4};\n");
%!     fputs(fd, "Line(3) = {4,5};\n");
%!     fputs(fd, "Circle(4) = {5, 1, 6};\n");
%!     fputs(fd, "Circle(5) = {6, 1, 7};\n");
%!     fputs(fd, "Line(6) = {7, 2};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fprintf(fd, "tmp[] = Extrude {{0, 0, 1},{0,0,0}, Pi/2}{ Surface{6}; };\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {6, tmp[0]};\n");
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
%!   load_case.locked_dof(mesh.groups.tria6.nodes, 1:3) = true;
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.M, ...
%!    mat_ass.K] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_STIFFNESS], ...
%!                                load_case);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N);
%!   for i=1:size(sol_eig.def, 3)
%!     sol_eig.def(:, :, i) /= max(max(abs(sol_eig.def(:, 1:3, i))));
%!   endfor
%!   if (animate)
%!     opt_anim.scale_def = scale_def;
%!     opt_anim.animation_delay = 1;
%!     opt_anim.print_and_exit = true;
%!     opt_anim.print_to_file = filename;
%!     opt_anim.rotation_angle = [90, 0, 0] * pi / 180;
%!     opt_anim.skin_only = true;
%!     opt_anim.show_element = false;
%!     unwind_protect
%!       fem_post_sol_external(mesh, sol_eig, opt_anim);
%!       [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title("Gmsh - deformed mesh / continuous stress tensor");
%!     unwind_protect_cleanup
%!       unlink([opt_anim.print_to_file, "_001.jpg"]);
%!     end_unwind_protect
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
