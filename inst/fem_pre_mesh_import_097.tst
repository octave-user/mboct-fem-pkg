## fem_pre_mesh_import.m:97
%!test
%! try
%! ## TEST 97: thermal strain tet10h
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
%!     d1 = 50e-3;
%!     d2 = 150e-3;
%!     L = 0.5e-3;
%!     h = 0.5e-3;
%!     Phi = 5 * pi / 180;
%!     alpha1 = 50;
%!     alpha2 = 100;
%!     lambda = 40;
%!     Theta2 = 100;
%!     q1 = 1000000;
%!     q2 = q1 * d1 / d2;
%!     Q1 = d1 * pi * L * q1;
%!     kR = pi / (1 / (alpha1 * d1) + 1 / (2 * lambda) * log(d2 / d1) + 1 / (alpha2 * d2));
%!     Theta1 = Theta2 + Q1 / (L * kR);
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "d1 = %g;\n", d1);
%!     fprintf(fd, "d2 = %g;\n", d2);
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "Phi = %g;\n", Phi);
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0, h};\n");
%!     fputs(fd, "Point(2) = {0.5 * d1, 0, 0, h};\n");
%!     fputs(fd, "Point(3) = {0.5 * d1 * Cos(Phi), 0.5 * d1 * Sin(Phi), 0, h};\n");
%!     fputs(fd, "Point(4) = {0.5 * d2 * Cos(Phi), 0.5 * d2 * Sin(Phi), 0, h};\n");
%!     fputs(fd, "Point(5) = {0.5 * d2, 0, 0, h};\n");
%!     fputs(fd, "Circle(1) = {2,1,3};\n");
%!     fputs(fd, "Line(2) = {3,4};\n");
%!     fputs(fd, "Circle(3) = {4,1,5};\n");
%!     fputs(fd, "Line(4) = {5,2};\n");
%!     fputs(fd, "Curve Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {0, 0, L}{ Surface{6}; };\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"inside-surface\", 2) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"outside-surface\", 3) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"top-surface\", 4) = {tmp[0]};\n");
%!     fputs(fd, "Physical Surface(\"bottom-surface\", 5) = {6};\n");
%!     fputs(fd, "Physical Surface(\"side-surface1\", 6) = {tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"side-surface2\", 7) = {tmp[4]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=3;\n");
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
%!   opt_mesh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.material_data.gamma = 1.26e-5;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.cp = 465;
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.k = eye(3) * lambda;
%!   grp_idx_inside = find([mesh.groups.tria6h.id] == 2);
%!   grp_idx_outside = find([mesh.groups.tria6h.id] == 3);
%!   load_case.heat_source.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_inside).elements, :);
%!   load_case.heat_source.tria6h.q = repmat(q1, numel(mesh.groups.tria6h(grp_idx_inside).elements), 6);
%!   mesh.elements.convection.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_outside).elements, :);
%!   mesh.elements.convection.tria6h.h = repmat(alpha2, numel(mesh.groups.tria6h(grp_idx_outside).elements), 6);
%!   load_case.convection.tria6h.theta = repmat(Theta2, numel(mesh.groups.tria6h(grp_idx_outside).elements), 6);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.Kk, ...
%!    mat_ass.Q] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_THERMAL_COND, ...
%!                                 FEM_VEC_LOAD_THERMAL], ...
%!                                load_case);
%!   sol.theta = fem_sol_factor(mat_ass.Kk) \ mat_ass.Q;
%!   idx_inside = mesh.groups.tria6h(grp_idx_inside).nodes;
%!   idx_outside = mesh.groups.tria6h(grp_idx_outside).nodes;
%!   Theta1s = sol.theta(idx_inside, :);
%!   Theta2s = sol.theta(idx_outside, :);
%!   Theta1a = Theta1s + q1 / alpha1;
%!   Theta2a = Theta2s - q2 / alpha2;
%!   tol = eps^0.5;
%!   assert_simple(Theta1a, repmat(Theta1, size(Theta1a)), tol * abs(Theta1 - Theta2));
%!   assert_simple(Theta2a, repmat(Theta2, size(Theta2a)), tol * abs(Theta1 - Theta2));
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
