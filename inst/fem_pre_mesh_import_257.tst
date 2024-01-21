## fem_pre_mesh_import.m:257
%!test
%! ## TEST 257
%! ## M. Maeder, R. D'Auria, E. Grasso, G. Petrone b, S. De Rosa, M. Klaerner, L. Kroll, S. Marburg
%! ## Numerical analysis of sound radiation from rotating discs
%! ## Journal of Sound and Vibration 468 (2020) 115085
%! close all;
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
%!     unit_meters = 1e-3;
%!     unit_second = 1e4;
%!     unit_kilograms = 1e-3;
%!     unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!     unit_pascal = unit_newton / unit_meters^2;
%!     N = 13;
%!     d1 = 800e-3 / unit_meters;
%!     d2 = 120e-3 / unit_meters;
%!     d3 = 60e-3 / unit_meters;
%!     t = 3.5e-3 / unit_meters;
%!     h1 = 10 * t;
%!     E1 = 210000e6 / unit_pascal;
%!     rho1 = 7800 / (unit_kilograms / unit_meters^3);
%!     nu1 = 0.3;
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "d1 = %g;\n", d1);
%!     fprintf(fd, "d2 = %g;\n", d2);
%!     fprintf(fd, "d3 = %g;\n", d3);
%!     fprintf(fd, "t = %g;\n", t);
%!     fprintf(fd, "h1 = %g;\n", h1);
%!     fputs(fd, "v1 = newv;\n");
%!     fputs(fd, "Cylinder(v1) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d1};\n");
%!     fputs(fd, "v2 = newv;\n");
%!     fputs(fd, "Cylinder(v2) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d2};\n");
%!     fputs(fd, "v3 = newv;\n");
%!     fputs(fd, "Cylinder(v3) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d3};\n");
%!     fputs(fd, "v5 = BooleanFragments{Volume{v1};Delete;}{Volume{v2,v3};Delete;};\n");
%!     fputs(fd, "Physical Volume(\"solid\", 1) = {8, 9};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 3) = {14, 15};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 1;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{8, 9}; } } = h1;\n");
%!     fputs(fd, "Mesh.HighOrderIterMax = 100;\n");
%!     fputs(fd, "Mesh.HighOrderPassMax = 100;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5",  [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v_solid = find([mesh.groups.tet20.id] == 1);
%!   grp_idx_s_clamp = find([mesh.groups.tria10.id] == 3);
%!   mesh.material_data.E = E1;
%!   mesh.material_data.rho = rho1;
%!   mesh.material_data.nu = nu1;
%!   mesh.materials.tet20 = zeros(rows(mesh.elements.tet20), 1, "int32");
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_v_solid).elements) = 1;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.tria10(grp_idx_s_clamp).nodes, 1:3) = true;
%!   load_case_dof.domain = FEM_DO_STRUCTURAL;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS, ...
%!                                         FEM_MAT_MASS], ...
%!                                        load_case_dof);
%!   rho = 0;
%!   alg = "shift-invert";
%!   solver = "pastix";
%!   num_threads = mbdyn_solver_num_threads_default();
%!   tol_modes = 1e-3;
%!   fref = [22.1, 22.1, 25.6, 32.6, 32.6, 68.6, 68.6, 119.2, 119.2, 155.8, 167.8, 167.8, 182.9];
%!   tol_freq_rel = 0.5e-2;
%!   tol_freq_abs = 0.5;
%!   tol_freq = tol_freq_abs + tol_freq_rel * fref;
%!   sol = fem_sol_modal(mesh, dof_map, mat_ass, N, rho, tol_modes, alg, solver, num_threads);
%!   assert_simple(all(abs(sol.f * unit_second^-1 - fref) < tol_freq));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
