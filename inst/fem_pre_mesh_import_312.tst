## fem_pre_mesh_import.m:312
%!test
%! try
%! ## TEST 312
%! ## M. Maeder, R. D'Auria, E. Grasso, G. Petrone b, S. De Rosa, M. Klaerner, L. Kroll, S. Marburg
%! ## Numerical analysis of sound radiation from rotating discs
%! ## Journal of Sound and Vibration 468 (2020) 115085
%! f_enable_plot = false;
%! if (f_enable_plot)
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
%!     unit_meters = 1e-3;
%!     unit_second = 1e4;
%!     unit_kilograms = 1e-3;
%!     unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!     unit_pascal = unit_newton / unit_meters^2;
%!     unit_watt = unit_newton * unit_meters / unit_second;
%!     P0 = 1e-12 / unit_watt;
%!     d1 = 800e-3 / unit_meters;
%!     d2 = 120e-3 / unit_meters;
%!     d3 = 60e-3 / unit_meters;
%!     t = 3.5e-3 / unit_meters;
%!     r = 1500e-3 / unit_meters;
%!     h1 = 20 * t;
%!     h2 = 400e-3 / unit_meters;
%!     h3 = 400e-3 / unit_meters;
%!     E1 = 210000e6 / unit_pascal;
%!     rho1 = 7800 / (unit_kilograms / unit_meters^3);
%!     nu1 = 0.3;
%!     alpha1 = 0.1826 / (unit_second^-1);
%!     beta1 = 5.0125e-6 / unit_second;
%!     c2 = 340 / (unit_meters / unit_second);
%!     rho2 = 1.225 / (unit_kilograms / unit_meters^3);
%!     Fz = (1 + 0j) / unit_newton;
%!     f = (100) / (unit_second^-1);
%!     deltaPML = 200e-3 / unit_meters;
%!     solver = "precond";
%!     f_enable_PML = true;
%!     num_threads = mbdyn_solver_num_threads_default();
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "d1 = %.16g;\n", d1);
%!     fprintf(fd, "d2 = %.16g;\n", d2);
%!     fprintf(fd, "d3 = %.16g;\n", d3);
%!     fprintf(fd, "t = %.16g;\n", t);
%!     fprintf(fd, "r = %.16g;\n", r);
%!     fprintf(fd, "deltaPML = %.16g;\n", deltaPML);
%!     fprintf(fd, "h1 = %.16g;\n", h1);
%!     fprintf(fd, "h2 = %.16g;\n", h2);
%!     fprintf(fd, "h3 = %.16g;\n", h3);
%!     fputs(fd, "v1 = newv;\n");
%!     fputs(fd, "Cylinder(v1) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d1};\n");
%!     fputs(fd, "v2 = newv;\n");
%!     fputs(fd, "Cylinder(v2) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d2};\n");
%!     fputs(fd, "v3 = newv;\n");
%!     fputs(fd, "Cylinder(v3) = {0, 0, -0.5 * t, 0, 0, t, 0.5 * d3};\n");
%!     fputs(fd, "v4 = newv;\n");
%!     fputs(fd, "Sphere(v4) = {0, 0, 0, r};\n");
%!     fputs(fd, "v5 = newv;\n");
%!     fputs(fd, "Sphere(v5) = {0, 0, 0, r + deltaPML};\n");
%!     fputs(fd, "v6 = BooleanFragments{Volume{v5};Delete;}{Volume{v1,v2,v3,v4};Delete;};\n");
%!     fputs(fd, "Physical Volume(\"solid\", 1) = {9, 10};\n");
%!     fputs(fd, "Physical Volume(\"fluid\", 2) = {7, 11};\n");
%!     fputs(fd, "Physical Volume(\"PML\", 3) = {8};\n");
%!     fputs(fd, "Physical Surface(\"fluid_struct\", 1) = {7, 12, 13, 14, 16, 17};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 2) = {16, 17};\n");
%!     fputs(fd, "Physical Surface(\"fluid-boundary\", 3) = {11};\n");
%!     fputs(fd, "Physical Surface(\"PML-boundary\", 4) = {10};\n");
%!     fputs(fd, "Physical Point(\"load\", 1) = {12};\n");
%!     fputs(fd, "ReverseMesh Surface{7,11,15};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 1;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{8}; } } = h3;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{7,11}; } } = h2;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{9,10}; } } = h1;\n");
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
%!   grp_idx_v_fluid = find([mesh.groups.tet20.id] == 2);
%!   grp_idx_v_PML = find([mesh.groups.tet20.id] == 3);
%!   grp_idx_s_fsi = find([mesh.groups.tria10.id] == 1);
%!   grp_idx_s_clamp = find([mesh.groups.tria10.id] == 2);
%!   grp_idx_s_boundary = find([mesh.groups.tria10.id] == 3);
%!   grp_idx_s_boundaryPML = find([mesh.groups.tria10.id] == 4);
%!   node_idx_constr = mesh.groups.tria10(grp_idx_s_boundaryPML).nodes;
%!   empty_cell = cell(1, 2);
%!   mesh.material_data = struct("E", empty_cell, ...
%!                               "nu", empty_cell, ...
%!                               "rho", empty_cell, ...
%!                               "c", empty_cell, ...
%!                               "alpha", empty_cell, ...
%!                               "beta", empty_cell);
%!   mesh.material_data(1).E = E1;
%!   mesh.material_data(1).rho = rho1;
%!   mesh.material_data(1).nu = nu1;
%!   mesh.material_data(1).alpha = alpha1;
%!   mesh.material_data(1).beta = beta1;
%!   mesh.material_data(2).c = c2;
%!   mesh.material_data(2).rho = rho2;
%!   mesh.materials.tet20 = zeros(rows(mesh.elements.tet20), 1, "int32");
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_v_solid).elements) = 1;
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_v_fluid).elements) = 2;
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_v_PML).elements) = 2;
%!   mesh.elements.fluid_struct_interface.tria10 = mesh.elements.tria10(mesh.groups.tria10(grp_idx_s_fsi).elements, :);
%!   mesh.elements.acoustic_boundary.tria10 = mesh.elements.tria10(mesh.groups.tria10(grp_idx_s_boundary).elements, :);
%!   mesh.materials.acoustic_boundary.tria10 = repmat(int32(2), rows(mesh.elements.acoustic_boundary.tria10), 1);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   load_case_dof.locked_dof(mesh.groups.tria10(grp_idx_s_clamp).nodes, 1:3) = true;
%!   if (f_enable_PML)
%!     load_case_dof.locked_dof(mesh.groups.tria10(grp_idx_s_boundaryPML).nodes, 7) = true;
%!   endif
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   empty_cell = cell(1, 2);
%!   load_case = struct("loads", empty_cell, "loaded_nodes", empty_cell);
%!   tol = sqrt(eps) * d1;
%!   node_id_load = find((abs(mesh.nodes(:, 1) - 0.5 * d1) < tol) & ...
%!                       (abs(mesh.nodes(:, 2)) < tol) & ...
%!                       (abs(mesh.nodes(:, 3) - 0.5 * t) < tol));
%!   node_id_load = node_id_load(1);
%!   load_case(1).loaded_nodes = node_id_load;
%!   load_case(1).loads = [zeros(1, 2), real(Fz), zeros(1, 3)];
%!   load_case(2).loaded_nodes = node_id_load;
%!   load_case(2).loads = [zeros(1, 2), imag(Fz), zeros(1, 3)];
%!
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = num_threads;
%!
%!   if (f_enable_PML)
%!     mesh.elements.perfectly_matched_layers.tet20.f = zeros(3, columns(mesh.elements.tet20), rows(mesh.elements.tet20));
%!     mesh.elements.perfectly_matched_layers.tet20.e1 = zeros(3, columns(mesh.elements.tet20), rows(mesh.elements.tet20));
%!     mesh.elements.perfectly_matched_layers.tet20.e2 = zeros(3, columns(mesh.elements.tet20), rows(mesh.elements.tet20));
%!     xi = mesh.nodes(:, 1)(mesh.elements.tet20.');
%!     yi = mesh.nodes(:, 2)(mesh.elements.tet20.');
%!     zi = mesh.nodes(:, 3)(mesh.elements.tet20.');
%!     elem_idx_PML = mesh.groups.tet20(grp_idx_v_PML).elements;
%!     ri = sqrt(xi(:, elem_idx_PML).^2 + yi(:, elem_idx_PML).^2 + zi(:, elem_idx_PML).^2) - r;
%!     mesh.elements.perfectly_matched_layers.tet20.e1(1, :, :) = xi;
%!     mesh.elements.perfectly_matched_layers.tet20.e1(2, :, :) = yi;
%!     mesh.elements.perfectly_matched_layers.tet20.e1(3, :, :) = zi;
%!     mesh.elements.perfectly_matched_layers.tet20.e2(1, :, :) = 1;
%!     for i=1:3
%!       mesh.elements.perfectly_matched_layers.tet20.f(i, :, :) = 1;
%!     endfor
%!   endif
%!   P = T = zeros(1, numel(f));
%!   ITER = repmat(intmax(), 1, 2);
%!   MAXITER = 20;
%!   TOL = eps^0.5;
%!   RESTART = 10;
%!   MAXITERS = int32(100);
%!   TOLS = eps^0.5;
%!   perm = [];
%!   Phi = zeros(dof_map.totdof, 1);
%!   cpu_factor = 0;
%!   cpu_solve = 0;
%!   FLAG = -1;
%!   alpha = 0.75;
%!   unwind_protect
%!     for i=1:numel(f)
%!       omega = 2 * pi * f(i);
%!       k = omega / c2;
%!       if (f_enable_PML)
%!         mesh.elements.perfectly_matched_layers.tet20.f(1, :, mesh.groups.tet20(grp_idx_v_PML).elements) = 1j * k * (deltaPML - ri);
%!       endif
%!       [mat_ass.Mfs_re, ...
%!        mat_ass.Mfs_im, ...
%!        mat_ass.Kfs_re, ...
%!        mat_ass.Kfs_im, ...
%!        mat_ass.Dfs_re, ...
%!        mat_ass.Dfs_im, ...
%!        mat_ass.Rfs, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                          FEM_MAT_MASS_FLUID_STRUCT_IM, ...
%!                                          FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                          FEM_MAT_STIFFNESS_FLUID_STRUCT_IM, ...
%!                                          FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                          FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                          FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                         load_case);
%!       Keff = -omega^2 * complex(mat_ass.Mfs_re, mat_ass.Mfs_im) ...
%!            + 1j * omega * complex(mat_ass.Dfs_re, mat_ass.Dfs_im) ...
%!            + complex(mat_ass.Kfs_re, mat_ass.Kfs_im);
%!       Reff = complex(mat_ass.Rfs(:, 1), mat_ass.Rfs(:, 2));
%!       opt_sol.number_of_threads = num_threads;
%!       opt_sol.solver = "pastix";
%!       opt_sol.compress_when = int32(0);
%!       opt_sol.compress_min_ratio = 1;
%!       opt_sol.compress_tolerance = 1e-2;
%!       opt_sol.verbose = int32(0);
%!       switch (solver)
%!         case "direct"
%!           opt_sol.refine_max_iter = MAXITER;
%!           opt_sol.epsilon_refinement = TOL;
%!           opt_sol.pre_scaling = true;
%!           Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!         case "precond"
%!           opt_sol.refine_max_iter = int32(0);
%!           opt_sol.pre_scaling = false;
%!           if (isempty(perm))
%!             perm = csymamd(Keff);
%!           endif
%!           do
%!             do_fact = (FLAG ~= 0 || cpu_solve > alpha * cpu_factor);
%!             if (do_fact)
%!               start = cputime();
%!               [KS, D1, D2] = fem_sol_matrix_scale(Keff(perm, perm), TOLS, MAXITERS);
%!               clear KSfact;
%!               KSfact = fem_sol_factor(KS, opt_sol);
%!               cpu_factor = cputime() - start;
%!             else
%!               KS = diag(D1) * Keff(perm, perm) * diag(D2);
%!             endif
%!             start = cputime();
%!             Z0 = KSfact \ (diag(D1) * Reff(perm));
%!             [Z, FLAG, RELRES, ITER] = gmres(@(Z) KSfact \ (KS * Z), Z0, RESTART, TOL, MAXITER, [], [], Z0);
%!             Phi(perm) = diag(D2) * Z;
%!             cpu_solve = cputime() - start;
%!           until (FLAG == 0 || do_fact)
%!         otherwise
%!           error("unknown solver: \"%s\"", solver);
%!       endswitch
%!       KeffPhi = Keff * Phi;
%!       f1 = norm(KeffPhi - Reff);
%!       f2 = norm(KeffPhi + Reff);
%!       sol.Phi = sol.PhiP = zeros(rows(mesh.nodes), 1);
%!       idx = dof_map.ndof(:, 7);
%!       iact = find(idx > 0);
%!       sol.Phi(iact) = Phi(idx(iact));
%!       sol.PhiP(iact) = 1j * omega * Phi(idx(iact));
%!       sol.p = -sol.PhiP;
%!       [sol.particle_velocity, ...
%!        sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                                 dof_map, ...
%!                                                 [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                                  FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                                 load_case, ...
%!                                                 sol);
%!       k = omega / c2;
%!       z = rho2 * c2 * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according Jont Allen equation 5.11.9
%!       p = sol.p(mesh.elements.acoustic_boundary.tria10);
%!       vx = sol.particle_velocity.vn.tria10;
%!       R = abs(mean(mean((p - z * vx) ./ (p + z * vx))));
%!       T(i) = 1 ./ (1 - R.^2);
%!       P(i) = sum(sol.acoustic_intensity.P.tria10);
%!       fprintf(stderr, "%d/%d %.2fHz P=%.2fdB TL=%.2fdB res=%.2e iter=%d/%d\n", i, numel(f), f(i) * unit_second^-1, 10 * log10(P(i) / P0), 10 * log10(T(i)), f1 / f2, ITER(1), ITER(2));
%!       solt.t = linspace(0, 2 * pi / omega, 36);
%!       solt.p = real(-sol.PhiP .* exp(1j * omega * solt.t));
%!       solt.def = zeros(rows(mesh.nodes), 6, numel(solt.t));
%!       solt.v = zeros(rows(mesh.nodes), 3, numel(solt.t));
%!       for j=1:6
%!         idx = dof_map.ndof(:, j);
%!         iact = find(idx > 0);
%!         for k=1:numel(solt.t)
%!           solt.def(iact, j, k) = real(Phi(idx(iact)) * exp(1j * omega * solt.t(k)));
%!         endfor
%!       endfor
%!       for j=1:columns(mesh.elements.tet20)
%!         for k=1:3
%!           for l=1:numel(solt.t)
%!             solt.v(mesh.elements.tet20(:, j), k, l) = real(sol.particle_velocity.v.tet20(:, j, k) .* exp(1j * omega * solt.t(l)));
%!           endfor
%!         endfor
%!       endfor
%!       if (f_enable_plot)
%!         Zc = linspace(0, r + deltaPML , 300);
%!         Xc = Yc = zeros(size(Zc));
%!         deltaC = 0.4 / unit_meters;
%!         delta0 = 0.1 / unit_meters;
%!         idxc = find((mesh.nodes(:, 3) >= 0) & (sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2) < abs(mesh.nodes(:, 3)) * (deltaC / (r + deltaPML)) + delta0));
%!         pc = zeros(numel(Zc), columns(solt.p));
%!         for j=1:columns(pc)
%!           pc(:, j) = griddata3(mesh.nodes(idxc, 1), mesh.nodes(idxc, 2), mesh.nodes(idxc, 3), solt.p(idxc, j), Xc, Yc, Zc);
%!         endfor
%!         for j=1:columns(pc)
%!           figure("visible","off");hold on
%!           plot(Zc * unit_meters, pc(:, j) * unit_pascal, "-;p(X=0,Y=0,Z);r");
%!           xlabel("r [m]");
%!           ylabel("p [Pa]");
%!           grid on;
%!           grid minor on;
%!           title(sprintf("sound pressure Phi=%.1fdeg f=%.1fHz", 180 / pi * (solt.t(i) * omega), f(i) * unit_second^-1));
%!           ylim([min(min(pc)),max(max(pc))] * unit_pascal);
%!         endfor
%!       endif
%!     endfor
%!   unwind_protect_cleanup
%!     clear KSfact;
%!   end_unwind_protect
%!   figure("visible", "off");
%!   hold on;
%!   plot(f * unit_second^-1, 10 * log10(P / P0), "-;P(f);r");
%!   xlabel("f [Hz]");
%!   ylabel("Lw [dB]");
%!   grid on;
%!   grid minor on;
%!   title("sound power level");
%!   figure("visible", "off");
%!   hold on;
%!   plot(f * unit_second^-1, 10 * log10(T), "-;TL(f);r");
%!   xlabel("f [Hz]");
%!   ylabel("TL [dB]");
%!   title("transmission loss");
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
