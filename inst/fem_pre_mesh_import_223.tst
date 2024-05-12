## fem_pre_mesh_import.m:223
%!test
%! try
%! ### TEST 223
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5
%! ####################################################
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
%!     unit_watt = unit_newton * unit_meters / unit_second;
%!     lambda = 10e-3 / unit_meters;
%!     deltaPML = 3 * lambda;
%!     r0 = 100e-3 / unit_meters;
%!     r1 = r0 + 3 * lambda;
%!     r2 = r1 + deltaPML;
%!     dx = lambda / 30;
%!     rho = 1000 / (unit_kilograms / unit_meters^3);
%!     c = 1440 / (unit_meters / unit_second);
%!     nPML = 1;
%!     alphaPML = 1;
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     f_use_PML = true;
%!     f_use_impedance = false;
%!     f_enable_plot = true;
%!     Aref = (1 + 0j) / (unit_pascal * unit_meters);
%!     pref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ r; ## according equation 5.11.6
%!     zref = @(r) rho * c * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according equation 5.11.9
%!     vnref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ (r * zref(r)); ## according equation 5.11.17
%!     alpha = atan2(2 * dx, r0);
%!     beta = atan2(2 * dx, r0);
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "r0=%g;\n", r0);
%!     fprintf(fd, "r1=%g;\n", r1);
%!     fprintf(fd, "r2=%g;\n", r2);
%!     fprintf(fd, "rref = (r0 + r2) / 2;\n");
%!     fprintf(fd, "alpha=%g;\n", alpha);
%!     fprintf(fd, "beta=%g;\n", beta);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {r0*Cos(alpha/2),r0*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(2) = {r1*Cos(alpha/2),r1*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(3) = {r2*Cos(alpha/2),r2*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "tmp1[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{1};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "tmp2[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{2};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Coherence;\n");
%!     fputs(fd, "tmp4[] = Extrude{{0,0,1},{0,0,0},-alpha/2} {\n");
%!     fputs(fd, "  Surface{tmp1[1],tmp2[1]};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{tmp4[1],tmp4[7]};}} = dx;\n");
%!     fputs(fd, "Physical Volume(\"v1\",1) = {tmp4[1]};\n");
%!     fputs(fd, "Physical Volume(\"v2\",2) = {tmp4[7]};\n");
%!     fputs(fd, "Physical Surface(\"s1\", 4) = {tmp4[2]};\n");
%!     fputs(fd, "Physical Surface(\"s2\", 5) = {tmp4[3]};\n");
%!     fputs(fd, "Physical Surface(\"s3\", 6) = {tmp4[9]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_s1 = find([mesh.groups.tria6h.id] == 4);
%!   grp_idx_s2 = find([mesh.groups.tria6h.id] == 5);
%!   grp_idx_s3 = find([mesh.groups.tria6h.id] == 6);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.elements.particle_velocity.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s1).elements, :);
%!   mesh.materials.particle_velocity.tria6h = ones(rows(mesh.elements.particle_velocity.tria6h.nodes), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s2).elements, :);
%!   mesh.materials.acoustic_boundary.tria6h = ones(rows(mesh.elements.acoustic_boundary.tria6h), 1, "int32");
%!   if (f_use_impedance)
%!     mesh.elements.acoustic_impedance.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s3).elements, :);
%!     mesh.elements.acoustic_impedance.tria6h.z = repmat(zref(r2), size(mesh.elements.acoustic_impedance.tria6h.nodes));
%!     mesh.materials.acoustic_impedance.tria6h = ones(rows(mesh.elements.acoustic_impedance.tria6h.nodes), 1, "int32");
%!   endif
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.tria6h.vn = repmat(-real(vnref(r0,0)), size(mesh.elements.particle_velocity.tria6h.nodes));
%!   load_case(2).particle_velocity.tria6h.vn = repmat(-imag(vnref(r0,0)), size(mesh.elements.particle_velocity.tria6h.nodes));
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   if (f_use_PML)
%!     load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s3).nodes, :) = true;
%!   endif
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   if (f_use_PML)
%!     mesh.elements.perfectly_matched_layers.tet10h.f = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!     mesh.elements.perfectly_matched_layers.tet10h.e1 = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!     mesh.elements.perfectly_matched_layers.tet10h.e2 = zeros(3, columns(mesh.elements.tet10h), rows(mesh.elements.tet10h));
%!     for j=1:3
%!       mesh.elements.perfectly_matched_layers.tet10h.f(j, :, :) = 1;
%!     endfor
%!     xi = mesh.nodes(:, 1)(mesh.elements.tet10h.');
%!     yi = mesh.nodes(:, 2)(mesh.elements.tet10h.');
%!     zi = mesh.nodes(:, 3)(mesh.elements.tet10h.');
%!     ri = sqrt(xi.^2 + yi.^2 + zi.^2) - r1;
%!     mesh.elements.perfectly_matched_layers.tet10h.e1(1, :, :) = xi;
%!     mesh.elements.perfectly_matched_layers.tet10h.e1(2, :, :) = yi;
%!     mesh.elements.perfectly_matched_layers.tet10h.e1(3, :, :) = zi;
%!     mesh.elements.perfectly_matched_layers.tet10h.e2(2, :, :) = 1;
%!   endif
%!   P = zeros(1, numel(alphaPML));
%!   T = zeros(1, numel(alphaPML));
%!   for i=1:numel(alphaPML)
%!     if (f_use_PML)
%!       sigmax = (ri >= 0) .* (1 ./ (1 - ri / deltaPML).^nPML - 1) * (2 * pi * alphaPML(i) * c / deltaPML);
%!       mesh.elements.perfectly_matched_layers.tet10h.f(1, :, :) = 1 ./ (1 - 1j * sigmax / omega);
%!     endif
%!     [mat_ass.Ka_re, ...
%!      mat_ass.Ka_im, ...
%!      mat_ass.Da_re, ...
%!      mat_ass.Da_im, ...
%!      mat_ass.Ma_re, ...
%!      mat_ass.Ma_im, ...
%!      mat_ass.Ra, ...
%!      mat_ass.mat_info, ...
%!      mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                           FEM_MAT_STIFFNESS_ACOUSTICS_IM, ...
%!                                           FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                           FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                           FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                           FEM_MAT_MASS_ACOUSTICS_IM, ...
%!                                           FEM_VEC_LOAD_ACOUSTICS], ...
%!                                          load_case);
%!     opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!     opt_sol.solver = "pastix";
%!     opt_sol.refine_max_iter = int32(50);
%!     Keff = -omega^2 * complex(mat_ass.Ma_re, mat_ass.Ma_im) + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + complex(mat_ass.Ka_re, mat_ass.Ka_im);
%!     Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!     Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!     Psi = linspace(0, 2 * pi, 37);
%!     idx = dof_map.ndof(:, 1);
%!     iact = find(idx > 0);
%!     solC.Phi = solC.PhiP = zeros(rows(dof_map.ndof), columns(Phi));
%!     solC.Phi(iact, :) = Phi(idx(iact), :);
%!     solC.PhiP(iact, :) = 1j * omega * Phi(idx(iact), :);
%!     [solC.particle_velocity, ...
%!      solC.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                                dof_map, ...
%!                                                [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                                 FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                                load_case, ...
%!                                                solC);
%!
%!     p = -solC.PhiP(mesh.elements.acoustic_boundary.tria6h);
%!     vx = solC.particle_velocity.vn.tria6h;
%!     R = abs(mean(mean((p - zref(r1) * vx) ./ (p + zref(r1) * vx))));
%!     T(i) = 1 / (1 - R^2);
%!     P(i) = sum(solC.acoustic_intensity.P.tria6h);
%!     P0 = 1e-12 / unit_watt;
%!     fprintf(stderr, "alphaPML=%.2f LW=%.2fdB TL=%.2fdB\n", alphaPML(i), 10 * log10(P(i)/P0), 10 * log10(T(i)));
%!     solR.t = Psi / omega;
%!     solR.p = real(-solC.PhiP * exp(1j * Psi));
%!     rg = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!     [rg, idx] = sort(rg);
%!     if (f_enable_plot)
%!       for j=1:columns(solR.p)
%!         figure("visible", "off");
%!         hold on;
%!         plot(rg * unit_meters, solR.p(idx, j) * unit_pascal, "-;p(r);r");
%!         plot(rg * unit_meters, real(pref(rg, solR.t(j))) * unit_pascal, "-;pref(r);k");
%!         xlabel("x [m]");
%!         ylabel("p [Pa]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("alpha=%.2f pressure t=%.2fs", alphaPML(i), solR.t(j)));
%!         ylim([-1, 1] * max(max(abs(solR.p))) * unit_pascal);
%!       endfor
%!     endif
%!     idx2 = find(rg < r1);
%!     rg = rg(idx2);
%!     idx = idx(idx2);
%!     tol = 1e-2;
%!     for j=1:numel(solR.t)
%!       prefj = real(pref(rg, solR.t(j)));
%!       assert_simple(solR.p(idx, j), prefj, tol * max(abs(prefj)));
%!     endfor
%!   endfor
%!   alphaPMLopt = alphaPML(find(T == min(T)));
%!   fprintf(stderr, "optimum alphaPML=%.2f\n", alphaPMLopt);
%!   figure("visible", "off");
%!   plot(alphaPML, 10 * log10(P / P0), "-x;P;r");
%!   grid minor on;
%!   xlabel("alphaPML [1]");
%!   ylabel("P [dB]");
%!   title("sound power level");
%!   figure("visible", "off");
%!   plot(alphaPML, 10 * log10(T), "-x;TL;r");
%!   grid minor on;
%!   xlabel("alphaPML [1]");
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
