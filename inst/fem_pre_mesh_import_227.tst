## fem_pre_mesh_import.m:227
%!test
%! ### TEST 227
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
%!     unit_meters = 1;
%!     unit_second = 1;
%!     unit_kilograms = 1;
%!     unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!     unit_pascal = unit_newton / unit_meters^2;
%!     unit_watt = unit_newton * unit_meters / unit_second;
%!     lambda = 1 / unit_meters;
%!     c = 1 / (unit_meters / unit_second);
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
%!     deltaPML = 3 * lambda;
%!     r0 = lambda / unit_meters;
%!     r1 = r0 + 2 * lambda;
%!     r2 = r1 + deltaPML;
%!     dx = lambda / 50;
%!     rho = 1 / (unit_kilograms / unit_meters^3);
%!     alphaPML = 1;
%!     f_use_PML = true;
%!     f_use_impedance = false;
%!     f_enable_plot = true;
%!     Aref = (1 + 0j) / (unit_pascal * unit_meters);
%!     pref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ r; ## according equation 5.11.6
%!     zref = @(r) rho * c * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according equation 5.11.9
%!     vnref = @(r, t) Aref * exp(1j * (omega * t - k * r)) ./ (r * zref(r)); ## according equation 5.11.17
%!     alpha = atan2(dx, r0);
%!     beta = atan2(dx, r0);
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "r0=%.16g;\n", r0);
%!     fprintf(fd, "r1=%.16g;\n", r1);
%!     fprintf(fd, "r2=%.16g;\n", r2);
%!     fprintf(fd, "rref = (r0 + r1) / 2;\n");
%!     fprintf(fd, "alpha=%g;\n", alpha);
%!     fprintf(fd, "beta=%g;\n", beta);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {r0*Cos(alpha/2),r0*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(2) = {r1*Cos(alpha/2),r1*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(3) = {r2*Cos(alpha/2),r2*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "tmp1[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{1}; Layers{Ceil(beta/2 * rref / dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "tmp2[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{2}; Layers{Ceil(beta/2 * rref / dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Coherence;\n");
%!     fputs(fd, "tmp4[] = Extrude{{0,0,1},{0,0,0},-alpha/2} {\n");
%!     fputs(fd, "  Surface{tmp1[1],tmp2[1]}; Layers{Ceil(alpha/2 * rref / dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
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
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v1 = find([mesh.groups.iso8.id] == 1);
%!   grp_idx_v2 = find([mesh.groups.iso8.id] == 2);
%!   elem_idx_v2 = mesh.groups.iso8(grp_idx_v2).elements;
%!   grp_idx_s1 = find([mesh.groups.iso4.id] == 4);
%!   grp_idx_s2 = find([mesh.groups.iso4.id] == 5);
%!   grp_idx_s3 = find([mesh.groups.iso4.id] == 6);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.elements.particle_velocity.iso4.nodes = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s1).elements, :);
%!   mesh.materials.particle_velocity.iso4 = ones(rows(mesh.elements.particle_velocity.iso4.nodes), 1, "int32");
%!   mesh.elements.acoustic_boundary.iso4 = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s2).elements, :);
%!   mesh.materials.acoustic_boundary.iso4 = ones(rows(mesh.elements.acoustic_boundary.iso4), 1, "int32");
%!   if (f_use_impedance)
%!     mesh.elements.acoustic_impedance.iso4.nodes = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s3).elements, :);
%!     mesh.elements.acoustic_impedance.iso4.z = repmat(zref(r2), size(mesh.elements.acoustic_impedance.iso4.nodes));
%!     mesh.materials.acoustic_impedance.iso4 = ones(rows(mesh.elements.acoustic_impedance.iso4.nodes), 1, "int32");
%!   endif
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.iso4.vn = repmat(-real(vnref(r0,0)), size(mesh.elements.particle_velocity.iso4.nodes));
%!   load_case(2).particle_velocity.iso4.vn = repmat(-imag(vnref(r0,0)), size(mesh.elements.particle_velocity.iso4.nodes));
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   if (f_use_PML)
%!     load_case_dof.locked_dof(mesh.groups.iso4(grp_idx_s3).nodes, :) = true;
%!   endif
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   if (f_use_PML)
%!     mesh.elements.perfectly_matched_layers.iso8.f = zeros(3, columns(mesh.elements.iso8), rows(mesh.elements.iso8));
%!     for j=1:3
%!       mesh.elements.perfectly_matched_layers.iso8.f(j, :, :) = 1;
%!     endfor
%!     mesh.elements.perfectly_matched_layers.iso8.e1 = zeros(3, columns(mesh.elements.iso8), rows(mesh.elements.iso8));
%!     mesh.elements.perfectly_matched_layers.iso8.e2 = zeros(3, columns(mesh.elements.iso8), rows(mesh.elements.iso8));
%!     xi = mesh.nodes(:, 1)(mesh.elements.iso8.');
%!     yi = mesh.nodes(:, 2)(mesh.elements.iso8.');
%!     zi = mesh.nodes(:, 3)(mesh.elements.iso8.');
%!     ri = sqrt(xi.^2 + yi.^2 + zi.^2) - r1;
%!     mesh.elements.perfectly_matched_layers.iso8.e1(1, :, :) = xi;
%!     mesh.elements.perfectly_matched_layers.iso8.e1(2, :, :) = yi;
%!     mesh.elements.perfectly_matched_layers.iso8.e1(3, :, :) = zi;
%!     mesh.elements.perfectly_matched_layers.iso8.e2(2, :, :) = 1;
%!   endif
%!   P = zeros(1, numel(alphaPML));
%!   T = zeros(1, numel(alphaPML));
%!   err = zeros(1, numel(alphaPML));
%!   for i=1:numel(alphaPML)
%!     if (f_use_PML)
%!       ## A parameter-free perfectly matched layer formulation
%!       ## for the finite-element-based solution of the Helmholtz
%!       ## equation
%!       ## Radu Cimpeanu a, Anton Martinsson b , Matthias Heil c
%!       ## a Department
%!       ## of Mathematics, Imperial College London, SW7 2AZ, London, United
%!       ## Kingdom
%!       ## (https://core.ac.uk/download/pdf/77019151.pdf)
%!       sigmax = alphaPML(i) ./ (deltaPML - ri(:, elem_idx_v2));
%!       mesh.elements.perfectly_matched_layers.iso8.f(1, :, elem_idx_v2) = 1 ./ (1 - 1j * sigmax / k);
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
%!     opt_sol.verbose = int32(0);
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
%!     p = -solC.PhiP(mesh.elements.acoustic_boundary.iso4);
%!     vx = solC.particle_velocity.vn.iso4;
%!     R = abs(mean(mean((p - zref(r1) * vx) ./ (p + zref(r1) * vx))));
%!     T(i) = 1 / (1 - R^2);
%!     P(i) = sum(solC.acoustic_intensity.P.iso4);
%!     P0 = 1e-12 / unit_watt;
%!     solR.t = Psi / omega;
%!     solR.p = real(-solC.PhiP * exp(1j * Psi));
%!     rg = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!     [rg, idx] = sort(rg);
%!     if (f_enable_plot)
%!       for j=1:columns(solR.p)
%!         figure("visible", "off");
%!         hold on;
%!         plot(rg * unit_meters, solR.p(idx, j) * unit_pascal, "-;p(r);1");
%!         plot(rg * unit_meters, real(pref(rg, solR.t(j))) * unit_pascal, "-;pref(r);0");
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
%!     tol = 5e-3;
%!     for j=1:numel(solR.t)
%!       prefj = real(pref(rg, solR.t(j)));
%!       err(i) = max([err(i), max(abs(solR.p(idx, j) - prefj)) / max(abs(prefj))]);
%!       assert_simple(solR.p(idx, j), prefj, tol * max(abs(prefj)));
%!     endfor
%!     fprintf(stderr, "alphaPML=%.2f LW=%.2fdB TL=%.2fdB err=%.2e\n", alphaPML(i), 10 * log10(P(i)/P0), 10 * log10(T(i)), err(i));
%!   endfor
%!   fprintf(stderr, "optimum alphaPML=%.2f based on reflections\n", alphaPML(find(T == min(T))));
%!   fprintf(stderr, "optimum alphaPML=%.2f based on analytical solution\n", alphaPML(find(err == min(err))));
%!   figure("visible", "off");
%!   plot(alphaPML, 10 * log10(P / P0), "-x;P;1");
%!   grid minor on;
%!   xlabel("alphaPML [1]");
%!   ylabel("P [dB]");
%!   title("sound power level");
%!   figure("visible", "off");
%!   plot(alphaPML, 10 * log10(T), "-x;TL;1");
%!   grid minor on;
%!   xlabel("alphaPML [1]");
%!   ylabel("TL [dB]");
%!   title("transmission loss");
%!   figure("visible", "off");
%!   hold on;
%!   plot(alphaPML, 100 * err, "-;err(alphaPML);1");
%!   xlabel("alphaPML [1]");
%!   ylabel("err [%]");
%!   grid minor on;
%!   title("numerical error");
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
