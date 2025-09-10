## fem_pre_mesh_extrude_surf.m:04
%!test
%! try
%! ### TEST 4
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
%!     lambda = 20e-3 / unit_meters;
%!     deltaPML = 1e-4 * lambda;
%!     r0 = 10e-3 / unit_meters;
%!     r1 = r0 + lambda / 2;
%!     r2 = r1 + 50 * lambda;
%!     r3 = r2 + deltaPML;
%!     dx = lambda / 10;
%!     rho = 1000 / (unit_kilograms / unit_meters^3);
%!     c = 1440 / (unit_meters / unit_second);
%!     layersPML = max([1, ceil(deltaPML / dx)]);
%!     f = c / lambda;
%!     omega = 2 * pi * f;
%!     k = omega / c;
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
%!     fprintf(fd, "alpha=%g;\n", alpha);
%!     fprintf(fd, "beta=%g;\n", beta);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {r0*Cos(alpha/2),r0*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Point(2) = {r1*Cos(alpha/2),r1*Sin(alpha/2),0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "tmp1[] = Extrude{{-Sin(alpha/2),Cos(alpha/2),0},{0,0,0},-beta/2} {\n");
%!     fputs(fd, "  Line{1};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "tmp4[] = Extrude{{0,0,1},{0,0,0},-alpha/2} {\n");
%!     fputs(fd, "  Surface{tmp1[1]};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"v1\",1) = {tmp4[1]};\n");
%!     fputs(fd, "Physical Surface(\"s1\", 2) = {2};\n");
%!     fputs(fd, "Physical Surface(\"s2\", 3) = {3};\n");
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
%!   N = ceil((r2 - r1) / dx);
%!   [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "tria6h", 3, [repmat((r2 - r1) / N, 1, N), r3 - r2]);
%!   r = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!   mesh.groups.tria6h(end + 1).id = 4;
%!   mesh.groups.tria6h(end).nodes = find(r >= r3 - deltaPML / 2);
%!   grp_idx_s1 = find([mesh.groups.tria6h.id] == 2);
%!   grp_idx_s2 = find([mesh.groups.tria6h.id] == 3);
%!   grp_idx_s3 = find([mesh.groups.tria6h.id] == 4);
%!   node_idx_constr = mesh.groups.tria6h(grp_idx_s1).nodes;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_s2).elements, :);
%!   mesh.materials.acoustic_boundary.tria6h = ones(rows(mesh.elements.acoustic_boundary.tria6h), 1, "int32");
%!   mesh.elements.acoustic_constr = struct("C", mat2cell(ones(1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))), ...
%!                                          "nodes", mat2cell(node_idx_constr, 1, ones(1, numel(node_idx_constr))), ...
%!                                          "scale", mat2cell(repmat(1/omega, 1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))));
%!   load_case = struct("acoustic_constr", cell(1, 2));
%!   p_constr = repmat(pref(r0, 0), size(node_idx_constr));
%!   load_case(1).acoustic_constr = struct("p", mat2cell(real(p_constr), 1, ones(1, numel(p_constr))));
%!   load_case(2).acoustic_constr = struct("p", mat2cell(imag(p_constr), 1, ones(1, numel(p_constr))));
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_s3).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   mesh.elements.perfectly_matched_layers.penta15.f = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   mesh.elements.perfectly_matched_layers.penta15.e1 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   mesh.elements.perfectly_matched_layers.penta15.e2 = zeros(3, columns(mesh.elements.penta15), rows(mesh.elements.penta15));
%!   for j=1:3
%!     mesh.elements.perfectly_matched_layers.penta15.f(j, :, :) = 1;
%!   endfor
%!   xi = mesh.nodes(:, 1)(mesh.elements.penta15.');
%!   yi = mesh.nodes(:, 2)(mesh.elements.penta15.');
%!   zi = mesh.nodes(:, 3)(mesh.elements.penta15.');
%!   vi = sqrt(xi.^2 + yi.^2 + zi.^2) - r2;
%!   [ridx, cidx] = find(vi >= 0);
%!   mesh.elements.perfectly_matched_layers.penta15.e1(1, :, :) = xi;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(2, :, :) = yi;
%!   mesh.elements.perfectly_matched_layers.penta15.e1(3, :, :) = zi;
%!   mesh.elements.perfectly_matched_layers.penta15.e2(2, :, :) = 1;
%!   mesh.elements.perfectly_matched_layers.penta15.f(1, ridx, cidx) = 1j * k * (deltaPML - vi(ridx, cidx));
%!   [mat_ass.Ka_re, ...
%!    mat_ass.Ka_im, ...
%!    mat_ass.Da_re, ...
%!    mat_ass.Da_im, ...
%!    mat_ass.Ma_re, ...
%!    mat_ass.Ma_im, ...
%!    mat_ass.Ra, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_STIFFNESS_ACOUSTICS_IM, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_IM, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pardiso";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * complex(mat_ass.Ma_re, mat_ass.Ma_im) + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + complex(mat_ass.Ka_re, mat_ass.Ka_im);
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   idx = dof_map.ndof(:, 1);
%!   iact = find(idx > 0);
%!   solC.Phi = solC.PhiP = zeros(rows(mesh.nodes), columns(Phi));
%!   solC.Phi(iact, :) = Phi(idx(iact), :);
%!   solC.PhiP(iact, :) = 1j * omega * Phi(idx(iact), :);
%!   [solC.particle_velocity, ...
%!    solC.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                              dof_map, ...
%!                                              [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                               FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                              load_case, ...
%!                                              solC);
%!   p = -solC.PhiP(mesh.elements.acoustic_boundary.tria6h);
%!   vx = solC.particle_velocity.vn.tria6h;
%!   R = abs(mean(mean((p - zref(r1) * vx) ./ (p + zref(r1) * vx))));
%!   T = 1 / (1 - R^2);
%!   P = sum(solC.acoustic_intensity.P.tria6h);
%!   P0 = 1e-12 / unit_watt;
%!   fprintf(stderr, "LW=%.2fdB TL=%.2fdB\n", 10 * log10(P/P0), 10 * log10(T));
%!   solR.t = Psi / omega;
%!   solR.p = real(-solC.PhiP * exp(1j * Psi));
%!   rg = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2 + mesh.nodes(:, 3).^2);
%!   [rg, idx] = sort(rg);
%!   if (f_enable_plot)
%!     for j=1:columns(solR.p)
%!       figure("visible", "off");
%!       hold on;
%!       plot(rg * unit_meters, solR.p(idx, j) * unit_pascal, "-;p(r);r");
%!       plot(rg * unit_meters, real(pref(rg, solR.t(j))) * unit_pascal, "-;pref(r);k");
%!       xlabel("x [m]");
%!       ylabel("p [Pa]");
%!       grid on;
%!       grid minor on;
%!       title(sprintf("pressure t=%.2fs", solR.t(j)));
%!       ylim([-1, 1] * max(max(abs(solR.p))) * unit_pascal);
%!     endfor
%!   endif
%!   idx2 = find(rg < r1);
%!   rg = rg(idx2);
%!   idx = idx(idx2);
%!   tol = 2e-2;
%!   pref = real(pref(rg, solR.t));
%!   assert_simple(solR.p(idx, :), pref, tol * max(max(abs(pref))));
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
