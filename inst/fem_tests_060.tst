## fem_tests.m:60
%!test
%! ## TEST 60
%! ## 24 TH INTERNATIONAL CONGRESS OF THE AERONAUTICAL SCIENCES
%! ## NATURAL FREQUENCIES OF ROTATING
%! ## CANTILEVER FLEXIBLE BEAMS BY USING THE
%! ## p-VERSION OF THE FINITE ELEMENT METHOD
%! ## Hamza cherif S. M.*, Houmat A.*
%! ## *Department of Mechanical Engineering , University of Tlemcen, Algeria
%! lambda_ref = 70; ## Table 2. page 6
%! mu_ref = [2,4,6,8,10,50];
%! omegax_ref_flapwise = [4.1373,5.5850,7.3603,9.2568,11.2023,51.0805];
%! omegax_ref_chordwise = [3.6195,3.888,4.2393,4.6105,4.97,7.3362];
%! close all;
%! rho = 2700;
%! E = 70000e6;
%! W = 10e-3;
%! H = 10e-3;
%! A = W * H;
%! I = W * H^3 / 12;
%! mu = 8;
%! lambda = lambda_ref;
%! L = lambda * sqrt(I/A);
%! material.E = 70000e6;
%! material.nu = 0.3;
%! material.rho = 2700;
%! num_modes = int32(10);
%! h = min(W, H)/4;
%! geometry.l = L;
%! geometry.w = W;
%! geometry.h = H;
%! r = linspace(0.2 * geometry.l, 0.9 * geometry.l, 100);
%! omegay = mu / sqrt(rho * A * L^4 / (E * I));
%! sigmaxx_ref = material.rho * omegay^2 * (geometry.l^2 - r.^2) / 2;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = zeros(3, 1);
%! [mesh, load_case_dof] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! mesh.nodes(:, 2) -= 0.5 * geometry.w;
%! mesh.nodes(:, 3) -= 0.5 * geometry.h;
%! xi = mesh.nodes(:, 1)(mesh.elements.iso8);
%! yi = mesh.nodes(:, 2)(mesh.elements.iso8);
%! zi = mesh.nodes(:, 3)(mesh.elements.iso8);
%! e1 = [1; 0.7; 0.9];
%! e2 = [0.8; -0.5; 0.6];
%! e3 = cross(e1, e2);
%! e2 = cross(e3, e1);
%! R1 = [e1, e2, e3];
%! R1 *= diag(1 ./ norm(R1, "cols"));
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! omega = R1 * [0; omegay; 0];
%! omegaq = [omega(1)^2;
%!           omega(2)^2;
%!           omega(3)^2;
%!           omega(1) * omega(2);
%!           omega(2) * omega(3);
%!           omega(3) * omega(1)];
%! load_case = fem_pre_load_case_create_empty(numel(omegaq));
%! for i=1:numel(load_case)
%!   load_case(i).omegaq = zeros(6, 1);
%!   load_case(i).omegaq(i) = 1;
%! endfor
%! load_case_tot.omega = omega;
%! mesh.nodes = [mesh.nodes(:, 1:3) * R1.', zeros(rows(mesh.nodes), 3)];
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! mat_ass_tot.M = mat_ass.M;
%! mat_ass_tot.K = mat_ass.K;
%! [mat_ass_tot.KOMEGA, ...
%!  mat_ass_tot.DOMEGA, ...
%!  mat_ass_tot.R] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS_OMEGA, ...
%!                                  FEM_MAT_DAMPING_OMEGA, ...
%!                                  FEM_VEC_LOAD_CONSISTENT], ...
%!                                  load_case_tot);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! sol_stat_tot = fem_sol_static(mesh, dof_map, mat_ass_tot);
%! sol_stat_tot.stress = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case_tot, ...
%!                                      sol_stat_tot);
%! load_case_pre_stress_tot.tau0 = sol_stat_tot.stress.tau;
%! mat_ass_tot.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress_tot);
%! for i=1:numel(load_case)
%!   load_case_pre_stress.omegaq = load_case(i).omegaq;
%!   load_case_pre_stress.tau0.iso8 = sol_stat.stress.tau.iso8(:, :, :, i);
%!   [KTAU0, KOMEGA] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0, ...
%!                                     FEM_MAT_STIFFNESS_OMEGA], ...
%!                                    load_case_pre_stress);
%!   KTAU0 *= omegaq(i);
%!   KOMEGA *= omegaq(i);
%!   if (i == 1)
%!     mat_ass.KTAU0 = KTAU0;
%!     mat_ass.KOMEGA = KOMEGA;
%!   else
%!     mat_ass.KTAU0 += KTAU0;
%!     mat_ass.KOMEGA += KOMEGA;
%!   endif
%! endfor
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
%! mat_ass.K += mat_ass.KTAU0 + mat_ass.KOMEGA;
%! sol_eig_preload = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
%! mat_ass_tot.K += mat_ass_tot.KTAU0 + mat_ass_tot.KOMEGA;
%! mat_ass_tot.D = mat_ass_tot.DOMEGA;
%! sol_eig_preload_tot = fem_sol_modal(mesh, dof_map, mat_ass_tot, num_modes);
%! sol_eig_rot = fem_sol_modal_damped(mesh, dof_map, mat_ass_tot, num_modes);
%! TAU = zeros(3, 3, rows(sol_stat.stress.taum.iso8), columns(sol_stat.stress.taum.iso8));
%! for k=1:size(TAU, 4)
%!   for j=1:size(TAU, 3)
%!     for i=1:3
%!       TAU(i, i, j, k) = sum(sol_stat.stress.taum.iso8(j, k, i, :)(:) .* omegaq);
%!     endfor
%!     TAU(1, 2, j, k) = TAU(2, 1, j, k) = sum(sol_stat.stress.taum.iso8(j, k, 4, :)(:) .* omegaq);
%!     TAU(2, 3, j, k) = TAU(3, 2, j, k) = sum(sol_stat.stress.taum.iso8(j, k, 5, :)(:) .* omegaq);
%!     TAU(1, 3, j, k) = TAU(3, 1, j, k) = sum(sol_stat.stress.taum.iso8(j, k, 6, :)(:) .* omegaq);
%!     TAU(:, :, j, k) = R1.' * TAU(:, :, j, k) * R1;
%!   endfor
%! endfor
%! sigmaxx = griddata3(xi, yi, zi, reshape(TAU(1, 1, :, :), size(xi)), r, zeros(size(r)), zeros(size(r)));
%! tol_tau = 0.2e-1;
%! assert_simple(sigmaxx(:), sigmaxx_ref(:), tol_tau * max(abs(sigmaxx_ref)));
%! tol_R = 1e-6;
%! assert_simple(mat_ass.R * omegaq, mat_ass_tot.R, tol_R * norm(mat_ass_tot.R));
%! tol_f = 1e-6;
%! assert_simple(sol_eig_preload.f, sol_eig_preload_tot.f, tol_f * max(sol_eig_preload_tot.f));
%! tol_KTAU0 = 1e-5;
%! assert_simple(max(max(abs(mat_ass.KTAU0 - mat_ass_tot.KTAU0))) / max(max(abs(mat_ass_tot.KTAU0))) < tol_KTAU0);
%! tol_tau = 1e-5;
%! for i=1:size(sol_stat.stress.tau.iso8, 1)
%!   for j=1:size(sol_stat.stress.tau.iso8, 2)
%!     for k=1:size(sol_stat.stress.tau.iso8, 3)
%!       assert_simple(sum(sol_stat.stress.tau.iso8(i, j, k, :)(:) .* omegaq), sol_stat_tot.stress.tau.iso8(i, j, k), tol_tau * max(max(max(abs(sol_stat_tot.stress.tau.iso8)))));
%!     endfor
%!   endfor
%! endfor
%! omegax_pre = sqrt(rho * A * L^4 / (E * I)) * 2 * pi * sol_eig_preload_tot.f;
%! assert_simple(omegax_pre(1), interp1(mu_ref, omegax_ref_chordwise, mu, "extrap"), 0.3);
%! assert_simple(omegax_pre(2), interp1(mu_ref, omegax_ref_flapwise, mu, "extrap"), 0.3);
%! omegax_rot = sqrt(rho * A * L^4 / (E * I)) * 2 * pi * sol_eig_rot.f(find(sol_eig_rot.f >= 0));
%! assert_simple(omegax_rot(1), interp1(mu_ref, omegax_ref_chordwise, mu, "extrap"), 0.3);
%! assert_simple(omegax_rot(2), interp1(mu_ref, omegax_ref_flapwise, mu, "extrap"), 0.3);
