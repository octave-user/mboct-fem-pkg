## fem_pre_mesh_struct_create.m:13
%!test
%! try
%! ## TEST13
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.N1 = 20;
%! param.N2 = 200;
%! param.fmin = 0;
%! param.fmax = 1000 / (SI_unit_second^-1);
%! param.num_freq = 1000;
%! helspr1.L = 25.8e-3 / SI_unit_meter;
%! helspr1.Di = 12.12e-3 / SI_unit_meter;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.n = 5.3;
%! helspr1.ni = 3;
%! helspr1.ng = 0.75;
%! helspr1.m = 40;
%! helspr1.material.D = 0.001578;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! helspr1.material.tan_delta = 2 * helspr1.material.D;
%! helspr1.material.nu = helspr1.material.E / (2 * helspr1.material.G) - 1;
%! helspr1.material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! helspr1.D = helspr1.Di + helspr1.d;
%! section1.A = helspr1.d^2 * pi / 4;
%! section1.Ay = 9 / 10 * section1.A;
%! section1.Az = section1.Ay;
%! section1.Iy = helspr1.d^4 * pi / 64;
%! section1.Iz = section1.Iy;
%! section1.It = section1.Iy + section1.Iz;
%! helspr1.Phi = linspace(0, 2 * pi * helspr1.n, ceil(helspr1.m * helspr1.n)).' + 2 * pi * helspr1.ni;
%! helspr1.x = 0.5 * helspr1.D * cos(helspr1.Phi);
%! helspr1.y = 0.5 * helspr1.D * sin(helspr1.Phi);
%! helspr1.z = (helspr1.L - helspr1.d * (2 * (helspr1.ni - helspr1.ng) + 1)) * linspace(0, 1, numel(helspr1.Phi))(:) + helspr1.d * (helspr1.ni - helspr1.ng + 0.5);
%! helspr1.e2 = [0, 0, 1];
%! elnodes = int32([(1:numel(helspr1.Phi) - 1).', (2:numel(helspr1.Phi)).']);
%! mesh.nodes = [[helspr1.x, helspr1.y, helspr1.z], zeros(numel(helspr1.Phi), 3)];
%! mesh.elements.beam2 = struct("nodes", mat2cell(elnodes, ones(numel(helspr1.Phi) - 1, 1, "int32"), 2), ...
%!                              "section", mat2cell(repmat(section1, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32")), ...
%!                              "e2", mat2cell(repmat(helspr1.e2, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32"), 3));
%! mesh.material_data(1).E = helspr1.material.E;
%! mesh.material_data(1).nu = helspr1.material.nu;
%! mesh.material_data(1).tan_delta = helspr1.material.tan_delta;
%! mesh.material_data(1).rho = helspr1.material.rho;
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.beam2 = repmat(int32(1), numel(helspr1.Phi) - 1, 1);
%! empty_cell = cell(1, 2);
%! mesh.elements.joints = struct("nodes", empty_cell, "C", empty_cell);
%! mesh.elements.joints(1).nodes = int32(1);
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.joints(2).nodes = rows(mesh.nodes);
%! mesh.elements.joints(2).C = eye(6);
%! Udyn = eye(6);
%! omega = linspace(0, 2 * pi * 15000, 10000);
%! g = eye(3);
%! load_case_dyn = fem_pre_load_case_create_empty(columns(Udyn) + columns(g));
%! for i=1:numel(load_case_dyn)
%!   load_case_dyn(i).joints = struct("U", repmat({zeros(6, 1)}, numel(mesh.elements.joints), 1));
%!   load_case_dyn(i).g = zeros(3, 1);
%! endfor
%! for i=1:columns(Udyn)
%!   load_case_dyn(i).joints(2).U = Udyn(:, i);
%! endfor
%! for i=1:columns(g)
%!   load_case_dyn(i + columns(Udyn)).g = g(:, i);
%! endfor
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.K_im, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_STIFFNESS_IM, ...
%!                                       FEM_MAT_DAMPING, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case_dyn);
%! opt_linsol.solver = "umfpack";
%! opt_linsol.pre_scaling = true;
%! opt_linsol.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_linsol.refine_max_iter = int32(50);
%! opt_linsol.epsilon_refinement = 1e-10;
%! opt_linsol.verbose = int32(0);
%! opt_eig.p = 5 * param.N2;
%! opt_eig.maxit = int32(100);
%! [sol_eig, PhiD] = fem_sol_modal(mesh, dof_map, mat_ass, param.N2, opt_linsol, opt_eig);
%! PhiD *= diag(1 ./ norm(PhiD, "cols"));
%! [sol_stat, PhiS] = fem_sol_static(mesh, dof_map, mat_ass, opt_linsol);
%! PhiS *= diag(1 ./ norm(PhiS, "cols"));
%! PhiR = PhiS(:, columns(Udyn) + (1:columns(g)));
%! PhiS = PhiS(:, 1:columns(Udyn));
%! idx_retained = false(rows(dof_map.ndof), 1);
%! idx_retained(mesh.elements.joints(2).nodes) = true;
%! idx_retained = dof_map.ndof(find(idx_retained), :).'(:);
%! idx_retained = idx_retained(idx_retained > 0);
%! idx_dropped = true(rows(dof_map.ndof), 1);
%! idx_dropped(idx_retained) = false;
%! idx_dropped = dof_map.ndof(find(idx_dropped), :).'(:);
%! idx_dropped = idx_dropped(idx_dropped > 0);
%! phi = PhiD(idx_dropped, :);
%! psi = PhiS(idx_retained, :);
%! Xir = PhiR(idx_dropped, :);
%! M = mat_ass.M(idx_dropped, idx_dropped);
%! m = fem_cms_matrix_trans(phi, M, "Upper");
%! Xir -= phi * (m \ (phi.' * M * Xir));
%! Xir *= diag(1 ./ norm(Xir, "cols"));
%! Mrr = fem_cms_matrix_trans(Xir, mat_ass.M(idx_dropped, idx_dropped), "Upper");
%! Krr = fem_cms_matrix_trans(Xir, mat_ass.K(idx_dropped, idx_dropped), "Upper");
%! [Vrr, lambdarr] = eig(Krr, Mrr);
%! phir = Xir * Vrr;
%! phir *= diag(1 ./ norm(phir, "cols"));
%! PhiR = zeros(dof_map.totdof, columns(Xir));
%! PhiR(idx_dropped, :) = phir;
%! Phi = [PhiD, PhiS, PhiR];
%! phiii = zeros(columns(Phi), columns(Phi));
%! for i=1:columns(Phi)
%!   for j=1:columns(Phi)
%!     phiii(i, j) = Phi(idx_dropped, i).' * Phi(idx_dropped, j);
%!   endfor
%! endfor
%! mat_ass.K += 1j * mat_ass.K_im;
%! Phi *= diag(1 ./ norm(Phi, "cols"));
%! Mred = fem_cms_matrix_trans(Phi, mat_ass.M, "Upper");
%! Dred = fem_cms_matrix_trans(Phi, mat_ass.D, "Upper");
%! Kred = fem_cms_matrix_trans(Phi, mat_ass.K, "Upper");
%! Rred = Phi.' * mat_ass.R(:, 1:numel(load_case_dyn));
%! omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%! opt_eig.p = 5 * param.N1;
%! opt_eig.maxit = int32(100);
%! [sol_eigd] = fem_sol_modal_damped(mesh, dof_map, mat_ass, param.N1, opt_linsol, opt_eig);
%! opt_post.elem_types = {"beam2"};
%! Freact = Freact_red = complex(zeros(columns(dof_map.edof.joints), columns(mat_ass.R), rows(dof_map.edof.joints), numel(omega)));
%! opt_factor.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_factor.verbose = int32(0);
%! for i=1:numel(omega)
%!   fprintf(stderr, "%3d:%.2fHz\n", i, omega(i) / (2 * pi));
%!   A = -omega(i)^2 * mat_ass.M + 1j * omega(i) * mat_ass.D + mat_ass.K;
%!   Uij = fem_sol_factor(A, opt_linsol) \ mat_ass.R(:, 1:numel(load_case_dyn));
%!   Uij_red = Phi * ((-omega(i)^2 * Mred + 1j * omega(i) * Dred + Kred) \ Rred);
%!   for k=1:size(Freact, 3)
%!     for j=1:size(Freact, 2)
%!       Freact(:, j, k, i) = Uij(dof_map.edof.joints(k, :), j) * mat_ass.mat_info.beta(2);
%!       Freact_red(:, j, k, i) = Uij_red(dof_map.edof.joints(k, :), j) * mat_ass.mat_info.beta(2);
%!     endfor
%!   endfor
%! endfor
%! tol = 1e-4;
%! tol2 = 1e-2;
%! assert_simple(max(abs(sol_eigd.D - helspr1.material.D)) < tol * helspr1.material.D);
%! assert_simple(max(abs(Freact_red(:) - Freact(:))) < tol2 * max(abs(Freact(:))));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
