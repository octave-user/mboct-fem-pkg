## fem_pre_mesh_struct_create.m:12
%!test
%! try
%! ## TEST12
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.N = 200;
%! param.fmin = 0;
%! param.fmax = 15000 / (SI_unit_second^-1);
%! param.num_freq = 1000;
%! damp.D = [1e-2; 0.03e-2];
%! damp.f = [20, 15000] / (SI_unit_second^-1);
%! helspr1.L = 25.8e-3 / SI_unit_meter;
%! helspr1.Di = 12.12e-3 / SI_unit_meter;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.n = 5.3;
%! helspr1.ni = 3;
%! helspr1.ng = 0.75;
%! helspr1.m = 400;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! [helspr1.material.alpha, helspr1.material.beta] = fem_pre_mat_rayleigh_damping(damp.D, damp.f);
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
%! mesh.material_data(1).beta = helspr1.material.beta;
%! mesh.material_data(1).alpha = helspr1.material.alpha;
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
%! load_case_dyn = fem_pre_load_case_create_empty(columns(Udyn));
%! for i=1:numel(load_case_dyn)
%!   load_case_dyn(i).joints = struct("U", repmat({zeros(6, 1)}, numel(mesh.elements.joints), 1));
%!   load_case_dyn(i).joints(2).U = Udyn(:, i);
%! endfor
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_DAMPING, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case_dyn);
%! opt_linsol.solver = "umfpack";
%! opt_linsol.pre_scaling = true;
%! opt_linsol.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_linsol.refine_max_iter = int32(50);
%! opt_linsol.epsilon_refinement = 1e-10;
%! opt_linsol.verbose = int32(0);
%! opt_eig.p = 5 * param.N;
%! opt_eig.maxit = int32(100);
%! omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%! [sol_eig, Phi, ~] = fem_sol_modal(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eig);
%! opt_post.elem_types = {"beam2"};
%! Freact = complex(zeros(columns(dof_map.edof.joints), columns(mat_ass.R), rows(dof_map.edof.joints), numel(omega)));
%! opt_factor.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_factor.verbose = int32(0);
%! for i=1:numel(omega)
%!   fprintf(stderr, "%3d:%.2fHz\n", i, omega(i) / (2 * pi));
%!   A = -omega(i)^2 * mat_ass.M + 1j * omega(i) * mat_ass.D + mat_ass.K;
%!   Uij = fem_sol_factor(A, opt_linsol) \ mat_ass.R;
%!   for k=1:size(Freact, 3)
%!     for j=1:size(Freact, 2)
%!       Freact(:, j, k, i) = Uij(dof_map.edof.joints(k, :), j) * mat_ass.mat_info.beta(2);
%!     endfor
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
