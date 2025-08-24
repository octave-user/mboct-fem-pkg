%!test
%! try
%! function varargout = callback(stage, opt_sol, idx, omega, Ui, status, Mi, Di, Ki, Ri)
%!   switch (stage)
%!   case "init"
%!     U = complex(zeros(columns(Mi), numel(omega)));
%!     varargout = {U, opt_sol};
%!   case "pre"
%!     varargout = {Mi, Di, Ki, Ri};
%!   case "post"
%!     varargout = {Ui};
%!     fprintf(stderr, "%d: %d:%d: %e\n", idx, status.gmres.flag, status.do_fact, status.gmres.relres);
%!   endswitch
%! endfunction
%!
%! k = 100;
%! d = 0.1;
%! m = 10;
%! K = [k,    -k,  0, 1, 0;
%!     -k, 2 * k, -k, 0, 0;
%!      0,    -k,  k, 0, 1;
%!      1,     0,  0, 0, 0;
%!      0,     0,  1, 0, 0];
%!
%! M = [2 / 3 * m, 1 / 3 * m,         0, 0, 0;
%!      1 / 3 * m, 4 / 3 * m, 1 / 3 * m, 0, 0;
%!              0, 1 / 3 * m, 2 / 3 * m, 0, 0;
%!              0,         0,         0, 0, 0;
%!              0,         0,         0, 0, 0];
%!
%! D = [d,    -d,  0, 0, 0;
%!     -d, 2 * d, -d, 0, 0;
%!      0,    -d,  d, 0, 0;
%!      0,     0,  0, 0, 0;
%!      0,     0,  0, 0, 0];
%!
%! R = [0; 0; 0; 1; 0];
%!
%! omega = linspace(0, 2 * pi * 2, 1000);
%! opt_sol.scale.maxiter = int32(100);
%! opt_sol.scale.tol = sqrt(eps);
%! opt_sol.factor.solver = "lu";
%! opt_sol.gmres.restart = columns(M);
%! opt_sol.gmres.tol = eps^0.9;
%! opt_sol.gmres.maxiter = 3;
%! opt_sol.gmres.alpha = 0.5;
%! cb = @(stage, opt, idx, omega, Ui, status) callback(stage, opt_sol, idx, omega, Ui, status, M, D, K, R);
%! U1 = fem_sol_harmonic_gmres(omega, cb);
%! U2 = fem_sol_harmonic_gmres(omega, M, D, K, R, 1:columns(K), opt_sol);
%! Uref = complex(zeros(columns(M)), numel(omega));
%!
%! for i=1:columns(omega)
%!   Uref(:, i) = (-omega(i)^2 * M + 1j * omega(i) * D + K) \ R;
%! endfor
%!
%! tol = eps^0.7;
%!
%! assert_simple(Uref, U1, tol * max(norm(Uref, "cols")));
%! assert_simple(Uref, U2, tol * max(norm(Uref, "cols")));
%!
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%! function varargout = callback(stage, opt_sol, idx, omega, Ui, status, mat_ass, dof_map)
%!   switch (stage)
%!   case "init"
%!    Freact = complex(zeros(columns(dof_map.edof.joints) * columns(mat_ass.R) * rows(dof_map.edof.joints), numel(omega)));
%!    opt_sol.scale.maxiter = int32(100);
%!    opt_sol.scale.tol = sqrt(eps);
%!    opt_sol.factor.solver = "umfpack";
%!    opt_sol.factor.pre_scaling = true;
%!    opt_sol.factor.refine_max_iter = int32(0);
%!    opt_sol.gmres.restart = 20;
%!    opt_sol.gmres.tol = eps^0.8;
%!    opt_sol.factor.number_of_threads = mbdyn_solver_num_threads_default();
%!    opt_sol.factor.verbose = int32(0);
%!    opt_sol.gmres.alpha = 3;
%!    opt_sol.gmres.maxiter = 20;
%!    varargout = {Freact, opt_sol};
%!   case "pre"
%!    varargout = {mat_ass.M, mat_ass.D, mat_ass.K, mat_ass.R};
%!   case "post"
%!     if (mod(idx, 10) == 1)
%!       fprintf(stderr, "%3d:%.2fHz flag=%d, factor=%d, relres=%e iter=%d/%d\n", idx, omega / (2 * pi), status.gmres.flag, status.do_fact, status.gmres.relres, status.gmres.iter(1), status.gmres.iter(2));
%!     endif
%!     Freact = complex(zeros(columns(dof_map.edof.joints) * columns(mat_ass.R) * rows(dof_map.edof.joints), 1));
%!     for k=1:rows(dof_map.edof.joints)
%!       for j=1:columns(mat_ass.R)
%!         Freact((1:columns(dof_map.edof.joints)) + columns(dof_map.edof.joints) * ((j - 1) + columns(mat_ass.R)) * (k - 1)) = Ui(dof_map.edof.joints(k, :), j) * mat_ass.mat_info.beta(2);
%!       endfor
%!     endfor
%!     varargout = {Freact};
%!   endswitch
%! endfunction
%!
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.N = 200;
%! param.fmin = 0;
%! param.fmax = 15000 / (SI_unit_second^-1);
%! param.num_freq = 1000;
%! damp.D = [1e-2; 2e-2];
%! damp.f = [20, 15000] / (SI_unit_second^-1);
%! helspr1.L = 25.8e-3 / SI_unit_meter;
%! helspr1.Di = 12.12e-3 / SI_unit_meter;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.n = 5.3;
%! helspr1.ni = 3;
%! helspr1.ng = 0.75;
%! helspr1.m = 20;
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
%! omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%! opt_post.elem_types = {"beam2"};
%! cb = @(stage, opt_sol, idx, omega, Ui, status) callback(stage, opt_sol, idx, omega, Ui, status, mat_ass, dof_map);
%! start = tic();
%! Freact = fem_sol_harmonic_gmres(omega, cb);
%! time_direct = toc(start);
%! start = tic();
%! Freact_ref = complex(zeros(columns(dof_map.edof.joints) * columns(mat_ass.R) * rows(dof_map.edof.joints), numel(omega)));
%! opt_sol.solver = "umfpack";
%! opt_sol.refine_max_iter = int32(10);
%! for i=1:columns(omega)
%!   Ui = fem_sol_factor(-omega(i)^2 * mat_ass.M + 1j * omega(i) * mat_ass.D + mat_ass.K) \ mat_ass.R;
%!   if (mod(i, 10) == 1)
%!     fprintf(stderr, "%3d:%.2fHz\n", i, omega(i) / (2 * pi));
%!   endif
%!   for k=1:rows(dof_map.edof.joints)
%!     for j=1:columns(mat_ass.R)
%!       Freact_ref((1:columns(dof_map.edof.joints)) + columns(dof_map.edof.joints) * ((j - 1) + columns(mat_ass.R)) * (k - 1), i) = Ui(dof_map.edof.joints(k, :), j) * mat_ass.mat_info.beta(2);
%!     endfor
%!   endfor
%! endfor
%! tol = 1e-8;
%! assert_simple(max(abs(Freact(:) - Freact_ref(:))) < tol * norm(Freact_ref(:)));
%! time_full = toc(start);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
