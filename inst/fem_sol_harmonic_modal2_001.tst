## fem_sol_harmonic_modal2.m:01
%!test
%! try
%! m1 = 0.5;
%! m2 = 30;
%! m3 = 0;
%! m4 = 0;
%! m5 = 0;
%! m6 = 0;
%! k1 = 1e-3;
%! k2 = 10000;
%! k3 = 1e10;
%! k4 = 1e-10;
%! d1 = 50;
%! d2 = 1e2;
%! d3 = 0;
%! d4 = 0;
%! r1 = 1e3;
%! r2 = -0.5;
%! r3 = 1e5;
%! r4 = -1e4;
%! r5 = 1e3;
%! r6 = 1e2;
%! c5 = c6 = norm([k1, k2, k3, k4]);
%! omega = linspace(0, 2 * pi * 10, 100000);
%! num_modes = 2;
%! mat_ass.M = [m1,  0,  0,  0,  0,  0;
%!              0,  m2,  0,  0,  0,  0;
%!              0,   0, m3,  0,  0,  0;
%!              0,   0,  0, m4,  0,  0;
%!              0,   0,  0,  0, m5,  0;
%!              0,   0,  0,  0,  0, m6];
%! mat_ass.K = [k1, 0,  0,  0,  0,  0;
%!              0, k2,  0,  0,  0,  0;
%!              0,  0, k3,  0, c5,  0;
%!              0,  0,  0, k4,  0, c6;
%!              0,  0, c5,  0,  0,  0;
%!              0,  0,  0, c6,  0,  0];
%! mat_ass.D = [d1,   0,   0,  0, 0, 0;
%!              0,  d2,   0,  0, 0, 0;
%!              0,   0,  d3,  0, 0, 0;
%!              0,   0,   0, d4, 0, 0;
%!              0,   0,   0,  0, 0, 0;
%!              0,   0,   0,  0, 0, 0];
%! mat_ass.R = [r1;
%!              r2;
%!              r3;
%!              r4;
%!              r5;
%!              r6];
%! mat_ass.R = [mat_ass.R, diag(mat_ass.R)];
%! mat_ass.M = sparse(mat_ass.M);
%! mat_ass.D = sparse(mat_ass.D);
%! mat_ass.K = sparse(mat_ass.K);
%! mat_ass.R = sparse(mat_ass.R);
%! mesh.nodes = zeros(2, 6);
%! dof_map.domain = FEM_DO_STRUCTURAL;
%! dof_map.totdof = rows(mat_ass.M);
%! dof_map.ndof = zeros(4, 6, "int32");
%! dof_map.ndof(1, 1:2) = [1:2];
%! dof_map.ndof(2, 1:2) = [3:4];
%! solvers = {"pastix", "pardiso", "mumps", "chol", "umfpack", "lu", "mldivide"};
%! for i=1:numel(solvers)
%!   fprintf(stderr, "linear solver: \"%s\"\n", solvers{i});
%!   if (~fem_sol_check_func(solvers{i}))
%!     fprintf(stderr, "solver %s is not available\n", solvers{i});
%!     continue;
%!   endif
%!   opt_sol.solver = solvers{i};
%!   opt_sol.verbose = int32(0);
%!   opt_sol.refine_max_iter = int32(100);
%!   opt_sol.pre_scaling = true;
%!   opt_sol.scaling = true;
%!   opt_eig.disp = int32(0);
%!   opt_eig.maxit = int32(10);
%!   switch (solvers{i})
%!   case "pardiso"
%!     opt_sol.symmetric = false; ## FIXME: Pardiso is not able to solve it in symmetric mode with weighted matching enabled
%!     opt_sol.weighted_matching = false;
%!   otherwise
%!     opt_sol.symmetric = true;
%!     opt_sol.weighted_matching = true;
%!   endswitch
%!   opt_sol.algorithm = "generic";
%!   [sol, Phi] = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, opt_sol, opt_eig);
%!   [dgen, kgen, rgen] = fem_sol_modes_scale2(mat_ass.M, mat_ass.D, mat_ass.K, Phi, mat_ass.R);
%!   U = fem_sol_harmonic_modal2(dgen, kgen, rgen, Phi(dof_map.ndof(1, 1:2), :), omega);
%!   Uref = [r1 ./ (-omega.^2 * m1 + 1j * omega * d1 + k1);
%!           r2 ./ (-omega.^2 * m2 + 1j * omega * d2 + k2)];
%!   omega01 = sqrt(k1 / m1);
%!   omega02 = sqrt(k2 / m2);
%!   delta1 = d1 / (2 * m1);
%!   delta2 = d2 / (2 * m2);
%!   D1 = delta1 / omega01;
%!   D2 = delta2 / omega02;
%!   Dref = [D1, D1, D2, D2];
%!   lambdaref = 1j * [omega01, omega02];
%!   [~, idx] = sortrows([imag(lambdaref)(:), real(lambdaref)(:)]);
%!   lambdaref = lambdaref(idx);
%!   Dref = min(1, Dref(idx));
%!   tol = eps^0.9;
%!   assert_simple(U(:, :, 1), Uref, tol * norm(Uref));
%!   assert_simple(sum(U(:, :, 2:end), 3), Uref, tol * norm(Uref));
%!   assert_simple(sum(mat_ass.R(:, 2:end), 2), mat_ass.R(:, 1), tol * norm(mat_ass.R(:,1)));
%!   assert_simple(imag(sol.lambda), imag(lambdaref), tol * norm(imag(lambdaref)));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
