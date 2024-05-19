## fem_sol_eigs.m:01
%!test
%! try
%! tol = 1e-4;
%! N = 300;
%! n = 10;
%! toleigs = 0;
%! solvers = {"pastix", "pardiso", "mumps", "umfpack", "lu", "chol", "mldivide"};
%! alg={"symmetric-inverse","shift-invert","unsymmetric"};
%! problems = {"structural", "thermal", "acoustic", "pressure", "buckling"};
%! t = zeros(numel(alg), numel(solvers));
%! for k=1:numel(solvers)
%!   for a=1:numel(alg)
%!     for j=1:2
%!       rand("seed", 0);
%!       for i=1:3
%!         K = sprand(N, N, 0.01) + 10*abs(diag(0.1+rand(N, 1)));
%!         M = sprand(rows(K), columns(K), 0.01) + 10*abs(diag(0.1+rand(rows(K), 1)));
%!         K *= K.';
%!         M *= M.';
%!         Q = symrcm(M + K);
%!         M = M(Q, Q);
%!         K = K(Q, Q);
%!         if j == 2
%!           rho = -0.1 * max(abs(diag(K))) / max(abs(diag(M)));
%!         else
%!           rho = 0;
%!         endif
%!         start = tic();
%!         for l=1:7
%!           switch (l)
%!             case 1
%!               [U, lambda, err] = fem_sol_eigs(K, M, n, rho, toleigs, alg{a}, solvers{k});
%!             otherwise
%!               opt = struct();
%!               opt.disp = 0;
%!               opt.maxit = 100;
%!               if (l > 2)
%!                 opt.rho = rho;
%!               endif
%!               if (l > 3)
%!                 opt.tolerance = toleigs;
%!               endif
%!               if (l > 4)
%!                 opt.algorithm = alg{a};
%!               endif
%!               if (l > 5)
%!                 opt.solver = solvers{k};
%!               endif
%!               if (l > 6)
%!                 opt.number_of_threads = mbdyn_solver_num_threads_default();
%!               endif
%!               if (l > 7)
%!                 for m=1:numel(problems)
%!                   opt.problem = opt.problem = problems{m};
%!                   [U, lambda, err] = fem_sol_eigs(K, M, n, opt);
%!                   assert_simple(max(err) < tol);
%!                 endfor
%!               else
%!                   [U, lambda, err] = fem_sol_eigs(K, M, n, opt);
%!                   assert_simple(max(err) < tol);
%!               endif
%!           endswitch
%!           t(a, k) += toc(start);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! for i=1:numel(alg)
%!   for j=1:numel(solvers)
%!     fprintf(stderr, "algorithm \"%s:%s\": %.2fs\n", alg{i}, solvers{j}, t(i, j));
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
