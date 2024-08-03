## fem_sol_eigs.m:05
%!test
%! try
%! state = rand("state");
%! opts.refine_max_iter = int32(100);
%! opts.epsilon_refinement = eps^0.7;
%! opts.pre_scaling = true;
%! gamma = 1e-2;
%! solvers={"lu", "chol", "mldivide"};
%! unwind_protect
%!   rand("seed", 0);
%!   idx = int32(0);
%!   for isol=1:numel(solvers)
%!     opts.solver = solvers{isol};
%!     for N=[2, 10, 50, 100]
%!       num_modes = ceil(0.1 * N);
%!       for i=1:5
%!         K = rand(N, N);
%!         K *= K.';
%!         M = rand(N, N);
%!         M *= M.';
%!         opts.rho = -gamma * norm(K) / norm(M);
%!         [U, lambda, err] = fem_sol_eigs(K, M, num_modes, opts);
%!         fprintf(stderr, "%s: %d: %d: %e\n", opts.solver, N, i, max(err));
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
