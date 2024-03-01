## fem_sol_eigs.m:05
%!test
%! state = rand("state");
%! opts.refine_max_iter = int32(3);
%! opts.epsilon_refinement = eps^0.5;
%! opts.pre_scaling = true;
%! opts.rho = -1;
%! solvers={"lu", "chol", "mldivide"};
%! unwind_protect
%!   rand("seed", 0);
%!   for isol=1:numel(solvers)
%!     opts.solver = solvers{isol};
%!     for N=[2, 10, 50, 100]
%!       num_modes = ceil(0.1 * N);
%!       for i=1:5
%!         K = rand(N, N);
%!         K *= K.';
%!         M = rand(N, N);
%!         M *= M.';
%!         fprintf(stderr, "%s: %d: %d\n", opts.solver, N, i);
%!         [U, lambda] = fem_sol_eigs(K, M, num_modes, opts);
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
