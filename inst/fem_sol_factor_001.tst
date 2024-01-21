## fem_sol_factor.m:01
%!test
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"pastix", "pardiso", "mumps", "umfpack", "lu", "chol", "mldivide"};
%!   M = 20;
%!   options.refine_max_iter = 10;
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
%!   tol = sqrt(eps);
%!   for s=[false, true]
%!     options.pre_scaling = s;
%!     for N=2.^(1:8)
%!       A = gallery("poisson", N);
%!       A += A.';
%!       b = rand(columns(A), M);
%!       for i=1:numel(solvers)
%!         options.solver = solvers{i};
%!         Afact = fem_sol_factor(A, options);
%!         x = Afact \ b;
%!         f = max(norm(A * x - b, "cols") ./ norm(A * x + b, "cols"));
%!         assert_simple(f < tol);
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
