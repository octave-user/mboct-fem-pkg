## fem_sol_linsolve.m:05
%!test
%! ## TEST 5
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"pastix", "mumps", "pardiso", "umfpack", "chol", "lu"};
%!   tol = sqrt(eps);
%!   max_f = 0;
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     for i=1:10
%!       A = 2 * rand(N, N) - 1;
%!       A += 2j * rand(N, N) - 1;
%!       A += A';
%!       b = complex(2 * rand(N, 20) - 1, 2 * rand(N, 20) - 1);
%!       for i=1:numel(solvers)
%!         options.solver = solvers{i};
%!         options.refine_max_iter = int32(100);
%!         options.number_of_threads = int32(1);
%!         options.pre_scaling = s;
%!         options.symmetric = false;
%!         options.verbose = int32(0);
%!         x = fem_sol_linsolve(A, b, options);
%!         f = max(norm(A * x - b, "cols") ./ norm(A * x + b, "cols"));
%!         assert_simple(f < tol);
%!         max_f = max(f, max_f);
%!       endfor
%!     endfor
%!   endfor
%!   endfor
%!   assert_simple(max_f < tol);
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
