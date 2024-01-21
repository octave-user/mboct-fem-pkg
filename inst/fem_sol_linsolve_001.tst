## fem_sol_linsolve.m:01
%!test
%! ### TEST 1
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!   for s=[false,true]
%!     for N=[2,10,20,50,100]
%!       for i=1:10
%!         A = gallery("Poisson", N);
%!         b = 2 * rand(rows(A), 3) - 1;
%!         for i=1:numel(solvers)
%!           options.solver = solvers{i};
%!           options.refine_max_iter = int32(10);
%!           options.number_of_threads = int32(1);
%!           options.pre_scaling = s;
%!           options.epsilon_refinement = eps^0.9;
%!           options.verbose = int32(0);
%!           x = fem_sol_linsolve(A, b, options);
%!           f = max(norm(A * x - b, "cols") ./ norm(A * x + b, "cols"));
%!           assert_simple(f <= options.epsilon_refinement);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
