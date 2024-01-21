## fem_sol_linsolve.m:03
%!test
%! ## TEST 3
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     K = gallery("Poisson", N);
%!     K = complex(K, rand() * K);
%!     K += K';
%!     R = complex(rand(columns(K), 10), rand(columns(K), 10));
%!     solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!     tol = sqrt(eps);
%!     for i=1:numel(solvers)
%!       options.solver = solvers{i};
%!       options.refine_max_iter = int32(30);
%!       options.verbose = int32(0);
%!       options.pre_scaling = s;
%!       options.symmetric = false;
%!       U = fem_sol_linsolve(K, R, options);
%!       f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%!       assert_simple(f < tol);
%!     endfor
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
