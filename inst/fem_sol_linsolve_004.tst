## fem_sol_linsolve.m:04
%!test
%! try
%! ## TEST 4
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     K = gallery("Poisson", N);
%!     R = complex(rand(columns(K), 10), rand(columns(K), 10));
%!     solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!     tol = sqrt(eps);
%!     for i=1:numel(solvers)
%!       options.solver = solvers{i};
%!       options.refine_max_iter = int32(10);
%!       options.pre_scaling = s;
%!       U = fem_sol_linsolve(K, R, options);
%!       f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%!       assert_simple(f < tol);
%!     endfor
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
