## fem_sol_linsolve.m:07
%!test
%! try
%! N = 100;
%! K = gallery("Poisson", N);
%! R = linspace(0, 1, columns(K)).';
%! options.solver = "chol";
%! U = fem_sol_linsolve(K, R, options);
%! f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%! fprintf(stderr, "backward error: %.1e\n", f);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
