## fem_pre_mat_rayleigh_damping.m:01
%!test
%! try
%! D = [1e-2;  0.1e-2];
%! f = [10; 1000];
%! omega = 2 * pi * f;
%! [alpha, beta] = fem_pre_mat_rayleigh_damping(D, f);
%! Dref = 1/2* (alpha ./ omega + beta * omega);
%! assert_simple(D, Dref, eps * norm(Dref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
