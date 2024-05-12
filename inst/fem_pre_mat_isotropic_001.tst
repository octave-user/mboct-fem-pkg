## fem_pre_mat_isotropic.m:01
%!test
%! try
%! E = 210000;
%! nu = 0.3;
%! C = fem_pre_mat_isotropic(E, nu);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
