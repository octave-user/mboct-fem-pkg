## fem_tests.m:78
%!test
%! try
%! assert_simple(isinteger(FEM_SCA_TOT_MASS));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
