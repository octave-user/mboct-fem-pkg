## fem_tests.m:66
%!test
%! try
%! assert_simple(isinteger(FEM_MAT_INERTIA_INV5));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
