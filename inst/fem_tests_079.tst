## fem_tests.m:79
%!test
%! try
%! assert_simple(isinteger(FEM_VEC_INERTIA_M1));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
