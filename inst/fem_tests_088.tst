## fem_tests.m:88
%!test
%! try
%! assert_simple(isinteger(FEM_VEC_LOAD_THERMAL));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
