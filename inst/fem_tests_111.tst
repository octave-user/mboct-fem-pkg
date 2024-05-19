## fem_tests.m:111
%!test
%! try
%! assert_simple(isinteger(FEM_VEC_COLL_THERMAL_COND));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
