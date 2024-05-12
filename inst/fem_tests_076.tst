## fem_tests.m:76
%!test
%! try
%! assert_simple(isinteger(FEM_MAT_STIFFNESS_SYM_L));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
