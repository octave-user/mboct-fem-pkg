## fem_tests.m:61
%!test ## test 61
%! try
%! assert_simple(isinteger(FEM_MAT_DAMPING));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
