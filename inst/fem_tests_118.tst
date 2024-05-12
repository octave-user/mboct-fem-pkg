## fem_tests.m:118
%!test
%! try
%! assert_simple(isinteger(FEM_DO_ACOUSTICS));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
