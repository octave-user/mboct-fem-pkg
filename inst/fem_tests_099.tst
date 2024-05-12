## fem_tests.m:99
%!test
%! try
%! assert_simple(isinteger(FEM_SCA_ACOUSTIC_INTENSITY_C));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
