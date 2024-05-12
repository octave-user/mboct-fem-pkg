## fem_tests.m:77
%!test
%! try
%! assert_simple(isinteger(FEM_SCA_STRESS_VMIS));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
