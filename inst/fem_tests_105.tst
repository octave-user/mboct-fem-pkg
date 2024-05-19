## fem_tests.m:105
%!test
%! try
%! assert_simple(isinteger(FEM_MAT_DAMPING_FLUID_STRUCT_RE));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
