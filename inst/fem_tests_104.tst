## fem_tests.m:104
%!test
%! try
%! assert_simple(isinteger(FEM_MAT_MASS_FLUID_STRUCT_IM));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
