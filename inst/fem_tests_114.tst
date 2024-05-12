## fem_tests.m:114
%!test
%! try
%! assert_simple(isinteger(FEM_VEC_COLL_MASS_FLUID_STRUCT));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
