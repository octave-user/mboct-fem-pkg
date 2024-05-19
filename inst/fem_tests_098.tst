## fem_tests.m:98
%!test
%! try
%! assert_simple(isinteger(FEM_VEC_PARTICLE_VELOCITY_C));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
