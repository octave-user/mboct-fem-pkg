## fem_tests.m:96
%!test
%! try
%! assert_simple(isinteger(FEM_VEC_PARTICLE_VELOCITY));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
