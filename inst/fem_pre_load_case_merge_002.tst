## fem_pre_load_case_merge.m:02
%!test
%! try
%! load_case = fem_pre_load_case_merge();
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
