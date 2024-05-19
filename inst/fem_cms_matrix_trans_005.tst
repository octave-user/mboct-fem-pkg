## fem_cms_matrix_trans.m:05
%!test
%! try
%! A = rand(10, 10);
%! A *= A.';
%! T = rand(10, 3);
%! TAT1 = fem_cms_matrix_trans(T, A, "Upper");
%! TAT2 = fem_cms_matrix_trans(T, A, "Lower");
%! assert_simple(isdefinite(A));
%! assert_simple(isdefinite(TAT1));
%! assert_simple(TAT2, TAT1, eps * norm(TAT1));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
