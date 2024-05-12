## fem_cms_matrix_trans.m:02
%!test
%! try
%! for t={"Upper", "Lower"}
%! N = 100;
%! for i=1:100
%! A = sprand(N, N, 0.1);
%! A += A.';
%! T = eye(N, N)(:, 1:floor(0.9*N));
%! TAT = fem_cms_matrix_trans(T, A, t);
%! assert_simple(TAT, T.' * A * T, eps);
%! assert_simple(issymmetric(TAT));
%! endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
