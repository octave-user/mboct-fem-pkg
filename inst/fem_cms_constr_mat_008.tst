## fem_cms_constr_mat.m:08
%!test
%! try
%! C = [eye(3), -eye(3), zeros(3, 6)];
%! [T, res] = fem_cms_constr_mat(C);
%! f = C * T;
%! assert_simple(norm(f) < eps);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
