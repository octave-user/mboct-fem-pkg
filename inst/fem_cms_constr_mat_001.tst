## fem_cms_constr_mat.m:01
%!test
%! try
%! C = [eye(3), -eye(3), zeros(3, 6)];
%! [T, res] = fem_cms_constr_mat(C);
%! assert_simple(size(T),[12, 9]);
%! assert_simple(res <= eps^0.9);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
