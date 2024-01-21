## fem_cms_constr_mat.m:02
%!test
%! C = [zeros(3, 6), eye(3), -eye(3)];
%! [T, res] = fem_cms_constr_mat(C);
%! assert_simple(size(T),[12, 9]);
%! assert_simple(res <= eps^0.9);
