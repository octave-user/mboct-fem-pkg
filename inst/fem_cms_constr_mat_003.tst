## fem_cms_constr_mat.m:03
%!test
%! C = [eye(3), zeros(3, 6), -eye(3)];
%! [T, res] = fem_cms_constr_mat(C);
%! assert_simple(size(T),[12, 9]);
%! assert_simple(res <= eps^0.9);
