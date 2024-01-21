## fem_cms_constr_mat.m:05
%!test
%! C = [1, 2, 3, 0, 0, 0];
%! [T, res] = fem_cms_constr_mat(C);
%! assert_simple(size(T),[6, 5]);
%! assert_simple(res <= eps^0.9);
