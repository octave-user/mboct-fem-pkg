## fem_cms_constr_mat.m:06
%!test
%! C = [1, 2, 3, 0, 0, 0;
%!      0, 0, 0, 1, 2, 3];
%! [T, res] = fem_cms_constr_mat(C);
%! assert_simple(size(T),[6, 4]);
%! assert_simple(res <= eps^0.9);
