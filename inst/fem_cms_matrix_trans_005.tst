## fem_cms_matrix_trans.m:05
%!test
%! A = rand(10, 10);
%! A *= A.';
%! T = rand(10, 3);
%! TAT1 = fem_cms_matrix_trans(T, A, "Upper");
%! TAT2 = fem_cms_matrix_trans(T, A, "Lower");
%! assert_simple(isdefinite(A));
%! assert_simple(isdefinite(TAT1));
%! assert_simple(TAT2, TAT1, eps * norm(TAT1));
