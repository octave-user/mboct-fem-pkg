## fem_cms_matrix_trans.m:01
%!test
%! N = 10;
%! for t={"Upper", "Lower"}
%! for i=1:100
%! A = rand(N, N);
%! A += A.';
%! T = eye(N, N)(:, 1:floor(0.75*N));
%! TAT = fem_cms_matrix_trans(T, A, t);
%! assert_simple(TAT, T.' * A * T, eps);
%! assert_simple(issymmetric(TAT));
%! endfor
%! endfor
