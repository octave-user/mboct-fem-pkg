## fem_cms_matrix_trans.m:04
%!test
%! for i=1:100
%! A = [rand(10, 10), rand(10, 5);
%!      rand(5, 10), zeros(5, 5)];
%! A += A.';
%! T = rand(rows(A), floor(0.75 * columns(A)));
%! TAT = fem_cms_matrix_trans(T, A, "Lower");
%! assert_simple(TAT, T.' * A * T, eps^0.7);
%! endfor
