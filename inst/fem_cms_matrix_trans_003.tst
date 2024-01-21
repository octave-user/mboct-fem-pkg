## fem_cms_matrix_trans.m:03
%!test
%! N = 1000;
%! for j=0:1
%! for i=1:100
%! A = sprand(N, N, 0.01) + abs(diag(rand(N, 1)));
%! A *= A.';
%! [r, c, d] = find(A);
%! if j
%! idx = find(r <= c);
%! mat_type = "Upper";
%! else
%! idx = find(r >= c);
%! mat_type = "Lower";
%! endif
%! r = r(idx);
%! c = c(idx);
%! d = d(idx);
%! Asym = sparse(r, c, d, N, N);
%! switch matrix_type(Asym)
%! case "Upper"
%! assert_simple(j, 1);
%! case "Lower"
%! assert_simple(j, 0);
%! otherwise
%! assert_simple(false);
%! endswitch
%! T = eye(N, N)(:, 1:floor(0.1 * N));
%! TAT = fem_cms_matrix_trans(T, Asym, mat_type);
%! assert_simple(TAT, T.' * A * T, eps);
%! assert_simple(issymmetric(TAT));
%! endfor
%! endfor
