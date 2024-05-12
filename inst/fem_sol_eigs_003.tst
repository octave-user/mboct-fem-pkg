## fem_sol_eigs.m:03
%!test
%! try
%! rand("seed", 0);
%! tol = eps^0.3;
%! N = 10;
%! M = 3;
%! for j=1:2
%! for i=1:10
%! A = rand(N, N) + abs(diag(rand(N, 1)));
%! if j == 1
%! B = eye(columns(A));
%! else
%! B = rand(rows(A), columns(A)) + abs(diag(rand(rows(A), 1)));
%! endif
%! A *= A.';
%! B *= B.';
%! opts.isreal = 1;
%! opts.issym = 1;
%! opts.tol = realmin();
%! opts.maxiter = 5000;
%! [v, lambda] = eig(A, B);
%! [v3, lambda3] = eig(B \ A);
%! [lambda, idx_lambda] = sort(diag(lambda));
%! v = v(:, idx_lambda);
%! [lambda3, idx_lambda3] = sort(diag(lambda3));
%! v3 = v3(:, idx_lambda3);
%! [L, P] = chol(B, "lower");
%! assert_simple(P,0);
%! A2 = (L \ A) / L.';
%! assert_simple(issymmetric(A2, tol*norm(A2)));
%! [L2, P2] = chol(A2, "lower");
%! assert_simple(P2,0);
%! [v2, lambda2, flag2] = eigs(@(y) L2.' \ (L2 \ y), columns(A2), M, "sm", opts);
%! [lambda2, idx_lambda2] = sort(diag(lambda2));
%! v2 = v2(:, idx_lambda2);
%! v2 = L.' \ v2;
%! assert_simple(flag2, 0);
%! assert_simple(lambda(1:3), lambda3(1:M), tol*max(max(abs(lambda(1:M)))));
%! assert_simple(lambda(1:M), lambda2, tol*max(max(abs(lambda(1:M)))));
%! for k=1:M
%!  assert_simple(A * v(:, k), lambda(k) * B * v(:, k), tol * norm(A * v(:, k)));
%!  assert_simple(A * v3(:, k), lambda3(k) * B * v3(:, k), tol * norm(A * v3(:, k)));
%!  assert_simple(A * v2(:, k), lambda2(k) * B * v2(:, k), tol * norm(A * v2(:, k)));
%! endfor
%! endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
