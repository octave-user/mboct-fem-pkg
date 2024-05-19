## fem_sol_eigs.m:04
%!test
%! try
%! rand("seed", 0);
%! tol1 = 1e-4;
%! tol2 = 1e-7;
%! N = 100;
%! M = 3;
%! time_lambda = 0;
%! time_lambda2 = 0;
%! for i=1:50
%! A = sprand(N, N, 0.01) + 1000*abs(diag(rand(N, 1)));
%! B = sprand(rows(A), columns(A), 0.01) + 1000*abs(diag(rand(rows(A), 1)));
%! A *= A.';
%! B *= B.';
%! start = tic();
%! [v, lambda] = eig(A, B);
%! time_lambda += tic() - start;
%! lambda = diag(lambda);
%! [lambda, idx_lambda] = sort(lambda);
%! lambda = lambda(1:M);
%! idx_lambda = idx_lambda(1:M);
%! v = v(:, idx_lambda);
%! [v3, lambda3] = eig(B \ A);
%! lambda3 = diag(lambda3);
%! [lambda3, idx_lambda3] = sort(lambda3);
%! lambda3 = lambda3(1:M);
%! idx_lambda3 = idx_lambda3(1:M);
%! v3 = v3(:, idx_lambda3);
%! [L] = chol(A, "lower");
%! x = rand(rows(A), 1);
%! opts.isreal = 1;
%! opts.issym = 0;
%! opts.maxit = 50000;
%! opts.tol = realmin();
%! [v4, kappa4, flag4] = eigs(@(x) (A \ B) * x, columns(A), M, "lm", opts);
%! kappa4 = diag(kappa4);
%! for k=1:columns(v4)
%!  assert_simple(kappa4(k) * A * v4(:, k), B * v4(:, k), tol2 * norm(B * v4(:, k)));
%! endfor
%! [v5, kappa5, flag5] = eigs(@(x) A \ (B * x), columns(A), M, "lm", opts);
%! kappa5 = diag(kappa5);
%! for k=1:columns(v5)
%!  assert_simple(kappa5(k) * A * v5(:, k), B * v5(:, k), tol2 * norm(B * v5(:, k)));
%! endfor
%! start = tic();
%! [L, P, Q] = chol(A, "lower", "vector");
%! Bperm = B(Q, Q);
%! [v2(Q, :), kappa2, flag2] = eigs(@(x) (L.' \ (L \ (Bperm * x))), columns(A), M, "lm", opts);
%! time_lambda2 += tic() - start;
%! assert_simple(flag2, 0);
%! kappa2 = diag(kappa2);
%! [lambda2, idx_lambda2] = sort(1./kappa2);
%! kappa2 = kappa2(idx_lambda2);
%! v2 = v2(:, idx_lambda2);
%! for k=1:columns(v2)
%!  assert_simple(kappa2(k) * A * v2(:, k), B * v2(:, k), tol2 * max([norm(kappa2(k) * A * v2(:, k)), norm(B * v2(:, k))]));
%!  assert_simple(A * v2(:, k), lambda2(k) * B * v2(:, k), tol2 * max([norm(lambda2(k) * B * v2(:, k)), norm(A * v2(:, k))]));
%! endfor
%! assert_simple(lambda, lambda3, tol1*max(max(abs(lambda))));
%! assert_simple(lambda, lambda2, tol1*max(max(abs(lambda))));
%! for k=1:columns(v)
%!  assert_simple(A * v(:, k), lambda(k) * B * v(:, k), tol1 * norm(A * v(:, k)));
%!  assert_simple(A * v3(:, k), lambda3(k) * B * v3(:, k), tol1 * norm(A * v3(:, k)));
%!  assert_simple(A * v2(:, k), lambda2(k) * B * v2(:, k), tol2 * norm(A * v2(:, k)));
%! endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
