## fem_sol_eigsd.m:01
%!test
%! try
%! rndstate = rand("state");
%! unwind_protect
%!   for j=1:1000
%!     m1 = 10 * rand() + 1;
%!     m2 = 20 * rand() + 1;
%!     m3 = 0;
%!     m4 = 0;
%!     k1 = 100 * rand() + 1;
%!     k2 = 5000 * rand() + 1;
%!     k3 = 1 + rand();
%!     k4 = 1 + rand();
%!     d1 = 10 * rand();
%!     d2 = 10 * rand();
%!     d3 = 0;
%!     d4 = 0;
%!     M = diag([m1, m2, m3, m4]);
%!     K = diag([k1, k2, k3, k4]);
%!     D = diag([d1, d2, d3, d4]);
%!     [~, idx] = sort(rand(columns(M), 1));
%!     M = M(idx, idx);
%!     K = K(idx, idx);
%!     D = D(idx, idx);
%!     N = 4;
%!     rho = 0;
%!     tol = 0;
%!     opt_linsol.solver = "umfpack";
%!     opt_linsol.number_of_threads = int32(1);
%!     opts.disp = 0;
%!     [U, lambda] = fem_sol_eigsd(K, D, M, N, opt_linsol, opts);
%!     lambda_ref = [-d1 / (2 * m1) + [1, -1] * sqrt((d1 / (2 * m1))^2 - k1 / m1), ...
%!                   -d2 / (2 * m2) + [1, -1] * sqrt((d2 / (2 * m2))^2 - k2 / m2)];
%!     tol_lambda = eps^0.8;
%!     tol_U = eps^0.7;
%!     assert_simple(sort(lambda), sort(lambda_ref), tol_lambda * norm(lambda_ref));
%!     for i=1:columns(U)
%!       v1 = lambda(i)^2 * (M * U(:, i)) + lambda(i) * (D * U(:, i));
%!       v2 = -K * U(:, i);
%!       assert_simple(v1, v2, max(norm(v1), norm(v2)) * tol_U);
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", rndstate);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
