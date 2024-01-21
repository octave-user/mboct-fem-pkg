## fem_sol_eigs.m:02
%!test
%! rand("seed", 0);
%! tol = sqrt(eps);
%! N = 10;
%! k = 3;
%! K = rand(N, N);
%! M = rand(N, N);
%! K *= K.';
%! M *= M.';
%! [L, P] = chol(sparse(M), "lower");
%! sigma = 1;
%! opts.isreal = 1;
%! opts.issym = 1;
%! A = L \ (K / L.');
%! eigsfunc = @(x) A * x;
%! [Psi, lambda, flag] = eigs(eigsfunc, columns(A), k, "lm", opts);
%! if (flag ~= 0)
%! error("eigs failed with status %d", flag);
%! endif
%! v = L.' \ Psi;
%! for i=1:columns(v)
%!  a = K * v(:, i);
%!  b = lambda(i,i) * M * v(:, i);
%!  assert_simple(a, b, tol * norm(abs(a)+abs(b)));
%! endfor
