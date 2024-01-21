## fem_sol_linsolve.m:07
%!test
%! N = 100;
%! K = gallery("Poisson", N);
%! R = linspace(0, 1, columns(K)).';
%! options.solver = "chol";
%! U = fem_sol_linsolve(K, R, options);
%! f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%! fprintf(stderr, "backward error: %.1e\n", f);
