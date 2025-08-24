%!test
%! M = 10;
%! K = 500;
%! D = 30;
%! Phi = 1;
%! R = 1;
%! zeta_ref = 1 / 2 * D / sqrt(K * M);
%! [dgen, kgen, rgen, lambda, zeta] = fem_sol_modes_scale2(M, D, K, Phi, R);
%! assert(zeta, zeta_ref);

%!test
%! tan_delta = 0.01e-2;
%! M = 10;
%! K = 500 * (1 + 1j * tan_delta);
%! D = 0;
%! Phi = 1;
%! R = 1;
%! [dgen, kgen, rgen, lambda, zeta] = fem_sol_modes_scale2(M, D, K, Phi, R);
%! assert(2 * zeta, tan_delta, sqrt(eps) * zeta);
