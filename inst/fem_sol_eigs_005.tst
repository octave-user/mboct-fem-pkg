## fem_sol_eigs.m:05
%!test
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for N=[2, 10, 50, 100]
%!     num_modes = ceil(0.1 * N);
%!     for i=1:20
%!       K = rand(N, N);
%!       K *= K.';
%!       M = rand(N, N);
%!       M *= M.';
%!       [U, lambda] = fem_sol_eigs(K, M, num_modes);
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!  rand("state", state);
%! end_unwind_protect
