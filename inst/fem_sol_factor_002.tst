## fem_sol_factor.m:02
%!test
%! s = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"chol", "lu", "mldivide", "pastix", "pardiso", "mumps", "umfpack"};
%!   N = 20;
%!   for k=1:numel(solvers)
%!     if (~fem_sol_check_func(solvers{k}))
%!       warning("linear solver %s is not installed", solvers{k});
%!       continue;
%!     endif
%!     for i=1:100
%!       for j=1:3
%!         A = gallery("Poisson", N);
%!         opts.refine_max_iter = int32(100);
%!         opts.epsilon_refinement = eps^0.8;
%!         opts.solver = solvers{k};
%!         switch (j)
%!           case 1
%!             Asym = A;
%!             opts.symmetric = false;
%!           case {2, 3}
%!             switch (opts.solver)
%!             case {"umfpack", "lu", "mldivide"}
%!               Asym = A;
%!             otherwise
%!               [r, c, d] = find(A);
%!               switch (opts.solver)
%!                 case "chol"
%!                   idx = find(r <= c); ## chol uses only the upper triangular part
%!                 otherwise
%!                 switch (j)
%!                 case 2
%!                   idx = find(r >= c);
%!                 case 3
%!                   idx = find(r <= c);
%!                 endswitch
%!               endswitch
%!               Asym = sparse(r(idx), c(idx), d(idx), rows(A), columns(A));
%!             endswitch
%!             opts.symmetric = true;
%!         endswitch
%!         Afact = fem_sol_factor(Asym, opts);
%!         B = rand(rows(A), 30);
%!         X = Afact \ B;
%!         assert_simple(max(norm(A * X - B, "cols") ./ norm(A * X + B, "cols")) <= opts.epsilon_refinement);
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", s);
%! end_unwind_protect
