## fem_sol_linsolve.m:06
%!test
%! try
%! ## TEST 6
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   func={"mldivide", ...
%!         "lu" , ...
%!         "chol", ...
%!         "umfpack",  ...
%!         "pastix", ...
%!         "mumps", ...
%!         "pardiso"};
%!   classes={@fem_fact, ...
%!            @fem_fact_lu, ...
%!            @fem_fact_chol, ...
%!            @fem_fact_umfpack, ...
%!            @fem_fact_pastix, ...
%!            @fem_fact_mumps, ...
%!            @fem_fact_pardiso};
%!   warnfunc = false(size(func));
%!   options.refine_max_iter = int32(100);
%!   options.verbose = int32(0);
%!   for k=1:2
%!     for j=1:numel(func)
%!       for i=1:100
%!         switch (k)
%!           case 1
%!             A = rand(10,10);
%!             M = rand(10, 10);
%!           case 2
%!             A = sprand(100,100,0.05) + 5*diag(rand(100,1));
%!             M = sprand(100,100,0.05) + 5*diag(rand(100,1));
%!         endswitch
%!         A *= A.';
%!         M *= M.';
%!         Q = symrcm(A);
%!         A = A(Q, Q);
%!         M = M(Q, Q);
%!         b = rand(columns(A), 5);
%!         if (~fem_sol_check_func(func{j}))
%!           if (~warnfunc(j))
%!             warning("function \"%s\" not found", func{j});
%!             warnfunc(j) = true;
%!           endif
%!           continue;
%!         endif
%!         Afact = feval(classes{j}, A, options);
%!         x1 = A \ b;
%!         x2 = Afact \ b;
%!         tol = eps^0.4;
%!         assert_simple(x2, x1, tol * norm(x1));
%!         assert_simple(A * x2, b, tol * norm(A*x2 + b));
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
