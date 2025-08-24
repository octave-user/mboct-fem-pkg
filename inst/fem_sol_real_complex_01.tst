%!test
%! try
%! s = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"chol", "lu", "mldivide", "pastix", "pardiso", "mumps", "umfpack"};
%!   n = 10;
%!   m = 3;
%!   for m=1:10
%!     for l=[false,true]
%!       opts.symmetric = l;
%!       A = sparse(rand(n, n));
%!       if (opts.symmetric)
%!         A *= A.';
%!       endif
%!       B = rand(n, m) + 1j * rand(n, m);
%!       for k=[false, true]
%!         opts.pre_scaling = k;
%!         for j=[0, 1, 3, 10]
%!           opts.refine_max_iter = int32(j);
%!           for i=1:numel(solvers)
%!             opts.solver = solvers{i};
%!             switch (opts.solver)
%!               case "chol"
%!                 switch (opts.symmetric)
%!                   case false
%!                     continue;
%!                 endswitch
%!             endswitch
%!             Afact = fem_sol_factor(A, opts);
%!             X1 = A \ B;
%!             X2 = A \ real(B) + 1j * (A \ imag(B));
%!             X3 = Afact \ real(B) + 1j * (Afact \ imag(B));
%!             X4 = Afact \ B;
%!             f1 = max(norm(A * X1 - B, "cols")) / max(1,max(norm(A * X1 + B, "cols")));
%!             f2 = max(norm(A * X2 - B, "cols")) / max(1,max(norm(A * X2 + B, "cols")));
%!             f3 = max(norm(A * X3 - B, "cols")) / max(1,max(norm(A * X3 + B, "cols")));
%!             f4 = max(norm(A * X4 - B, "cols")) / max(1,max(norm(A * X4 + B, "cols")));
%!             tol = eps^0.5;
%!             assert_simple(f1 < tol);
%!             assert_simple(f2 < tol);
%!             assert_simple(f3 < tol);
%!             assert_simple(f4 < tol);
%!           endfor
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", s);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
