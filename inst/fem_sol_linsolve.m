## Copyright (C) 2011(-2021) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{U} = fem_sol_linsolve(@var{K}, @var{R})
## @deftypefnx {} @dots{} = fem_sol_linsolve(@dots{}, @var{options})
##
## Choose an appropriate linear solver and solve  @var{K} * @var{U} = @var{R}
##
## @var{K} @dots{} stiffness matrix
##
## @var{R} @dots{} matrix of load vectors
##
## @var{U} @dots{} matrix of displacement vectors
##
## @var{options}.number_of_threads @dots{} set the number of threads for parallel solution
##
## @var{options}.refine_max_iter @dots{} set the maximum number of iterations for refinement of the solution
##
## @end deftypefn

function U = fem_sol_linsolve(K, R, options)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  Kfact = fem_sol_factor(K, options);

  U = Kfact \ R;
endfunction

%!test
%! ### TEST 1
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!   for s=[false,true]
%!     for N=[2,10,20,50,100]
%!       for i=1:10
%!         A = gallery("Poisson", N);
%!         b = 2 * rand(rows(A), 3) - 1;
%!         for i=1:numel(solvers)
%!           options.solver = solvers{i};
%!           options.refine_max_iter = int32(10);
%!           options.number_of_threads = int32(1);
%!           options.pre_scaling = s;
%!           options.epsilon_refinement = eps^0.9;
%!           options.verbose = true;
%!           x = fem_sol_linsolve(A, b, options);
%!           f = max(norm(A * x - b, "cols") ./ norm(A * x + b, "cols"));
%!           assert(f <= options.epsilon_refinement);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect

%!test
%! ## TEST 2
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     K = gallery("Poisson", N);
%!     R = rand(columns(K), 10);
%!     solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!     tol = sqrt(eps);
%!     for i=1:numel(solvers)
%!       options.solver = solvers{i};
%!       options.refine_max_iter = int32(10);
%!       options.pre_scaling = s;
%!       U = fem_sol_linsolve(K, R, options);
%!       f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%!       assert(f < tol);
%!     endfor
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect

%!test
%! ## TEST 3
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     K = gallery("Poisson", N);
%!     K = complex(K, rand() * K);
%!     K += K';
%!     R = complex(rand(columns(K), 10), rand(columns(K), 10));
%!     solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!     tol = sqrt(eps);
%!     for i=1:numel(solvers)
%!       options.solver = solvers{i};
%!       options.refine_max_iter = int32(30);
%!       options.verbose = int32(0);
%!       options.pre_scaling = s;
%!       options.symmetric = false;
%!       U = fem_sol_linsolve(K, R, options);
%!       f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%!       assert(f < tol);
%!     endfor
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect

%!test
%! ## TEST 4
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     K = gallery("Poisson", N);
%!     R = complex(rand(columns(K), 10), rand(columns(K), 10));
%!     solvers = {"pastix", "pardiso", "mumps", "umfpack", "chol", "lu"};
%!     tol = sqrt(eps);
%!     for i=1:numel(solvers)
%!       options.solver = solvers{i};
%!       options.refine_max_iter = int32(10);
%!       options.pre_scaling = s;
%!       U = fem_sol_linsolve(K, R, options);
%!       f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%!       assert(f < tol);
%!     endfor
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect

%!test
%! ## TEST 5
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   solvers = {"pastix", "mumps", "pardiso", "umfpack", "chol", "lu"};
%!   tol = sqrt(eps);
%!   max_f = 0;
%!   for s=[false, true]
%!   for N=[2,10,100,200,300]
%!     for i=1:10
%!       A = 2 * rand(N, N) - 1;
%!       A += 2j * rand(N, N) - 1;
%!       A += A';
%!       b = complex(2 * rand(N, 20) - 1, 2 * rand(N, 20) - 1);
%!       for i=1:numel(solvers)
%!         options.solver = solvers{i};
%!         options.refine_max_iter = int32(100);
%!         options.number_of_threads = int32(1);
%!         options.pre_scaling = s;
%!         options.symmetric = false;
%!         options.verbose = int32(0);
%!         x = fem_sol_linsolve(A, b, options);
%!         f = max(norm(A * x - b, "cols") ./ norm(A * x + b, "cols"));
%!         assert(f < tol);
%!         max_f = max(f, max_f);
%!       endfor
%!     endfor
%!   endfor
%!   endfor
%!   assert(max_f < tol);
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect

%!test
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
%!         assert(x2, x1, tol * norm(x1));
%!         assert(A * x2, b, tol * norm(A*x2 + b));
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect

%!demo
%! N = 100;
%! K = gallery("Poisson", N);
%! R = linspace(0, 1, columns(K)).';
%! options.solver = "chol";
%! U = fem_sol_linsolve(K, R, options);
%! f = max(norm(K * U - R, "cols") ./ norm(K * U + R, "cols"));
%! fprintf(stderr, "backward error: %.1e\n", f);
