## Copyright (C) 2011(-2022) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{U}, @var{lambda}, @var{err}] = fem_sol_eigs(@var{K}, @var{M}, @var{N}, @var{rho}, @var{tolerance}, @var{algorithm}, @var{solver}, @var{num_threads}, @var{problem})
## @deftypefnx {} [@var{@dots{}}] = fem_sol_eigs(@var{K}, @var{M}, @var{N}, @var{options})
## Solve the general eigenvalue problem K * U = lambda^2 * M * U or K * U = -lambda * M * U.
##
## @var{K} @dots{} Assembled stiffness matrix. @var{K} is allowed to be singular provided that @var{K} - @var{rho} * @var{M} is regular.
##
## @var{M} @dots{} Assembled mass matrix. @var{M} is allowed to be singular.
##
## @var{N} @dots{} Number of modes to compute
##
## @var{tolerance} @dots{} Tolerance for verification of modes
##
## @var{rho} @dots{} Shift for eigenvalues (e.g. needed if @var{K} is singular)
##
## @var{algorithm} @dots{} Algorithm to use: One of ("generic", "unsymmetric", "shift-invert", "symmetric-inverse").
##
## @var{solver} @dots{} Linear solver to use: One of ("pastix", "pardiso", "mumps", "umfpack", "chol", "lu", "mldivide").
##
## @var{num_threads} @dots{} Number of threads to use for the linear solver.
##
## @var{problem} @dots{} May be one of ("structural", "acoustic", "thermal", "pressure", "buckling")
##
## @var{options} @dots{} Alternative form to pass all or some of the arguments above as a struct
##
## @seealso{fem_sol_modal}
## @end deftypefn

function [U, lambda, err] = fem_sol_eigs(K, M, N, varargin)
  if (nargin < 3 || nargin > 9 || nargout > 3)
    print_usage();
  endif

  if (~(issquare(K) && issquare(M)))
    error("K and M must be square matrices");
  endif

  if (~(isreal(K) && isreal(M)))
    error("K and M must be real");
  endif

  if (~all(size(K) == size(M)))
    error("K and M must have the same size");
  endif

  if (numel(varargin) >= 1)
    if (isstruct(varargin{1}))
      options = varargin{1};
    else
      options.rho = varargin{1};
      
      if (numel(varargin) >= 2)
        options.tolerance = varargin{2};
      endif

      if (numel(varargin) >= 3)
        options.algorithm = varargin{3};
      endif

      if (numel(varargin) >= 4)
        options.solver = varargin{4};
      endif

      if (numel(varargin) >= 5)
        options.number_of_threads = varargin{5};
      endif

      if (numel(varargin) >= 6)
        options.problem = varargin{6};
      endif
    endif
  else
    options = struct();
  endif
  
  if (~isfield(options, "tolerance"))
    options.tolerance = 1e-4;
  endif

  if (~isfield(options, "problem"))
    options.problem = "structural";
  endif

  if (~isfield(options, "algorithm"))
    switch (options.problem)
      case "buckling"
        options.algorithm = "symmetric-inverse";
      otherwise
        options.algorithm = "shift-invert";
    endswitch
  endif

  opt_eig.disp = 0;
  opt_eig.maxit = 50000;
  opt_eig.isreal = true;

  if (~isfield(options, "rho"))
    options.rho = 0;
    Ksh = K;
  else
    Ksh = K - options.rho * M;
  endif

  Kfact = fem_sol_factor(Ksh, options);

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);

    switch (options.algorithm)
      case {"generic", "unsymmetric"}
        SIGMA = "LM";

        opt_eig.issym = eigs_sym(Kfact);

        eigs_callback = @(x) eigs_func(Kfact, M, x);

        iter_eig = int32(0);
        max_iter_eig = int32(10);

        while (true)
          [U, mu, status] = eigs(eigs_callback, columns(M), N, SIGMA, opt_eig);

          if (status ~= 0)
            error("eigs failed to converge");
          endif

          if (isreal(U))
            break;
          endif

          ## If U is not real, call eigs again with a new real starting vector

          if (++iter_eig > max_iter_eig)
            error("eigs returned complex eigenvectors");
          endif

          opt_eig.v0 = real(U(:, 1));
        endwhile

        U = eigs_post(Kfact, U);
      case {"shift-invert", "symmetric-inverse"}
        opt_eig.v0 = rand(columns(M), 1);
        op{1} = @(x) M * x;
        op{2} = @(x) Kfact \ x;
        switch (options.algorithm)
          case "symmetric-inverse"
            SIGMA = "LM";
            op{3} = @(x) Ksh * x;
          otherwise
            SIGMA = 0;
        endswitch
        [U, mu] = eig_sym(op, columns(M), N, SIGMA, opt_eig);
      otherwise
        error("unknown algorithm: \"%s\"", options.algorithm);
    endswitch
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  mu = diag(mu).';
  
  switch (options.algorithm)
    case "shift-invert"
    otherwise
      mu = 1 ./ mu;
  endswitch

  kappa = mu + options.rho; # Apply a shift in order to work if rigid body modes are present

  switch (options.problem)
    case {"structural", "acoustic"}
      lambda = sqrt(-kappa);
      lambda2 = lambda.^2;
    case {"pressure", "buckling", "thermal"}
      lambda2 = lambda = -kappa;
    otherwise
      error("unknown value for parameter options.problem=\"%s\"", options.problem);
  endswitch
  
  [lambda, i_lambda] = sort(lambda);

  lambda2 = lambda2(i_lambda);
  
  s_U_d = zeros(columns(U), 1);

  for i=1:columns(U)
    s_U_d(i) = sqrt(abs(U(:, i).' * M * U(:, i))); ## abs needed for buckling type problems
  endfor

  U /= diag(s_U_d);

  U = U(:, i_lambda);

  err = zeros(columns(U), 1);

  for i=1:columns(U)
    v1 = K * U(:, i);
    v2 = M * U(:, i) * lambda2(i);
    err(i) = norm(v1 + v2) / max(norm(v1), norm(v2));
  endfor

  if (nargout < 3 && max(err) > options.tolerance)
    error("eigs failed to converge: max(err)=%g > tol=%g", max(err), tol);
  endif
endfunction

%!test
%! tol = 1e-4;
%! N = 300;
%! n = 10;
%! toleigs = 0;
%! solvers = {"pastix", "pardiso", "mumps", "umfpack", "lu", "chol", "mldivide"};
%! alg={"symmetric-inverse","shift-invert","unsymmetric"};
%! problems = {"structural", "thermal", "acoustic", "pressure", "buckling"};
%! t = zeros(numel(alg), numel(solvers));
%! for k=1:numel(solvers)
%!   for a=1:numel(alg)
%!     for j=1:2
%!       rand("seed", 0);
%!       for i=1:3
%!         K = sprand(N, N, 0.01) + 10*abs(diag(0.1+rand(N, 1)));
%!         M = sprand(rows(K), columns(K), 0.01) + 10*abs(diag(0.1+rand(rows(K), 1)));
%!         K *= K.';
%!         M *= M.';
%!         Q = symrcm(M + K);
%!         M = M(Q, Q);
%!         K = K(Q, Q);
%!         if j == 2
%!           rho = -0.1 * max(abs(diag(K))) / max(abs(diag(M)));
%!         else
%!           rho = 0;
%!         endif
%!         start = tic();
%!         for l=1:7
%!           switch (l)
%!             case 1
%!               [U, lambda, err] = fem_sol_eigs(K, M, n, rho, toleigs, alg{a}, solvers{k});
%!             otherwise
%!               opt = struct();
%!               if (l > 2)
%!                 opt.rho = rho;
%!               endif
%!               if (l > 3)
%!                 opt.tolerance = toleigs;
%!               endif
%!               if (l > 4)
%!                 opt.algorithm = alg{a};
%!               endif
%!               if (l > 5)
%!                 opt.solver = solvers{k};
%!               endif
%!               if (l > 6)
%!                 opt.number_of_threads = int32(4);
%!               endif
%!               if (l > 7)
%!                 for m=1:numel(problems)
%!                   opt.problem = opt.problem = problems{m};
%!                   [U, lambda, err] = fem_sol_eigs(K, M, n, opt);
%!                   assert(max(err) < tol);
%!                 endfor
%!               else
%!                   [U, lambda, err] = fem_sol_eigs(K, M, n, opt);
%!                   assert(max(err) < tol);
%!               endif
%!           endswitch
%!           t(a, k) += toc(start);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! for i=1:numel(alg)
%!   for j=1:numel(solvers)
%!     fprintf(stderr, "algorithm \"%s:%s\": %.2fs\n", alg{i}, solvers{j}, t(i, j));
%!   endfor
%! endfor

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
%!  assert(a, b, tol * norm(abs(a)+abs(b)));
%! endfor

%!test
%! rand("seed", 0);
%! tol = eps^0.3;
%! N = 10;
%! M = 3;
%! for j=1:2
%! for i=1:10
%! A = rand(N, N) + abs(diag(rand(N, 1)));
%! if j == 1
%! B = eye(columns(A));
%! else
%! B = rand(rows(A), columns(A)) + abs(diag(rand(rows(A), 1)));
%! endif
%! A *= A.';
%! B *= B.';
%! opts.isreal = 1;
%! opts.issym = 1;
%! opts.tol = realmin();
%! opts.maxiter = 5000;
%! [v, lambda] = eig(A, B);
%! [v3, lambda3] = eig(B \ A);
%! [lambda, idx_lambda] = sort(diag(lambda));
%! v = v(:, idx_lambda);
%! [lambda3, idx_lambda3] = sort(diag(lambda3));
%! v3 = v3(:, idx_lambda3);
%! [L, P] = chol(B, "lower");
%! assert(P,0);
%! A2 = (L \ A) / L.';
%! assert(issymmetric(A2, tol*norm(A2)));
%! [L2, P2] = chol(A2, "lower");
%! assert(P2,0);
%! [v2, lambda2, flag2] = eigs(@(y) L2.' \ (L2 \ y), columns(A2), M, "sm", opts);
%! [lambda2, idx_lambda2] = sort(diag(lambda2));
%! v2 = v2(:, idx_lambda2);
%! v2 = L.' \ v2;
%! assert(flag2, 0);
%! assert(lambda(1:3), lambda3(1:M), tol*max(max(abs(lambda(1:M)))));
%! assert(lambda(1:M), lambda2, tol*max(max(abs(lambda(1:M)))));
%! for k=1:M
%!  assert(A * v(:, k), lambda(k) * B * v(:, k), tol * norm(A * v(:, k)));
%!  assert(A * v3(:, k), lambda3(k) * B * v3(:, k), tol * norm(A * v3(:, k)));
%!  assert(A * v2(:, k), lambda2(k) * B * v2(:, k), tol * norm(A * v2(:, k)));
%! endfor
%! endfor
%! endfor

%!test
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
%!  assert(kappa4(k) * A * v4(:, k), B * v4(:, k), tol2 * norm(B * v4(:, k)));
%! endfor
%! [v5, kappa5, flag5] = eigs(@(x) A \ (B * x), columns(A), M, "lm", opts);
%! kappa5 = diag(kappa5);
%! for k=1:columns(v5)
%!  assert(kappa5(k) * A * v5(:, k), B * v5(:, k), tol2 * norm(B * v5(:, k)));
%! endfor
%! start = tic();
%! [L, P, Q] = chol(A, "lower", "vector");
%! Bperm = B(Q, Q);
%! [v2(Q, :), kappa2, flag2] = eigs(@(x) (L.' \ (L \ (Bperm * x))), columns(A), M, "lm", opts);
%! time_lambda2 += tic() - start;
%! assert(flag2, 0);
%! kappa2 = diag(kappa2);
%! [lambda2, idx_lambda2] = sort(1./kappa2);
%! kappa2 = kappa2(idx_lambda2);
%! v2 = v2(:, idx_lambda2);
%! for k=1:columns(v2)
%!  assert(kappa2(k) * A * v2(:, k), B * v2(:, k), tol2 * max([norm(kappa2(k) * A * v2(:, k)), norm(B * v2(:, k))]));
%!  assert(A * v2(:, k), lambda2(k) * B * v2(:, k), tol2 * max([norm(lambda2(k) * B * v2(:, k)), norm(A * v2(:, k))]));
%! endfor
%! assert(lambda, lambda3, tol1*max(max(abs(lambda))));
%! assert(lambda, lambda2, tol1*max(max(abs(lambda))));
%! for k=1:columns(v)
%!  assert(A * v(:, k), lambda(k) * B * v(:, k), tol1 * norm(A * v(:, k)));
%!  assert(A * v3(:, k), lambda3(k) * B * v3(:, k), tol1 * norm(A * v3(:, k)));
%!  assert(A * v2(:, k), lambda2(k) * B * v2(:, k), tol2 * norm(A * v2(:, k)));
%! endfor
%! endfor

%!demo
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
