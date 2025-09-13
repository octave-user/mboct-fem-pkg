## Copyright (C) 2011(-2025) Reinhard <octave-user@a1.net>
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

  if (~isfield(options, "disp"))
    options.disp = 0;
  endif

  if (~isfield(options, "maxit"))
    options.maxit = 50000;
  endif

  if (isfield(options, "p"))
    opt_eig.p = options.p;
  endif

  opt_eig.disp = options.disp;
  opt_eig.maxit = options.maxit;
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
    error("eigs failed to converge: max(err)=%g > tol=%g", max(err), options.tolerance);
  endif
endfunction
