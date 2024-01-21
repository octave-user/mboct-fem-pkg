## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{U}, @var{lambda}] = fem_sol_eigsd(@var{K}, @var{D}, @var{M}, @var{N}, @var{opt_linsol})
## Solve the general eigenvalue problem (lambda^2 * M + lambda * D + K) * U = 0
##
## @var{K} @dots{} Assembled stiffness matrix. @var{K} must not be singular
##
## @var{D} @dots{} Assembled damping matrix. @var{D} is allowed to be singular
##
## @var{M} @dots{} Assembled mass matrix. @var{M} is allowed to be singular
##
## @var{N} @dots{} Number of modes to compute
##
## @var{options} @dots{} Options passed to fem_sol_factor
##
## @seealso{fem_sol_modal_damped}
## @end deftypefn

function [U, lambda] = fem_sol_eigsd(K, D, M, N, opt_linsol, opts)
  if (nargin < 4 || nargin > 6 || nargout > 2)
    print_usage();
  endif

  if (nargin < 5)
    opt_linsol = struct();
  endif

  if (nargin < 6)
    opts = struct();
  endif

  if (~(ismatrix(K) && ismatrix(D) && ismatrix(M) && isscalar(N) && isstruct(opt_linsol)))
    print_usage();
  endif

  if (~(issquare(K) && issquare(D) && issquare(M)))
    error("K, D and M must be square");
  endif

  if (~all(size(K) == size(D) & size(K) == size(M)))
    error("K, D and must have the same size");
  endif

  Kfact = fem_sol_factor(K, opt_linsol);

  SIGMA = "LI";

  if (~isfield(opts, "disp"))
    opts.disp = 0;
  endif

  if (~isfield(opts, "maxit"))
    opts.maxit = 50000;
  endif

  if (~isfield(opts, "p"))
    opts.p = 2 * N;
  endif

  if (~isfield(opts, "output"))
    opts.output = "reduced";
  endif

  opts.isreal = isreal(K) && isreal(M) && isreal(D);
  opts.issym = false;

  callback = @(x) eigs_callback(Kfact, D, M, x);

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);
    [U, kappa, status] = eigs(callback, 2 * columns(M), N, SIGMA, opts);
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  lambda = -1 ./ diag(kappa).';

  switch (opts.output)
    case "reduced"
      U = U(columns(M) + (1:columns(M)), :);
    case "full"
    otherwise
      warning("unknown option output=\"%s\"", opts.output);
  endswitch
endfunction

function y = eigs_callback(Kfact, D, M, x)
  ndof = columns(M);
  idx1 = 1:ndof;
  idx2 = ndof + idx1;

  y = zeros(size(x));

  y(idx1, :) = -x(idx2, :);
  y(idx2, :) = Kfact \ (M * x(idx1, :) + D * x(idx2, :));
endfunction

