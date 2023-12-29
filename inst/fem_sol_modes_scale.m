## Copyright (C) 2023(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{Phi}, @var{h}] = fem_sol_modes_scale(@var{M}, @var{K}, @var{lambda}, @var{Phi}, @var{R})
## Scale mode shapes returned from fem_sol_modal_damped for use with fem_sol_harmonic_modal.
## Eigenvalues and mode shapes may be complex.
##
## @var{M} @dots{} Nodal symmetric mass matrix
##
## @var{K} @dots{} Nodal symmetric stiffness matrix
##
## @var{lambda} @dots{} Complex eigenvalues
##
## @var{Phi} @dots{} Complex mode shapes
##
## @var{R} @dots{} Nodal load vector
##
## @var{h} @dots{} Modal load vector
##
## @seealso{fem_sol_modal_damped, fem_sol_harmonic_modal}
## @end deftypefn

function [Phi, h] = fem_sol_modes_scale(M, K, lambda, Phi, R)
  if (~((nargin == 4 && nargout == 1) || (nargin == 5 && nargout == 2)))
    print_usage();
  endif

  if (~(ismatrix(M) && ismatrix(K) && isvector(lambda) && ismatrix(Phi)))
    error("input arguments must be matrices");
  endif

  if (~(issquare(M) && issquare(K) && columns(M) == columns(K)))
    error("M and K must be square matrices and must have the same size");
  endif

  if (~issymmetric(M))
    warning("M is not symmetric");
  endif

  if (~issymmetric(K))
    warning("K is not symmetric");
  endif

  if (~(numel(lambda) == columns(Phi) && rows(Phi) == rows(M)))
    error("invalid size for lambda and Phi");
  endif

  b = zeros(columns(Phi), 1);

  for i=1:columns(Phi)
    b(i) = -lambda(i)^2 * Phi(:, i).' * M * Phi(:, i) + Phi(:, i).' * K * Phi(:, i);
  endfor

  Phi *= diag(1 ./ sqrt(b));

  if (nargout >= 2)
    h = diag(lambda) * Phi.' * R;
  endif
endfunction
