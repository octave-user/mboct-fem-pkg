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
## @deftypefn {Function File} [@var{dgen}, @var{kgen}, @var{rgen}] = fem_sol_modes_scale2(@var{M}, @var{D}, @var{K}, @var{Phi}, @var{R})
## Scale mode shapes returned from fem_sol_modal for use with fem_sol_harmonic_modal2.
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
## @seealso{fem_sol_harmonic_modal2, fem_sol_modal}
## @end deftypefn

function [dgen, kgen, rgen] = fem_sol_modes_scale2(M, D, K, Phi, R)
  mgen = dgen = kgen = zeros(1, columns(Phi));

  for i=1:columns(Phi)
    mgen(i) = Phi(:, i).' * M * Phi(:, i);
    dgen(i) = Phi(:, i).' * D * Phi(:, i);
    kgen(i) = Phi(:, i).' * K * Phi(:, i);
  endfor

  ## divide by mgen, So there is no need to store mgen
  dgen ./= mgen;
  kgen ./= mgen;

  if (nargout >= 3)
    rgen = diag(1 ./ mgen) * (Phi.' * R);
  endif
endfunction
