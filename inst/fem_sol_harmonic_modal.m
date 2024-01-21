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
## @deftypefn {Function File} @var{U} = fem_sol_harmonic_modal(@var{h}, @var{lambda}, @var{Phi}, @var{omega})
## Compute the harmonic solution based on eigenvalues @var{lambda}, normalized mode shapes @var{Phi} and modal excitation @var{h}
##
## @var{h} @dots{} Modal excitation returned from fem_sol_modes_scale
##
## @var{lambda} @dots{} Complex eigenvalues returned from fem_sol_modal_damped
##
## @var{Phi} @dots{} Complex mode shapes returned from fem_sol_modal_damped
##
## @var{omega} @dots{} Angular velocity of harmonic excitation @var{h}
##
## @seealso{fem_sol_modal_damped}
## @end deftypefn

function U = fem_sol_harmonic_modal(h, lambda, Phi, omega)
  if (nargout > 1 || nargin ~= 4)
    print_usage();
  endif

  if (~(rows(h) == columns(Phi) && columns(Phi) == numel(lambda)))
    error("size of h, lambda and Phi is not consistent");
  endif

  if (rows(omega) ~= 1)
    error("omega must be a row vector");
  endif

  U = zeros(rows(Phi), numel(omega), columns(h));

  for j=1:columns(h)
    for i=1:columns(Phi)
      qi = h(i, j) ./ (lambda(i) - 1j * omega);
      U(:, :, j) += Phi(:, i) * qi;
    endfor
  endfor
endfunction

