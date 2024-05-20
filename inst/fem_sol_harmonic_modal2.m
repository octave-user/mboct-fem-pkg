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
## @deftypefn {Function File} @var{U} = fem_sol_harmonic_modal2(@var{dgen}, @var{kgen}, @var{rgen}, @var{Phi}, @var{omega})
## Compute the harmonic solution based on modal damping @var{dgen}, modal stiffness @var{kgen} modal excitation @var{rgen} and mode shapes @var{Phi}.
## It is assumed, that the modal mass @var{mgen} is equal to ones(size(kgen)).
##
## @var{dgen} @dots{} Modal damping returned from fem_sol_modes_scale2
##
## @var{kgen} @dots{} Modal stiffness returned from fem_sol_modes_scale2
##
## @var{rgen} @dots{} Modal excitation returned from fem_sol_modes_scale2
##
## @var{Phi} @dots{} Real mode shapes returned from fem_sol_modal
##
## @var{omega} @dots{} Angular velocity of harmonic excitation @var{rgen}
##
## @seealso{fem_sol_modal, fem_sol_modes_scale2}
## @end deftypefn

function U = fem_sol_harmonic_modal2(dgen, kgen, rgen, Phi, omega)
  U = zeros(rows(Phi), numel(omega), columns(rgen));

  for j=1:columns(rgen)
    qj = rgen(:, j) ./ (-omega.^2 + 1j * omega .* dgen(:) + kgen(:)); ## assume mgen == 1
    U(:, :, j) += Phi * qj;
  endfor
endfunction
