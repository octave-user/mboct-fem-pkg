## Copyright (C) 2021(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{U} = eigs_post(@var{Afact}, @var{X})
## Converts eigenvectors @var{X}, returned from eigs, to the solution @var{U} of K * U = lambda^2 * M * U.
## @seealso{fem_sol_eigs}
## @end deftypefn

function U = eigs_post(fact, X)
  narginchk(2, 2);
  
  U = X;
endfunction
