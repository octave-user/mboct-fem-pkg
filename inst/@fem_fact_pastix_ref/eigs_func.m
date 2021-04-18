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
## @deftypefn {Function File} @var{f} = eigs_func(@var{Kfact}, @var{M}, @var{x})
## Callback function for fem_sol_eigs.
## @seealso{fem_sol_eigs}
## @end deftypefn

function f = eigs_func(fact, M, x)
  narginchk(3, 3);
  
  f = fact \ (M * x);
endfunction
