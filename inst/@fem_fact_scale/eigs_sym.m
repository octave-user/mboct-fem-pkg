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
## @deftypefn {Function File} @var{bsym} = eigs_sym(@var{Afact})
## Return true if the function eigs_func(@var{Afact}, @var{M}, @var{x}) defines a symmetric problem.
## Currently only fem_fact_chol defines a symmetric problem.
## @seealso{fem_sol_eigs}
## @end deftypefn

function bsym = eigs_sym(Afact)
  narginchk(1, 1);
  bsym = eigs_sym(Afact.ASfact);
endfunction
