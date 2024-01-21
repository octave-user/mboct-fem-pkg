## Copyright (C) 2011(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{T}, @var{res}] = fem_cms_constr_mat(@var{C})
## 
## Returns a transformation matrix which can be used for the elimination of algebraic constraint equations.
##
## The set of algebraic constraint equations is defined as @var{C} * @var{U} == 0
##
## And the resulting reduced set of unknowns is defined as @var{U} = @var{T} * @var{q} and thus @var{C} * @var{T} * @var{q} == 0
## @seealso{fem_cms_constr_elim}
## @end deftypefn

function [T, res] = fem_cms_constr_mat(C)
  if (nargin ~= 1 || nargout > 2)
    print_usage();
  endif
  
  [q, r] = qr(C.');

  idx = find(sum(abs(r), 2) == 0);

  T = sparse(q(:, idx));

  if (nargout >= 2)
    res = full(max(max(abs(C * T))));
  endif
endfunction

