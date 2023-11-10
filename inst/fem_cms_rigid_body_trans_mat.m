## Copyright (C) 2011(-2022) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{A} = fem_cms_rigid_body_trans_mat(@var{Xn})
##
## Build transformation matrix for rigid body motion of a mesh with nodal coordinates @var{Xn}
##
## @end deftypefn

function A = fem_cms_rigid_body_trans_mat(Xn)
  A = zeros(3 * rows(Xn), 6);

  for j=1:3
    A(j:3:end, j) = 1;
  endfor

  A(1:3:end, 5) = Xn(:, 3);
  A(2:3:end, 4) = -Xn(:, 3);
  A(1:3:end, 6) = -Xn(:, 2);
  A(3:3:end, 4) = Xn(:, 2);
  A(2:3:end, 6) = Xn(:, 1);
  A(3:3:end, 5) = -Xn(:, 1);
endfunction
