## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{C} = fem_pre_mat_isotropic(@var{E}, @var{nu})
## Compute the constitutive law matrix @var{C} for isotropic elasticity
##
## @var{E} @dots{} Young's modulus
##
## @var{nu} @dots{} Poisson ratio
##
## @end deftypefn

function C = fem_pre_mat_isotropic(E, nu)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif
  
  a = nu / (1 - nu);
  b = (1 - 2 * nu) / (2 * (1 - nu));
  c = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu));
  
  C = c * [1, a, a, 0, 0, 0;
           a, 1, a, 0, 0, 0;
           a, a, 1, 0, 0, 0;
           0, 0, 0, b, 0, 0;
           0, 0, 0, 0, b, 0;
           0, 0, 0, 0, 0, b];
endfunction     

%!test ##demo
%! E = 210000;
%! nu = 0.3;
%! format("bank", "local");
%! C = fem_pre_mat_isotropic(E, nu)
