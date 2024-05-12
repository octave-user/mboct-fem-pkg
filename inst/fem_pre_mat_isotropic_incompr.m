## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{C} = fem_pre_mat_isotropic_incompr(@var{E}, @var{nu})
## Compute the constitutive law matrix @var{C} for incompressible isotropic elasticity
##
## @var{E} @dots{} Young's modulus
##
## @var{nu} @dots{} Poisson ratio
##
## @end deftypefn

function C = fem_pre_mat_isotropic_incompr(E, nu)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif
  
  mu = E / (2 * (1 + nu));
  kappa = E / (3 * (1 - 2 * nu));
  gamma = 3 * nu / (1 + nu);
  
  C =  [diag(2 * mu * ones(1, 3)),           zeros(3, 3), -gamma * ones(3, 1);
                      zeros(3, 3), diag(mu * ones(1, 3)),         zeros(3, 1);
                      -ones(1, 3),           zeros(1, 3),          -1 / kappa];
endfunction     

