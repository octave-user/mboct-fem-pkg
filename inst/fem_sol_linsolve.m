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
## @deftypefn {Function File} @var{U} = fem_sol_linsolve(@var{K}, @var{R})
## @deftypefnx {} @dots{} = fem_sol_linsolve(@dots{}, @var{options})
##
## Choose an appropriate linear solver and solve  @var{K} * @var{U} = @var{R}
##
## @var{K} @dots{} stiffness matrix
##
## @var{R} @dots{} matrix of load vectors
##
## @var{U} @dots{} matrix of displacement vectors
##
## @var{options}.number_of_threads @dots{} set the number of threads for parallel solution
##
## @var{options}.refine_max_iter @dots{} set the maximum number of iterations for refinement of the solution
##
## @end deftypefn

function U = fem_sol_linsolve(K, R, options)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  Kfact = fem_sol_factor(K, options);

  U = Kfact \ R;
endfunction

