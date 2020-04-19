## Copyright (C) 2011(-2019) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{sol}] = fem_sol_static(@var{mesh}, @var{dof_map}, @var{mat_ass})
## @deftypefnx {} [@dots{}, @var{U}] = fem_sol_static(@dots{}, @var{options})
##
## Solve a linear static problem and return a structure with nodal displacements
##
## @var{mesh} @dots{} finite element mesh data structure
##
## @var{dof_map} @dots{} degree of freedom mapping
##
## @var{mat_ass} @dots{} data structure containing the global stiffness matrix @var{K} and load vectors @var{R}
##
## @var{options} @dots{} options passed to the linear solver fem_sol_linsolve
##
## @end deftypefn

function [sol, U] = fem_sol_static(mesh, dof_map, mat_ass, options)
  if (nargin < 3 || nargin > 4 || nargout > 2)
    print_usage();
  endif

  if (nargin < 4)
    options = struct();
  endif

  U = fem_sol_linsolve(mat_ass.K, mat_ass.R, options);
  sol.def = fem_post_def_nodal(mesh, dof_map, U);
endfunction
