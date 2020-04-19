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
## @deftypefn {Function File} [@var{sol}, @var{err}, @var{U}] = fem_sol_modal(@var{mesh}, @var{dof_map}, @var{mat_ass}, @var{N}, @var{varargin})
## Compute undamped natural frequencies and mode shapes.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{dof_map} @dots{} Degree of freedom mapping
##
## @var{mat_ass}.K @dots{} Nodal stiffness matrix
##
## @var{mat_ass}.M @dots{} Nodal mass matrix
##
## @var{N} @dots{} Number of modes to compute
##
## @var{sol}.def @dots{} nodal solution for mode shapes
##
## @var{sol}.f @dots{} undamped natural frequencies
##
## @var{err} @dots{} backward error
##
## @var{U} @dots{} mode shapes
##
## @seealso{fem_sol_eigs}
## @end deftypefn

function [sol, err, U] = fem_sol_modal(mesh, dof_map, mat_ass, N, varargin)
  if (nargin < 4 || nargout > 3)
    print_usage();
  endif
  
  [U, lambda, err] = fem_sol_eigs(mat_ass.K, mat_ass.M, N, varargin{:});
  
  sol.lambda = lambda;
  sol.f = imag(lambda) / (2 * pi);
  sol.def = fem_post_def_nodal(mesh, dof_map, U);
endfunction  
