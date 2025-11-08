## Copyright (C) 2011(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{X} = fem_post_def_nodal(@var{mesh}, @var{dof_map}, @var{U})
## @deftypefnx {} [@var{def}, @var{Phi}] = fem_post_def_nodal(@dots{})
## Map a solution vector or load vector @var{U} back to nodal unknowns @var{X}
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{dof_map} @dots{} Mapping for degrees of freedom
##
## @var{U} @dots{} Solution vector or load vector
##
## @var{X} @dots{} Nodal unknowns (e.g. displacement, temperature, pressure)
##
## @var{Phi} @dots{} Complex velocity potential for fluid structure interaction problems
##
## @end deftypefn

function [X, Phi] = fem_post_def_nodal(mesh, dof_map, U)
  if (nargin ~= 3 || nargout > 2)
    print_usage();
  endif

  switch (dof_map.domain)
    case {FEM_DO_STRUCTURAL, FEM_DO_FLUID_STRUCT}
      num_cols = 6;
    case {FEM_DO_THERMAL, FEM_DO_ACOUSTICS}
      num_cols = 1;
    otherwise
      error("unkown value for dof_map.domain");
  endswitch

  X = fem_post_def_core(dof_map, U, 1, num_cols);
  
  if (nargout >= 2)
    switch (dof_map.domain)
      case FEM_DO_FLUID_STRUCT
        Phi = fem_post_def_core(dof_map, U, 7, 7);
      otherwise
        print_usage();
    endswitch
  endif
endfunction

function X = fem_post_def_core(dof_map, U, col_start, col_end)
  num_cols = col_end - col_start + 1;
  
  X = zeros(rows(dof_map.ndof), num_cols, columns(U));

  j = 0;
  
  for i=col_start:col_end
    dof_idx = dof_map.ndof(:, i);
    idx_dof_idx_ne_0 = find(dof_idx > 0);
    X(idx_dof_idx_ne_0, ++j, :) = U(dof_idx(idx_dof_idx_ne_0), :);
  endfor

  if (num_cols == 1)
    X = reshape(X, [size(X, 1), size(X, 3)]);
  endif  
endfunction
