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
## @deftypefn {Function File} @var{def} = fem_post_def_nodal(@var{mesh}, @var{dof_map}, @var{U})
## Map a solution vector or load vector @var{U} back to nodal displacements @var{def}
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{dof_map} @dots{} Mapping for degrees of freedom
##
## @var{U} @dots{} Solution vector or load vector
##
## @end deftypefn

function def = fem_post_def_nodal(mesh, dof_map, U)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  def = zeros(rows(dof_map.ndof), columns(dof_map.ndof), columns(U));
  
  for i=1:columns(dof_map.ndof)
    dof_idx = dof_map.ndof(:, i);
    idx_dof_idx_ne_0 = find(dof_idx > 0);
    def(idx_dof_idx_ne_0, i, :) = U(dof_idx(idx_dof_idx_ne_0), :);
  endfor
endfunction

