## Copyright (C) 2018(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{rbe2}] = fem_pre_mesh_rbe2_from_surf(@var{mesh}, @var{group_id}, @var{master_node_idx})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_rbe2_from_surf(@dots{}, @var{elem_type})
##
## Builds rbe2 elements from specified groups of surface elements
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_id} @dots{} array of group id-numbers used to identify surface elements with slave nodes
##
## @var{master_node_idx} @dots{} array of master node indices for rbe3 elements
##
## @var{elem_type} @dots{} the element type addressed by @var{group_id} (e.g. "tria6", "iso4", "quad8")
##
## @var{rbe2} @dots{} rigid body constraints
##
## @end deftypefn

function rbe2 = fem_pre_mesh_rbe2_from_surf(mesh, group_id, master_node_idx, elem_type)
  if (nargin < 3 || nargin > 4 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    elem_type = "tria6";
  endif

  if (~iscell(elem_type))
    elem_type = {elem_type};
  endif

  rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, group_id, master_node_idx, elem_type);

  num_rbe2 = int32(0);

  for j=1:numel(rbe3)
    num_rbe2 += columns(rbe3(j).nodes) - 1;
  endfor

  empty_cell = cell(1, num_rbe2);
  rbe2 = struct("nodes", empty_cell);
  num_rbe2 = int32(0);

  for j=1:numel(rbe3)
    for i=1:numel(rbe3(j).nodes) - 1
      rbe2(++num_rbe2).nodes = rbe3(j).nodes([1, i + 1]);
    endfor
  endfor
endfunction
