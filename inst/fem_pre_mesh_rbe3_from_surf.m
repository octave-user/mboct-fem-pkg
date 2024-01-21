## Copyright (C) 2018(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{rbe3}] = fem_pre_mesh_rbe3_from_surf(@var{mesh}, @var{group_id}, @var{master_node_idx})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_rbe3_from_surf(@dots{}, @var{elem_type})
##
## Builds rbe3 elements from specified groups of tria6 elements and assigns weighting factors proportional to equivalent surface area
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_id} @dots{} array of group id-numbers used to identify surface elements with slave nodes
##
## @var{master_node_idx} @dots{} array of master node indices for rbe3 elements
##
## @var{elem_type} @dots{} the element type addressed by @var{group_id} (e.g. "tria6" or "iso4")
##
## @end deftypefn

function rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, group_id, master_node_idx, elem_type)
  if (nargin < 3 || nargin > 4 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    elem_type = "tria6";
  endif

  if (~iscell(elem_type))
    elem_type = {elem_type};
  endif

  empty_cell = cell(1, numel(group_id));

  rbe3 = struct("nodes", empty_cell, "weight", empty_cell, "elements", empty_cell);

  if (numel(group_id) ~= numel(master_node_idx))
    error("numel(group_id) does not match numel(master_node_idx)");
  endif

  ##msh_groups = getfield(mesh.groups, elem_type);

  [F, group_idx, weight, node_id] = fem_pre_mesh_nodal_pressure_load(mesh, group_id, elem_type);

  for i=1:numel(node_id)
    rbe3(i).nodes = int32([master_node_idx(i), node_id{i}]);
    rbe3(i).weight = weight{i}.';
    rbe3(i).elements = struct();

    for j=1:numel(elem_type)
      if (~isfield(mesh.groups, elem_type{j}))
        continue;
      endif
      msh_groups = getfield(mesh.groups, elem_type{j});
      rbe3(i).elements = setfield(rbe3(i).elements, elem_type{j}, [msh_groups(group_idx{i, j}).elements]);
    endfor
  endfor
endfunction

