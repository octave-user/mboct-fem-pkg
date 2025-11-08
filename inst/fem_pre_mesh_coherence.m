## Copyright (C) 2019(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{dx}] = fem_pre_mesh_coherence(@var{mesh}, @var{tol_distance})
##
## Merge all coherent nodes in @var{mesh} with distances from node to node less than @var{tol_distance}
## @end deftypefn

function [mesh, dx] = fem_pre_mesh_coherence(mesh, tol_distance)
  if (nargin ~= 2 || nargout > 2)
    print_usage();
  endif

  node_map_dst = zeros(rows(mesh.nodes), 1, "int32");
  dx = zeros(rows(mesh.nodes), 1);
  inum_nodes = int32(0);

  nodes = zeros(size(mesh.nodes));

  for i=1:rows(mesh.nodes)
    node_idx_search = 1:i - 1;
    dx_i = norm(mesh.nodes(node_idx_search, 1:3) - mesh.nodes(i, 1:3), "rows");
    idx_coherence = find(dx_i < tol_distance);

    if (isempty(idx_coherence))
      node_map_dst(i) = ++inum_nodes;
      nodes(inum_nodes, :) = mesh.nodes(i, :);
    else
      node_map_dst(i) = node_map_dst(idx_coherence(1));
      dx(i) = dx_i(idx_coherence(1));
    endif
  endfor

  mesh.nodes = nodes(1:inum_nodes, :);

  eltype = fem_pre_mesh_elem_type();

  for i=1:numel(eltype)
    if (isfield(mesh.elements, eltype(i).name))
      elem_nodes = getfield(mesh.elements, eltype(i).name);

      for j=1:columns(elem_nodes)
        elem_nodes(:, j) = node_map_dst(elem_nodes(:, j));
      endfor

      mesh.elements = setfield(mesh.elements, eltype(i).name, elem_nodes);
    endif

    if (isfield(mesh.groups, eltype(i).name))
      groups = getfield(mesh.groups, eltype(i).name);

      for j=1:numel(groups)
        groups(j).nodes = unique(node_map_dst(groups(j).nodes));
      endfor

      mesh.groups = setfield(mesh.groups, eltype(i).name, groups);
    endif
  endfor
endfunction
