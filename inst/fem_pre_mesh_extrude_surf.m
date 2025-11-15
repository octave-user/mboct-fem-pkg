## Copyright (C) 2021(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{nodes}, @var{elem}] = fem_pre_mesh_extrude_surf(@var{mesh}, @var{elem_type}, @var{grp_id}, @var{h})
## Extrude a surface mesh along the surface normal in order to generate a volume mesh
##
## @var{mesh} @dots{} existing mesh with surface elements of type @var{elem_type}
##
## @var{elem_type} @dots{} element type of surface elements to be extruded
##
## @var{grp_id} @dots{} group id of surface elements to be extruded
##
## @var{h} @dots{} each element of the vector @var{h} defines the thickness of one extruded layer normal to the surface
##
## @end deftypefn

function [nodes, elem] = fem_pre_mesh_extrude_surf(mesh, elem_type_solid, grp_id, h)
  if (nargin ~= 4 || nargout ~= 2)
    print_usage();
  endif

  switch (elem_type_solid)
    case "iso8"
      elem_type = "iso4";
      corner_nodes = int32([1, 2, 3, 4]);
      num_nodes_elem = int32(8);
      bottom_node_idx = int32([5, 6, 7, 8]);
      top_node_idx = int32([1, 2, 3, 4]);
      interm_node_idx = [];
    case "penta15"
      elem_type = "tria6h";
      corner_nodes = int32([1, 2, 3]);
      num_nodes_elem = int32(15);
      bottom_node_idx = int32([1, 2, 3, 7, 8, 9]);
      top_node_idx = int32([4, 5, 6, 10, 11, 12]);
      interm_node_idx = int32([13, 14, 15]);
    case "iso20r"
      elem_type = "quad8r";
      corner_nodes = int32([1, 2, 3, 4]);
      num_nodes_elem = int32(20);
      bottom_node_idx = int32([1, 2, 3, 4, 9, 10, 11, 12]);
      top_node_idx = int32([5, 6, 7, 8, 13, 14, 15, 16]);
      interm_node_idx = int32([17, 18, 19, 20]);
    case "iso20"
      elem_type = "quad8";
      corner_nodes = int32([1, 2, 3, 4]);
      num_nodes_elem = int32(20);
      bottom_node_idx = int32([5, 6, 7, 8, 13, 14, 15, 16]);
      top_node_idx = int32([1, 2, 3, 4, 9, 10, 11, 12]);
      interm_node_idx = int32([17, 18, 19, 20]);
    case "iso27"
      elem_type = "quad9";
      corner_nodes = int32([1, 2, 3, 4, 5, 6, 7, 8, 9]);
      num_nodes_elem = int32(27);
      bottom_node_idx = int32([1, 2, 3, 4, 9, 10, 11, 12, 21]);
      top_node_idx = int32([5, 6, 7, 8, 17, 18, 19, 20, 26]);
      interm_node_idx = int32([13, 14, 15, 16, 22, 23, 24, 25, 27]);
    case "penta18"
      elem_type = "tria6h";
      corner_nodes = int32([1, 2, 3, 4, 5, 6]);
      num_nodes_elem = int32(18);
      bottom_node_idx = int32([1, 2, 3, 7, 8, 9]);
      top_node_idx = int32([4, 5, 6, 13, 14, 15]);
      interm_node_idx = int32([10, 11, 12, 16, 17, 18]);
    otherwise
      error("elem_type \"%s\" not supported", elem_type);
  endswitch
  
  layers = numel(h);

  elem_grp = getfield(mesh.groups, elem_type);

  grp_idx = find([elem_grp.id] == grp_id);

  if (isempty(grp_idx))
    error("group id %d not found in mesh", grp_id);
  endif

  elem_idx = [[elem_grp(grp_idx)].elements];

  elem_nodes = getfield(mesh.elements, elem_type)(elem_idx, :);

  load_case_dof_n.locked_dof = false(rows(mesh.nodes), 1);
  load_case_dof_n.domain = FEM_DO_ACOUSTICS;

  mesh_n.nodes = mesh.nodes;
  mesh_n.elements.particle_velocity = struct(elem_type, struct("nodes", elem_nodes));
  mesh_n.material_data.c = inf;
  mesh_n.material_data.rho = inf;
  mesh_n.materials.particle_velocity = struct(elem_type, ones(rows(elem_nodes), 1, "int32"));

  dof_map_n = fem_ass_dof_map(mesh_n, load_case_dof_n);

  load_case_n = struct();

  mat_ass.n = fem_ass_matrix(mesh_n, ...
                             dof_map_n, ...
                             [FEM_VEC_SURFACE_NORMAL_VECTOR], ...
                             load_case_n);

  top_nodes = unique(elem_nodes(:));
  interm_nodes = unique(elem_nodes(:, corner_nodes)(:));

  elem_n = getfield(mat_ass.n, elem_type);

  nodes = [mesh.nodes;
           zeros(layers * (numel(top_nodes) + numel(interm_nodes)), 6)];

  elem = zeros(layers * rows(elem_nodes), num_nodes_elem, "int32");

  elem(1:rows(elem_nodes), bottom_node_idx) = elem_nodes;

  for k=1:layers
    for i=1:rows(elem_nodes)
      for j=1:columns(elem_nodes)
        node_idx_ij = find(elem_nodes(i, j) == top_nodes) + rows(mesh.nodes) + (k - 1) * (numel(top_nodes) + numel(interm_nodes));
        nodes(node_idx_ij, 1:3) = mesh.nodes(elem_nodes(i, j), 1:3) + sum(h(1:k)) * reshape(elem_n(i, j, :), 1, 3);
        elem(i + (k - 1) * rows(elem_nodes), top_node_idx(j)) = node_idx_ij;
      endfor
      for j=1:numel(corner_nodes)
        node_idx_ij = find(elem_nodes(i, corner_nodes(j)) == interm_nodes) + rows(mesh.nodes) + numel(top_nodes)  + (k - 1) * (numel(top_nodes) + numel(interm_nodes));
        nodes(node_idx_ij, 1:3) = mesh.nodes(elem_nodes(i, corner_nodes(j)), 1:3) + (sum(h(1:k - 1)) + 0.5 * h(k)) * reshape(elem_n(i, corner_nodes(j), :), 1, 3);
        if (~isempty(interm_node_idx))
          elem(i + (k - 1) * rows(elem_nodes), interm_node_idx(j)) = node_idx_ij;
        endif
      endfor
    endfor
  endfor

  for k=2:layers
    elem((1:rows(elem_nodes)) + (k - 1) * rows(elem_nodes), bottom_node_idx) = elem((1:rows(elem_nodes)) + (k - 2) * rows(elem_nodes), top_node_idx);
  endfor
endfunction
