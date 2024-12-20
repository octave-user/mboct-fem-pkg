## Copyright (C) 2019(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{offset}] = fem_pre_mesh_merge(@var{mesh_data})
## Merge several meshes into a single mesh.
##
## @var{mesh_data} @dots{} Struct array of mesh data including fields @var{mesh_data}.mesh
##
## @var{mesh} @dots{} Output mesh containing all nodes, elements and groups of all input meshes.
##
## @var{offset} @dots{} A struct containing the offset of node-, element- and material numbers of the combined mesh
## @end deftypefn

function [mesh, offset] = fem_pre_mesh_merge(mesh_data)
  ## merge nodes
  offset.node_idx = zeros(numel(mesh_data) + 1, 1, "int32");

  for i=1:numel(mesh_data)
    offset.node_idx(i + 1) = offset.node_idx(i) + rows(mesh_data(i).mesh.nodes);
  endfor

  mesh.nodes = zeros(offset.node_idx(end), 6);

  for i=1:numel(mesh_data)
    mesh.nodes(offset.node_idx(i) + (1:rows(mesh_data(i).mesh.nodes)), :) = mesh_data(i).mesh.nodes;
  endfor

  ## merge material data
  offset.material_data_idx = zeros(numel(mesh_data) + 1, 1, "int32");

  for i=1:numel(mesh_data)
    offset.material_data_idx(i + 1) = offset.material_data_idx(i) + numel(mesh_data(i).mesh.material_data);
  endfor

  empty_cell = cell(1, offset.material_data_idx(end));

  mesh.material_data = struct("E", empty_cell, ...
                              "nu", empty_cell, ...
                              "C", empty_cell, ...
                              "rho", empty_cell, ...
                              "alpha", empty_cell, ...
                              "beta", empty_cell, ...
                              "gamma", empty_cell,
                              "tan_delta", empty_cell, ...
                              "k", empty_cell,
                              "cp", empty_cell);

  for i=1:numel(mesh_data)
    mat_prop_names = fieldnames(mesh_data(i).mesh.material_data);
    for j=1:numel(mat_prop_names)
      for k=1:numel(mesh_data(i).mesh.material_data)
        mat_prop_ijk = getfield(mesh_data(i).mesh.material_data(k), mat_prop_names{j});
        mesh.material_data(offset.material_data_idx(i) + k) = setfield(mesh.material_data(offset.material_data_idx(i) + k), mat_prop_names{j}, mat_prop_ijk);
      endfor
    endfor
  endfor

  ## merge elements and materials
  elem_types = fem_pre_mesh_elem_type();

  mesh.elements = struct();
  mesh.materials = struct();
  mesh.groups = struct();
  offset.elements = struct();

  for j=1:numel(elem_types)
    offset_elem_idx_j = zeros(numel(mesh_data) + 1, 1, "int32");
    for i=1:numel(mesh_data)
      num_elem_i = 0;
      if (isfield(mesh_data(i).mesh.elements, elem_types(j).name))
        elem_nodes_i = getfield(mesh_data(i).mesh.elements, elem_types(j).name);
        num_elem_i = rows(elem_nodes_i);
      endif
      offset_elem_idx_j(i + 1) = offset_elem_idx_j(i) + num_elem_i;
    endfor

    if (offset_elem_idx_j(end) == 0)
      continue
    endif

    elem_nodes = zeros(offset_elem_idx_j(end), numel(elem_types(j).norder), "int32");
    elem_mat = zeros(offset_elem_idx_j(end), 1, "int32");

    for i=1:numel(mesh_data)
      if (isfield(mesh_data(i).mesh.elements, elem_types(j).name))
        elem_nodes_i = getfield(mesh_data(i).mesh.elements, elem_types(j).name);
        elem_nodes(offset_elem_idx_j(i) + (1:rows(elem_nodes_i)), :) = elem_nodes_i + offset.node_idx(i);

        if (isfield(mesh_data(i).mesh.materials, elem_types(j).name))
          elem_mat_i = getfield(mesh_data(i).mesh.materials, elem_types(j).name);
          elem_mat(offset_elem_idx_j(i) + (1:rows(elem_nodes_i))) = elem_mat_i + offset.material_data_idx(i);
        endif
      endif
    endfor

    offset_group_j = zeros(numel(mesh_data) + 1, 1, "int32");

    for i=1:numel(mesh_data)
      num_groups_i = 0;

      if (isfield(mesh_data(i).mesh, "groups") && isfield(mesh_data(i).mesh.groups, elem_types(j).name))
        num_groups_i = numel(getfield(mesh_data(i).mesh.groups, elem_types(j).name));
      endif
      offset_group_j(i + 1) = offset_group_j(i) + num_groups_i;
    endfor

    if (offset_group_j(end) > 0)
      empty_cell = cell(1, offset_group_j(end));

      mesh_groups_j = struct("id", empty_cell, "name", empty_cell, "nodes", empty_cell, "elements", empty_cell);

      for i=1:numel(mesh_data)
        if (isfield(mesh_data(i).mesh, "groups") && isfield(mesh_data(i).mesh.groups, elem_types(j).name))
          mesh_group_i = getfield(mesh_data(i).mesh.groups, elem_types(j).name);

          for k=1:numel(mesh_group_i)
            mesh_groups_j(offset_group_j(i) + k).id = mesh_group_i(k).id;
            mesh_groups_j(offset_group_j(i) + k).name = mesh_group_i(k).name;
            mesh_groups_j(offset_group_j(i) + k).nodes = mesh_group_i(k).nodes + offset.node_idx(i);
            mesh_groups_j(offset_group_j(i) + k).elements = mesh_group_i(k).elements + offset_elem_idx_j(i);
          endfor
        endif
      endfor

      mesh.groups = setfield(mesh.groups, elem_types(j).name, mesh_groups_j);
    endif

    mesh.elements = setfield(mesh.elements, elem_types(j).name, elem_nodes);
    mesh.materials = setfield(mesh.materials, elem_types(j).name, elem_mat);
    offset.elements = setfield(offset.elements, elem_types(j).name, offset_elem_idx_j);
  endfor


endfunction

%!test
%! mesh_data(1).mesh.nodes = repmat(1000, 3, 6);
%! mesh_data(1).mesh.elements.iso8 = repmat(int32(1000), 3, 8);
%! mesh_data(1).mesh.materials.iso8 = int32([1;2;2]);
%! mesh_data(1).mesh.elements.line2 = repmat(int32(1100), 2, 2);
%! mesh_data(1).mesh.material_data(1).E = 1000;
%! mesh_data(1).mesh.material_data(1).nu = 0.1;
%! mesh_data(1).mesh.material_data(1).rho = 10000;
%! mesh_data(1).mesh.material_data(2).E = 1100;
%! mesh_data(1).mesh.material_data(2).nu = 0.11;
%! mesh_data(1).mesh.material_data(2).rho = 11000;
%! mesh_data(1).mesh.groups.iso8(1).id = 1000;
%! mesh_data(1).mesh.groups.iso8(1).name = "iso8 1000";
%! mesh_data(1).mesh.groups.iso8(1).nodes = int32(1:3)(:);
%! mesh_data(1).mesh.groups.iso8(1).elements = int32(1:3)(:);
%! mesh_data(1).mesh.groups.iso8(2).id = 1100;
%! mesh_data(1).mesh.groups.iso8(2).name = "iso8 1100";
%! mesh_data(1).mesh.groups.iso8(2).nodes = int32(1:2)(:);
%! mesh_data(1).mesh.groups.iso8(2).elements = int32(1:2)(:);

%! mesh_data(2).mesh.nodes = repmat(2000, 2, 6);
%! mesh_data(2).mesh.elements.iso8 = repmat(int32(2000), 2, 8);
%! mesh_data(2).mesh.materials.iso8 = int32([1;2]);
%! mesh_data(2).mesh.elements.line2 = repmat(int32(2100), 4, 2);
%! mesh_data(2).mesh.material_data(1).E = 2000;
%! mesh_data(2).mesh.material_data(1).nu = 0.2;
%! mesh_data(2).mesh.material_data(1).rho = 20000;
%! mesh_data(2).mesh.material_data(2).E = 2200;
%! mesh_data(2).mesh.material_data(2).nu = 0.22;
%! mesh_data(2).mesh.material_data(2).rho = 22000;

%! mesh_data(2).mesh.groups.iso8(1).id = 2000;
%! mesh_data(2).mesh.groups.iso8(1).name = "iso8 2000";
%! mesh_data(2).mesh.groups.iso8(1).nodes = int32(1:2)(:);
%! mesh_data(2).mesh.groups.iso8(1).elements = int32(1:2)(:);
%! mesh_data(2).mesh.groups.iso8(2).id = 2100;
%! mesh_data(2).mesh.groups.iso8(2).name = "iso8 2100";
%! mesh_data(2).mesh.groups.iso8(2).nodes = int32(1)(:);
%! mesh_data(2).mesh.groups.iso8(2).elements = int32(1)(:);

%! mesh_data(3).mesh.nodes = repmat(3000, 4, 6);
%! mesh_data(3).mesh.elements.iso8 = repmat(int32(3000), 2, 8);
%! mesh_data(3).mesh.materials.iso8 = int32([1;1]);
%! mesh_data(3).mesh.elements.line2 = repmat(int32(3100), 5, 2);
%! mesh_data(3).mesh.material_data.E = 3000;
%! mesh_data(3).mesh.material_data.nu = 0.3;
%! mesh_data(3).mesh.material_data.rho = 30000;
%! mesh_data(3).mesh.groups.iso8(1).id = 3000;
%! mesh_data(3).mesh.groups.iso8(1).name = "iso8 3000";
%! mesh_data(3).mesh.groups.iso8(1).nodes = int32(1:2)(:);
%! mesh_data(3).mesh.groups.iso8(1).elements = int32(1:2)(:);

%! mesh_data(4).mesh.nodes = repmat(4000, 4, 6);
%! mesh_data(4).mesh.elements.iso20 = repmat(int32(4000), 2, 20);
%! mesh_data(4).mesh.materials.iso20 = int32([1;1]);
%! mesh_data(4).mesh.elements.line2 = repmat(int32(4100), 5, 2);
%! mesh_data(4).mesh.material_data.E = 4000;
%! mesh_data(4).mesh.material_data.nu = 0.3;
%! mesh_data(4).mesh.material_data.rho = 40000;
%! mesh_data(4).mesh.groups.iso20(1).id = 4000;
%! mesh_data(4).mesh.groups.iso20(1).name = "iso8 4000";
%! mesh_data(4).mesh.groups.iso20(1).nodes = int32(1:2)(:);
%! mesh_data(4).mesh.groups.iso20(1).elements = int32(1:2)(:);

%! [mesh, offset] = fem_pre_mesh_merge(mesh_data);

%! grp_idx_1000 = find([mesh.groups.iso8.id] == 1000);
%! assert(all(mesh.elements.iso8(mesh.groups.iso8(grp_idx_1000).elements) == 1000));
%! grp_idx_2000 = find([mesh.groups.iso8.id] == 2000);
%! assert(all(mesh.elements.iso8(mesh.groups.iso8(grp_idx_2000).elements) == 2000 + rows(mesh_data(1).mesh.nodes)));
%! grp_idx_3000 = find([mesh.groups.iso8.id] == 3000);
%! assert(all(mesh.elements.iso8(mesh.groups.iso8(grp_idx_3000).elements) == 3000 + rows(mesh_data(1).mesh.nodes) + rows(mesh_data(2).mesh.nodes)));
%! grp_idx_4000 = find([mesh.groups.iso20.id] == 4000);
%! assert(all(mesh.elements.iso20(mesh.groups.iso8(grp_idx_4000).elements) == 4000 + rows(mesh_data(1).mesh.nodes) + rows(mesh_data(2).mesh.nodes) + rows(mesh_data(3).mesh.nodes)));
