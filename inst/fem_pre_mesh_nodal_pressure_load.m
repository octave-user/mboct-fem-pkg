## Copyright (C) 2023(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{F}, @var{group_idx}, @var{weight}] = fem_pre_mesh_nodal_pressure_load(@var{mesh}, @var{group_id}, @var{elem_type})
##
## Compute unit pressure loads for surface elements
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_id} @dots{} array of group id-numbers used to identify surface elements with slave nodes
##
## @var{elem_type} @dots{} the element type addressed by @var{group_id} (e.g. "tria6", "tria6h", "iso4", "quad8")
##
## @var{F} @dots{} resulting nodal load
##
## @var{group_idx} @dots{} index of all element groups identified by @var{group_id}
##
## @var{weight} @dots{} resulting weight factors based on the affected surface area
##
## @end deftypefn

function [F, group_idx, weight, node_id] = fem_pre_mesh_nodal_pressure_load(mesh, group_id, elem_type)
  if (nargin < 3 || nargout > 4)
    print_usage();
  endif

  if (~iscell(elem_type))
    elem_type = {elem_type};
  endif

  F = weight = node_id = cell(1, numel(group_id));
  group_idx = cell(numel(group_id), numel(elem_type));

  for i=1:numel(group_id)
    inum_nodes = int32(0);

    for j=1:numel(elem_type)
      if (~isfield(mesh.groups, elem_type{j}))
        continue;
      endif
      msh_groups = getfield(mesh.groups, elem_type{j});

      group_idx_ij = find([msh_groups.id] == group_id(i));

      if (isempty(group_idx_ij))
        continue;
      endif

      group_idx{i, j} = group_idx_ij;
      inum_nodes += numel([[msh_groups(group_idx{i, j})].nodes]);

      clear msh_groups group_idx_ij;
    endfor

    node_id{i} = zeros(1, inum_nodes, "int32");
    inum_nodes = int32(0);

    for j=1:numel(elem_type)
      if (~isfield(mesh.groups, elem_type{j}))
        continue;
      endif

      msh_groups = getfield(mesh.groups, elem_type{j});

      if (isempty(group_idx{i, j}))
        continue;
      endif

      nodes_j = [[msh_groups(group_idx{i, j})].nodes];

      node_id{i}((1:numel(nodes_j)) + inum_nodes) = nodes_j;
      inum_nodes += numel(nodes_j);

      clear msh_groups nodes_j;
    endfor

    node_id{i} = unique(node_id{i});

    clear inum_nodes;
  endfor

  for i=1:rows(group_idx)
    F{i} = zeros(numel(node_id{i}), 3);

    if (nargout >= 3)
      weight{i} = zeros(1, numel(node_id{i}));
    endif

    for j=1:columns(group_idx)
      inode = node_id{i};
      inode_tr = zeros(1, rows(mesh.nodes), "int32");
      inode_tr(inode) = 1:numel(inode);

      if (~isfield(mesh.groups, elem_type{j}))
        continue;
      endif

      msh_groups = getfield(mesh.groups, elem_type{j});
      ielno = getfield(mesh.elements, elem_type{j})([msh_groups(group_idx{i, j}).elements], :);
      press_elem.elements = inode_tr(ielno);
      press_elem.p = repmat(1, rows(press_elem.elements), columns(press_elem.elements));
      load_case_i.pressure = setfield(struct(), elem_type{j}, press_elem);
      dof_map_i.ndof = [reshape(1:(numel(inode) * 3), numel(inode), 3), zeros(numel(inode), 3)];
      dof_map_i.totdof = numel(inode) * 3;
      dof_map_i.domain = FEM_DO_STRUCTURAL;
      mesh_i.nodes = mesh.nodes(inode, :);
      mesh_i.elements = struct();
      mesh_i.materials = struct();
      mesh_i.material_data = struct("E", cell(0, 0), "nu", cell(0, 0), "C", cell(0, 0), "rho", cell(0, 0));
      [mat_ass_i.R, ...
       mat_ass_i.surface] = fem_ass_matrix(mesh_i, ...
                                           dof_map_i, ...
                                           [FEM_VEC_LOAD_CONSISTENT, ...
                                            FEM_VEC_SURFACE_AREA], ...
                                           load_case_i);
      if (nargout >= 3)
        elem_surf = getfield(mat_ass_i.surface, elem_type{j});

        for j=1:columns(press_elem.elements)
          for k=1:rows(press_elem.elements)
            weight{i}(press_elem.elements(k, j)) += elem_surf(k, j); ## Cannot be vectorized because we can have duplicate nodes
          endfor
        endfor
      endif

      F{i} += full(mat_ass_i.R(dof_map_i.ndof(:, 1:3)));

      clear inode inode_tr msh_groups ielno press_elem elem_surf load_case_i mesh_i dof_map_i mat_ass_i;
    endfor
  endfor
endfunction
