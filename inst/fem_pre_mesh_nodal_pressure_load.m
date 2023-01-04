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
## @deftypefn {Function File} [@var{F}, @var{group_idx}] = fem_pre_mesh_nodal_pressure_load(@var{mesh}, @var{group_id}, @var{elem_type})
##
## Compute unit pressure loads for surface elements
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_id} @dots{} array of group id-numbers used to identify surface elements with slave nodes
##
## @var{elem_type} @dots{} the element type addressed by @var{group_id} (e.g. "tria6", "tria6h", "iso4", "quad8")
##
## @end deftypefn

function [F, group_idx] = fem_pre_mesh_nodal_pressure_load(mesh, group_id, elem_type)
  if (nargin < 3 || nargout > 2)
    print_usage();
  endif

  F = group_idx = cell(1, numel(group_id));

  msh_groups = getfield(mesh.groups, elem_type);

  for i=1:numel(group_id)
    group_idx_i = find([msh_groups.id] == group_id(i));

    if (isempty(group_idx_i))
      error("no elements found with id %d", group_id(i));
    endif

    group_idx{i} = group_idx_i;
  endfor

  for i=1:numel(group_idx)
    inode = [msh_groups(group_idx{i}).nodes];
    inode_tr = zeros(1, rows(mesh.nodes), "int32");
    inode_tr(inode) = 1:numel(inode);
    ielno = getfield(mesh.elements, elem_type)([msh_groups(group_idx{i}).elements], :);
    press_elem.elements = inode_tr(ielno);
    press_elem.p = repmat(1, rows(press_elem.elements), columns(press_elem.elements));
    load_case_i.pressure = setfield(struct(), elem_type, press_elem);
    dof_map_i.ndof = [reshape(1:(numel(inode) * 3), numel(inode), 3), zeros(numel(inode), 3)];
    dof_map_i.totdof = numel(inode) * 3;
    dof_map_i.domain = FEM_DO_STRUCTURAL;
    mesh_i.nodes = mesh.nodes(inode, :);
    mesh_i.elements = struct();
    mesh_i.materials = struct();
    mesh_i.material_data = struct("E", cell(0, 0), "nu", cell(0, 0), "C", cell(0, 0), "rho", cell(0, 0));
    mat_ass_i.R = fem_ass_matrix(mesh_i, ...
                                 dof_map_i, ...
                                 [FEM_VEC_LOAD_CONSISTENT], ...
                                 load_case_i);

    F{i} = full(mat_ass_i.R(dof_map_i.ndof(:, 1:3)));
  endfor
endfunction
