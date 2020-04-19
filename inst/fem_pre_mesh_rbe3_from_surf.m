## Copyright Reinhard <octave-user@a1.net>
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

  empty_cell = cell(1, numel(group_id));
  
  rbe3 = struct("nodes", empty_cell, "weight", empty_cell, "elements", empty_cell);

  if (numel(group_id) ~= numel(master_node_idx))
    error("numel(group_id) does not match numel(master_node_idx)");
  endif
  
  group_idx = cell(1, numel(group_id));

  msh_groups = getfield(mesh.groups, elem_type);
  
  for i=1:numel(group_id)
    group_idx_i = find([msh_groups.id] == group_id(i));

    if numel(group_idx_i) == 0
      error("no elements found with id %d", group_id(i));
    endif

    group_idx{i} = group_idx_i;
  endfor
  
  for i=1:numel(group_idx)
    inode = [msh_groups(group_idx{i}).nodes];
    inode_tr = zeros(rows(mesh.nodes), 1, "int32");
    inode_tr(inode) = 1:numel(inode);
    press_elem.elements = inode_tr(getfield(mesh.elements, elem_type)([msh_groups(group_idx{i}).elements], :));
    press_elem.p = ones(rows(press_elem.elements), columns(press_elem.elements));
    load_case_i.pressure = setfield(struct(), elem_type, press_elem);
    dof_map_i.ndof = [reshape(1:(numel(inode) * 3), numel(inode), 3), zeros(numel(inode), 3)];
    dof_map_i.totdof = numel(inode) * 3;
    mesh_i.nodes = mesh.nodes(inode, :);
    mesh_i.elements = struct();
    mesh_i.materials = struct();
    mesh_i.material_data = struct("C", [], "rho", [])([]);
    mat_ass_i.R = fem_ass_matrix(mesh_i, ...
                                                 dof_map_i, ...
                                                 [FEM_VEC_LOAD_CONSISTENT], ...
                                                 load_case_i);

    F = norm(mat_ass_i.R(dof_map_i.ndof(:, 1:3)), "rows");
  
    rbe3(i).nodes = int32([master_node_idx(i), inode]);
    rbe3(i).weight = F.';
    rbe3(i).elements = setfield(struct(), elem_type, [msh_groups(group_idx{i}).elements]);
    
    if (any(rbe3(i).weight < 0))
      error("negative weight factor detected");
    endif
  endfor
endfunction
