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
## @deftypefn {Function File} [@var{C1}, @var{C2}, @var{mesh_constr}, @var{dof_map_constr}, @var{mat_ass_constr}] = fem_ass_constraint_matrix(@var{mesh}, @var{elem_type_master}, @var{grp_id_master}, @var{elem_type_slave}, @var{grp_id_slave}, @var{maxdist}, @var{constrtype})
## Assemble a sparse matrix containing the constraint equations of a single group of surface to node constraints
## @end deftypefn

function [C1, C2, mesh_constr] = fem_ass_constraint_matrix(mesh, elem_type_master, grp_id_master, elem_type_slave, grp_id_slave, maxdist)
  grp_type_master = getfield(mesh.groups, elem_type_master);
  grp_type_slave = getfield(mesh.groups, elem_type_slave);

  grp_idx_master = find([grp_type_master.id] == grp_id_master);
  grp_idx_slave = find([grp_type_slave.id] == grp_id_slave);

  grp_master = getfield(mesh.groups, elem_type_master)(grp_idx_master);
  grp_slave = getfield(mesh.groups, elem_type_slave)(grp_idx_slave);

  elem_master = getfield(mesh.elements, elem_type_master)(grp_master.elements, :);
  elem_slave = getfield(mesh.elements, elem_type_master)(grp_slave.elements, :);

  node_idx_master = grp_master.nodes(:);
  node_idx_slave = grp_slave.nodes(:);

  node_idx_trans = zeros(rows(mesh.nodes), 1, "int32");

  node_idx_trans(node_idx_master) = 1:rows(node_idx_master);
  node_idx_trans(node_idx_slave) = rows(node_idx_master) + (1:rows(node_idx_slave));

  mesh_constr.nodes = mesh.nodes([node_idx_master; node_idx_slave], :);
  mesh_constr.elements = struct();

  sfncon.slave = node_idx_trans(node_idx_slave);
  sfncon.master = node_idx_trans(elem_master);
  sfncon.maxdist = maxdist;
  sfncon.constraint = FEM_CT_SLIDING;
  elem_type_map = {"tria6", "sfncon6";
                   "tria6h", "sfncon6h";
                   "iso4", "sfncon4";
                   "quad8", "sfncon8";
                   "quad9", "sfncon9";
                   "tria10", "sfncon10"};

  elem_type_sfncon = [];

  for i=1:rows(elem_type_map)
    switch (elem_type_master)
      case elem_type_map{i, 1}
        elem_type_sfncon = elem_type_map{i, 2};
        break
    endswitch
  endfor

  if (isempty(elem_type_sfncon))
    error("element type \"%s\" is not yet supported", elem_type_master);
  endif

  mesh_constr.elements = setfield(mesh_constr.elements, elem_type_sfncon, sfncon);

  mesh_constr.materials = struct();
  mesh_constr.material_data = struct()([]);

  dof_map_constr = [];
  mat_ass_constr = [];

  mesh_constr.joints = fem_pre_mesh_constr_surf_to_node(mesh_constr.nodes, mesh_constr.elements, FEM_DO_STRUCTURAL).joints;

  nnz_C2 = nnz_C1 = int32(0);

  for i=1:numel(mesh_constr.joints)
    nnz_C1 += (columns(mesh_constr.joints(i).C) - 6) / 2;
    nnz_C2 += 3;
  endfor

  ridx_C1 = zeros(nnz_C1, 1, "int32");
  ridx_C2 = zeros(nnz_C2, 1, "int32");
  cidx_C1 = zeros(nnz_C1, 1, "int32");
  cidx_C2 = zeros(nnz_C2, 1, "int32");
  data_C1 = zeros(nnz_C1, 1);
  data_C2 = zeros(nnz_C2, 1);

  idx_C2 = idx_C1 = int32(0);

  nnz_C2 = nnz_C1 = int32(0);

  for i=1:numel(mesh_constr.joints)
    for j=2:columns(mesh_constr.joints(i).nodes)
      ridx_C1(nnz_C1 + (1:3)) = i;
      cidx_C1(nnz_C1 + (1:3)) = (mesh_constr.joints(i).nodes(j) - 1) * 3 + (1:3);
      data_C1(nnz_C1 + (1:3)) = mesh_constr.joints(i).C(1, (j - 1) * 6 + (1:3));
      nnz_C1 += 3;
    endfor
    ridx_C2(nnz_C2 + (1:3)) = i;
    cidx_C2(nnz_C2 + (1:3)) = (mesh_constr.joints(i).nodes(1) - rows(node_idx_master) - 1) * 3 + (1:3);
    data_C2(nnz_C2 + (1:3)) = mesh_constr.joints(i).C(1, 1:3);
    nnz_C2 += 3;
  endfor

  C1 = sparse(ridx_C1, cidx_C1, data_C1, rows(node_idx_slave), rows(node_idx_master) * 3);
  C2 = sparse(ridx_C2, cidx_C2, data_C2, rows(node_idx_slave), rows(node_idx_slave) * 3);

  ## C1 * U1 + C2 * U2 = 0
  ## C2 * U2 = -C1 * U1
  ## U2n = -C1 * U1
endfunction
