## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{dof_map}] = fem_post_mesh_merge(@var{mesh_data}, @var{options})
## Merge several meshes into a single mesh.
##
## @var{mesh_data} @dots{} Struct array of mesh data including fields @var{mesh_data}.mesh and @var{mesh_data}.dof_map
##
## @var{options}.group_id @dots{} If @var{options}.group_id == "preserve", then the output mesh will contain the same group_id's like the input meshes. If @var{options}.group_id == "unique", then a new unique group_id will be assigned to each group.
##
## @var{mesh} @dots{} Output mesh containing all nodes, elements and groups of all input meshes.
##
## @var{dof_map} @dots{} Degree of freedom mapping for the output mesh
## @end deftypefn

function [mesh, dof_map] = fem_post_mesh_merge(mesh_data, options)
  if (nargin < 1 || nargin > 2 || nargout < 1 || nargout > 2)
    print_usage();
  endif

  if (nargout >= 2)
    dof_map.totdof = int32(0);
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "group_id"))
    options.group_id = "preserve";
  endif

  num_nodes = int32(0);
  num_mat = int32(0);

  eltype = fem_pre_mesh_elem_type();

  elem_types = struct("name", {eltype.name}, ...
                      "num_elem", mat2cell(zeros(1, numel(eltype), "int32"), 1, ones(1, numel(eltype), "int32")), ...
                      "num_nodes", mat2cell(zeros(1, numel(eltype), "int32"), 1, ones(1, numel(eltype), "int32")), ...
                      "have_mat", mat2cell([eltype.dim] == 3, 1, ones(1, numel(eltype), "int32")));

  for i=1:numel(mesh_data)
    if (nargout >= 2)
      dof_map.totdof += mesh_data(i).dof_map.totdof;
    endif
    num_nodes += rows(mesh_data(i).mesh.nodes);
    num_mat += numel(mesh_data(i).mesh.material_data);
  endfor

  mesh.nodes = zeros(num_nodes, 6);
  mesh.elements = struct();
  mesh.materials = struct();
  empty_cell = cell(1, num_mat);
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

  if (nargout >= 2)
    dof_map.ndof = zeros(num_nodes, 6, "int32");
    dof_map.idx_lambda = [];
    dof_map.submesh.offset.dof = zeros(numel(mesh_data), 1, "int32");
  endif

  dof_map.submesh.offset.nodes = zeros(numel(mesh_data), 1, "int32");
  dof_map.submesh.offset.materials = zeros(numel(mesh_data), 1, "int32");

  for i=1:numel(mesh_data) - 1
    if (nargout >= 2)
      dof_map.submesh.offset.dof(i + 1) = dof_map.submesh.offset.dof(i) + mesh_data(i).dof_map.totdof;
    endif
    dof_map.submesh.offset.nodes(i + 1) = dof_map.submesh.offset.nodes(i) + rows(mesh_data(i).mesh.nodes);
    dof_map.submesh.offset.materials(i + 1) = dof_map.submesh.offset.materials(i) + numel(mesh_data(i).mesh.material_data);
  endfor

  for i=1:numel(mesh_data)
    mesh.nodes((1:rows(mesh_data(i).mesh.nodes)) + dof_map.submesh.offset.nodes(i), :) = mesh_data(i).mesh.nodes;
    if (nargout >= 2)
      dof_map.ndof((1:rows(mesh_data(i).dof_map.ndof)) + dof_map.submesh.offset.nodes(i), :) = mesh_data(i).dof_map.ndof ...
                                                                                               + dof_map.submesh.offset.dof(i) ...
                                                                                                 * (mesh_data(i).dof_map.ndof > 0);
    endif
  endfor

  mat_prop = fieldnames(mesh.material_data);

  for i=1:numel(mesh_data)
    for j=1:numel(mesh_data(i).mesh.material_data)
      for k=1:numel(mat_prop)
        if (isfield(mesh_data(i).mesh.material_data(j), mat_prop{k}))
          mesh.material_data(dof_map.submesh.offset.materials(i) + j) = setfield(mesh.material_data(dof_map.submesh.offset.materials(i) + j), ...
                                                                                 mat_prop{k}, ...
                                                                                 getfield(mesh_data(i).mesh.material_data(j), mat_prop{k}));
        endif
      endfor
    endfor
  endfor

  dof_map.submesh.offset.elements = struct();

  for i=1:numel(elem_types)
    for j=1:numel(mesh_data)
      if (isfield(mesh_data(j).mesh.elements, elem_types(i).name))
        elem_nodes = getfield(mesh_data(j).mesh.elements, elem_types(i).name);
        elem_types(i).num_elem += rows(elem_nodes);
        elem_types(i).num_nodes = columns(elem_nodes);
      endif
    endfor

    elem_offset = zeros(numel(mesh_data), 1, "int32");

    for j=1:numel(mesh_data) - 1
      num_elem = int32(0);

      if (isfield(mesh_data(j).mesh.elements, elem_types(i).name))
        num_elem = rows(getfield(mesh_data(j).mesh.elements, elem_types(i).name));
      endif

      elem_offset(j + 1) = elem_offset(j) + num_elem;
    endfor

    dof_map.submesh.offset.elements = setfield(dof_map.submesh.offset.elements, elem_types(i).name, elem_offset);
  endfor

  for i=1:numel(elem_types)
    if (elem_types(i).num_elem)
      elem_nodes_m = zeros(elem_types(i).num_elem, elem_types(i).num_nodes, "int32");

      if (elem_types(i).have_mat)
        elem_mat_m = zeros(elem_types(i).num_elem, 1, "int32");
      endif

      curr_elem_m = int32(0);

      for j=1:numel(mesh_data)
        if (isfield(mesh_data(j).mesh.elements, elem_types(i).name))
          elem_nodes = getfield(mesh_data(j).mesh.elements, elem_types(i).name);
          elem_nodes_m(curr_elem_m + (1:rows(elem_nodes)), :) = elem_nodes + dof_map.submesh.offset.nodes(j);

          if (elem_types(i).have_mat)
            elem_mat = getfield(mesh_data(j).mesh.materials, elem_types(i).name);
            elem_mat_m(curr_elem_m + (1:rows(elem_mat)), :) = elem_mat + dof_map.submesh.offset.materials(j);
          endif

          curr_elem_m += rows(elem_nodes);
        endif
      endfor

      mesh.elements = setfield(mesh.elements, elem_types(i).name, elem_nodes_m);

      if (elem_types(i).have_mat)
        mesh.materials = setfield(mesh.materials, elem_types(i).name, elem_mat_m);
      endif
    endif
  endfor

  num_rbe3 = int32(0);

  for i=1:numel(mesh_data)
    if (isfield(mesh_data(i).mesh.elements, "rbe3"))
      num_rbe3 += numel(mesh_data(i).mesh.elements.rbe3);
    endif
  endfor

  if (num_rbe3)
    mesh.elements.rbe3 = repmat(struct("nodes", [], "weight", []), 1, num_rbe3);
    dof_map.edof.rbe3 = zeros(num_rbe3, 6, "int32");
    curr_rbe3 = int32(0);

    for i=1:numel(mesh_data)
      if (isfield(mesh_data(i).mesh.elements, "rbe3"))
        for j=1:numel(mesh_data(i).mesh.elements.rbe3)
          ++curr_rbe3;
          mesh.elements.rbe3(curr_rbe3).nodes = mesh_data(i).mesh.elements.rbe3(j).nodes + dof_map.submesh.offset.nodes(i);
          mesh.elements.rbe3(curr_rbe3).weight = mesh_data(i).mesh.elements.rbe3(j).weight;
          if (nargout >= 2)
            dof_map.edof.rbe3(curr_rbe3, :) = mesh_data(i).dof_map.edof.rbe3(j, :) + dof_map.submesh.offset.dof(i);
          endif
        endfor
      endif
    endfor

    if (nargout >= 2)
      dof_map.idx_lambda = [dof_map.idx_lambda(:), dof_map.edof.rbe3(:)];
    endif
  endif

  num_joints = int32(0);
  max_joint_dof = int32(0);

  for i=1:numel(mesh_data)
    if (isfield(mesh_data(i).mesh.elements, "joints"))
      num_joints += numel(mesh_data(i).mesh.elements.joints);

      for j=1:numel(mesh_data(i).mesh.elements.joints)
        max_joint_dof = max(max_joint_dof, columns(mesh_data(i).dof_map.edof.joints));
      endfor
    endif
  endfor

  if (num_joints)
    mesh.elements.joints = repmat(struct("nodes", [], "C", []), 1, num_joints);
    dof_map.edof.joints = zeros(num_joints, max_joint_dof, "int32");
    curr_joint = int32(0);
    curr_joint_dof = int32(0);

    for i=1:numel(mesh_data)
      if (isfield(mesh_data(i).mesh.elements, "joints"))
        for j=1:numel(mesh_data(i).mesh.elements.joints)
          ++curr_joint;
          mesh.elements.joints(curr_joint).nodes = mesh_data(i).mesh.elements.joints(j).nodes + dof_map.submesh.offset.nodes(i);
          mesh.elements.joints(curr_joint).C = mesh_data(i).mesh.elements.joints(j).C;
        endfor

        if (nargout >= 2)
          dof_map.edof.joints(curr_joint_dof + (1:rows(mesh_data(i).dof_map.edof.joints)), 1:columns(mesh_data(i).dof_map.edof.joints)) = mesh_data(i).dof_map.edof.joints + (mesh_data(i).dof_map.edof.joints > 0) * dof_map.submesh.offset.dof(i);
          curr_joint_dof += rows(mesh_data(i).dof_map.edof.joints);
        endif
      endif
    endfor

    if (nargout >= 2)
      dof_map.idx_lambda = [dof_map.idx_lambda(:); dof_map.edof.joints(:)];
    endif
  endif

  for j=1:numel(elem_types)
    num_groups = int32(0);
    offset_elem = zeros(1, numel(mesh_data) + 1, "int32");

    for i=1:numel(mesh_data)
      if (isfield(mesh_data(i).mesh, "groups") && isfield(mesh_data(i).mesh.groups, elem_types(j).name))
        ++num_groups;
      endif
    endfor

    if (num_groups)
      curr_group_m = repmat(struct("id", [], "name", [], "nodes", [], "elements", []), 1, num_groups);
      curr_group_idx_m = int32(0);

      for i=1:numel(mesh_data)
        if (isfield(mesh_data(i).mesh, "groups") && isfield(mesh_data(i).mesh.groups, elem_types(j).name))
          curr_group = getfield(mesh_data(i).mesh.groups, elem_types(j).name);
          curr_elem = getfield(mesh_data(i).mesh.elements, elem_types(j).name);
          for k=1:numel(curr_group)
            ++curr_group_idx_m;
            switch (options.group_id)
              case "unique"
                curr_group_m(curr_group_idx_m).id = curr_group_idx_m;
              otherwise
                curr_group_m(curr_group_idx_m).id = curr_group(k).id;
            endswitch
            curr_group_m(curr_group_idx_m).name = curr_group(k).name;
            curr_group_m(curr_group_idx_m).nodes = curr_group(k).nodes + dof_map.submesh.offset.nodes(i);
            curr_group_m(curr_group_idx_m).elements = curr_group(k).elements + offset_elem(i);
          endfor
          num_elem = rows(curr_elem);
        else
          num_elem = 0;
        endif
        offset_elem(i + 1) = offset_elem(i) + num_elem;
      endfor

      if (~isfield(mesh, "groups"))
        mesh.groups = struct();
      endif

      mesh.groups = setfield(mesh.groups, elem_types(j).name, curr_group_m);
    endif
  endfor

  if (nargout >= 2)
    dof_map.idx_lambda = sort(dof_map.idx_lambda(find(dof_map.idx_lambda > 0)));

    dof_map.idx_node = dof_map.ndof(:);

    dof_map.idx_node = sort(dof_map.idx_node(find(dof_map.idx_node > 0)));
  endif
endfunction

%!test
%! mesh_data(1).mesh.nodes = ones(15, 6) * 1000;
%! mesh_data(1).mesh.elements.penta15 = int32([1:15]);
%! mesh_data(1).mesh.materials.penta15 = int32(1);
%! mesh_data(1).mesh.material_data = struct("E", 100000e6, "nu", 0.3, "rho", 7850);
%! mesh_data(1).mesh.groups.penta15.id = int32(1000);
%! mesh_data(1).mesh.groups.penta15.nodes = unique(mesh_data(1).mesh.elements.penta15);
%! mesh_data(1).mesh.groups.penta15.elements = int32(1);
%! mesh_data(1).mesh.groups.penta15.name = "penta15 group 1";
%! mesh_data(2).mesh.nodes = ones(20+15, 6) * 2000;
%! mesh_data(2).mesh.elements.iso20 = int32([1:20]);
%! mesh_data(2).mesh.elements.penta15 = int32([1:15] + 20);
%! mesh_data(2).mesh.materials.penta15 = int32(1);
%! mesh_data(2).mesh.materials.iso20 = int32(1);
%! mesh_data(2).mesh.material_data = struct("E", 200000e6, "nu", 0.3, "rho", 7850);

%! mesh_data(2).mesh.groups.penta15.id = int32(2000);
%! mesh_data(2).mesh.groups.penta15.nodes = unique(mesh_data(2).mesh.elements.penta15);
%! mesh_data(2).mesh.groups.penta15.elements = int32(1);
%! mesh_data(2).mesh.groups.penta15.name = "penta15 group 2";

%! mesh_data(2).mesh.groups.iso20.id = int32(2100);
%! mesh_data(2).mesh.groups.iso20.nodes = unique(mesh_data(2).mesh.elements.iso20);
%! mesh_data(2).mesh.groups.iso20.elements = int32(1);
%! mesh_data(2).mesh.groups.iso20.name = "iso20 group 2";

%! mesh_data(3).mesh.nodes = ones(20 + 15, 6) * 3000;
%! mesh_data(3).mesh.elements.iso20 = int32([1:20]);
%! mesh_data(3).mesh.elements.penta15 = int32([1:15] + 20);
%! mesh_data(3).mesh.materials.penta15 = int32(1);
%! mesh_data(3).mesh.materials.iso20 = int32(1);
%! mesh_data(3).mesh.material_data = struct("E", 300000e6, "nu", 0.3, "rho", 7850);

%! mesh_data(3).mesh.groups.penta15.id = int32(3000);
%! mesh_data(3).mesh.groups.penta15.nodes = unique(mesh_data(3).mesh.elements.penta15);
%! mesh_data(3).mesh.groups.penta15.elements = int32(1);
%! mesh_data(3).mesh.groups.penta15.name = "penta15 group 3";

%! mesh_data(3).mesh.groups.iso20.id = int32(3100);
%! mesh_data(3).mesh.groups.iso20.nodes = unique(mesh_data(3).mesh.elements.iso20);
%! mesh_data(3).mesh.groups.iso20.elements = int32(1);
%! mesh_data(3).mesh.groups.iso20.name = "iso20 group 3";

%! mesh_data(4).mesh.nodes = ones(10, 6) * 4000;
%! mesh_data(4).mesh.material_data = struct()([]);
%! mesh_data(4).mesh.elements.point1 = int32(1:10)(:);
%! for i=1:10
%!   mesh_data(4).mesh.groups.point1(i).id = 4000 + i;
%!   mesh_data(4).mesh.groups.point1(i).nodes = i;
%!   mesh_data(4).mesh.groups.point1(i).elements = i;
%!   mesh_data(4).mesh.groups.point1(i).name = sprintf("point1 %d", i);
%! endfor

%! mesh1 = fem_post_mesh_merge(mesh_data);
%! mesh2 = fem_post_mesh_merge(mesh_data([4, 3, 2, 1]));
%! mesh3 = fem_post_mesh_merge(mesh_data([2, 1, 4, 3]));
%! mesh4 = fem_post_mesh_merge(mesh_data([4, 3, 2, 1]));
%! mesh5 = fem_post_mesh_merge(mesh_data([1, 3, 4, 2]));
%! mesh6 = fem_post_mesh_merge(mesh_data([2, 3, 4, 1]));
%! mesh7 = fem_post_mesh_merge(mesh_data([3, 1, 2, 4]));
