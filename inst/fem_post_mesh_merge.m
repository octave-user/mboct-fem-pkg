## Copyright (C) 2019(-2020) Reinhard <octave-user@a1.net>
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

  elem_types = struct("name", cell(1, 6), "have_mat", cell(1, 6));
  
  elem_types(1).name = "iso8";
  elem_types(1).have_mat = true;
  elem_types(2).name = "tet10";
  elem_types(2).have_mat = true;
  elem_types(3).name = "iso4";
  elem_types(3).have_mat = false;
  elem_types(4).name = "tria6";
  elem_types(4).have_mat = false;
  elem_types(5).name = "tria3";
  elem_types(5).have_mat = false;
  elem_types(6).name = "iso20";
  elem_types(6).have_mat = true;
  elem_types(7).name = "quad8";
  elem_types(7).have_mat = false;
  
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
  mesh.material_data = repmat(struct("C", nan(6, 6), "rho", -1), 1, num_mat);

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

  for i=1:numel(mesh_data)
    for j=1:numel(mesh_data(i).mesh.material_data)
      mesh.material_data(dof_map.submesh.offset.materials(i) + j).C = mesh_data(i).mesh.material_data(j).C;
      mesh.material_data(dof_map.submesh.offset.materials(i) + j).rho = mesh_data(i).mesh.material_data(j).rho;
    endfor
  endfor

  dof_map.submesh.offset.elements = struct();
  
  for i=1:numel(elem_types)
    elem_types(i).num_elem = int32(0);
    elem_types(i).num_nodes = int32(0);
    
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
          offset_elem(i + 1) = offset_elem(i) + rows(curr_elem);
        endif
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

%!test ##demo
%! close all;
%! c = 10e-3;
%!
%! material1.E = 210000e6;
%! material1.nu = 0.3;
%! material1.rho = 7850;
%!
%! material2.E = material1.E;
%! material2.nu = material1.nu;
%! material2.rho = material1.rho;
%!
%! h1 = 4.5e-3;
%! h2 = 5.5e-3;
%!
%! geometry1.l = 20e-3;
%! geometry1.w = 20e-3;
%! geometry1.h = 40e-3;
%! geometry2.l = 60e-3;
%! geometry2.w = 30e-3;
%! geometry2.h = 50e-3;
%!
%! mesh_size1.num_elem_l = ceil(geometry1.l / h1);
%! mesh_size1.num_elem_w = ceil(geometry1.w / h1);
%! mesh_size1.num_elem_h = ceil(geometry1.h / h1);
%! mesh_size2.num_elem_l = ceil(geometry2.l / h2);
%! mesh_size2.num_elem_w = ceil(geometry2.w / h2);
%! mesh_size2.num_elem_h = ceil(geometry2.h / h2);
%!
%! [data(1).mesh] = fem_pre_mesh_cube_create(geometry1, mesh_size1, material1, zeros(3, 1));
%!
%! data(1).mesh.nodes(:, 2) -= 0.5 * geometry1.w;
%! data(1).mesh.nodes(:, 3) -= 0.5 * geometry1.h;
%!
%! [data(2).mesh] = fem_pre_mesh_cube_create(geometry2, mesh_size2, material2, zeros(3, 1));
%!
%! data(2).mesh.nodes(:, 1) += geometry1.l + c;
%! data(2).mesh.nodes(:, 2) -= 0.5 * geometry2.w;
%! data(2).mesh.nodes(:, 3) -= 0.5 * geometry2.h;
%!
%! for i=1:numel(data)
%!  data(i).load_case.locked_dof = false(size(data(i).mesh.nodes));
%!  data(i).dof_map = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%! endfor
%!
%! [mesh, dof_map] = fem_post_mesh_merge(data);

%! figure("visible", "off");
%! fem_post_sol_plot(mesh);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("merged mesh");
%! view(30, 30);
%! figure_list();
