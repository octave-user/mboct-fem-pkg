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
## @deftypefn {Function File} @var{groups} = fem_pre_mesh_groups_create(@var{mesh}, @var{group_defs}, @var{tol_rel})
## @deftypefnx {} @dots{} = fem_pre_mesh_groups_create(@dots{}, @var{tol_abs})
## Create groups of elements or nodes which are located within geometrical limits defined by simple shapes like boxes or cylinders
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_defs} @dots{} Struct array of geometrical shapes
##
## @var{group_defs}.id @dots{} id-number of the group to be created
##
## @var{group_defs}.name @dots{} name of the group to be created
##
## @var{group_defs}.R @dots{} rotation matrix of geometrical shapes (e.g. R(:, 3) is equal to the axis of a cylinder)
##
## @var{group_defs}.X0 @dots{} origin of geometrical shapes (e.g. a point at the axis of a cylinder)
##
## @var{group_defs}.type @dots{} selected shape type (e.g. "cylinder", "box")
##
## @var{group_defs}.geometry.xmin @dots{} minimum x-distance for shape type "box"
##
## @var{group_defs}.geometry.xmax @dots{} maximum x-distance for shape type "box"
##
## @var{group_defs}.geometry.ymin @dots{} minimum y-distance for shape type "box"
##
## @var{group_defs}.geometry.ymax @dots{} maximum y-distance for shape type "box"
##
## @var{group_defs}.geometry.zmin @dots{} minimum z-distance for shape type "box" or "cylinder"
##
## @var{group_defs}.geometry.zmax @dots{} maximum z-distance for shape type "box" or "cylinder"
##
## @var{group_defs}.geometry.rmin @dots{} minimum radius for shape type "cylinder"
##
## @var{group_defs}.geometry.rmax @dots{} maximum radius for shape type "cylinder"
##
## @var{tol_rel} @dots{} relative tolerance for node positions
##
## @var{tol_abs} @dots{} absolute tolerance for node positions
##
## @end deftypefn

function groups = fem_pre_mesh_groups_create(mesh, group_defs, tolrel, tolabs, elem_type)
  if (nargin < 3 || nargin > 5 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    tolabs = 0;
  endif

  if (nargin < 5)
    ## default global element type
    elem_type = "tria6";
  endif
  
  dX = zeros(3, 1);

  for i=1:3
    dX(i) = max(mesh.nodes(:, i)) - min(mesh.nodes(:, i));
  endfor

  tol = norm(dX) * tolrel + tolabs;

  elem_types = {"iso4", "tria6"};

  empty_cell = cell(1, 0);

  groups.nodes = struct("id", empty_cell, "name", empty_cell, "nodes", empty_cell);

  for i=1:numel(elem_types)
    groups = setfield(groups, elem_types{i}, struct("id", empty_cell, "name", empty_cell, "elements", empty_cell, "nodes", empty_cell));
  endfor

  for j=1:numel(group_defs)
    group_defs_num_nodes = int32(0);

    if (isfield(group_defs, "elem_type"))
      ## overwrite global element type
      elem_type = group_defs(j).elem_type;
    endif
    
    switch (elem_type)
      case {"iso4", "tria6"}
        if (isfield(mesh.groups, elem_type))
          grp_data = getfield(groups, elem_type);
          mgrp_data = getfield(mesh.groups, elem_type);
          idx_group = numel(grp_data) + 1;

          for i=1:numel(mgrp_data)
            X = mesh.nodes(mgrp_data(i).nodes, 1:3).';

            if (all_nodes_in_group(X, group_defs(j), tol))
              grp_data(idx_group).id = group_defs(j).id;
              grp_data(idx_group).name = group_defs(j).name;
              grp_data(idx_group).nodes(end + (1:numel(mgrp_data(i).nodes))) = mgrp_data(i).nodes;
              grp_data(idx_group).elements(end + (1:numel(mgrp_data(i).elements))) = mgrp_data(i).elements;
              group_defs_num_nodes += numel(mgrp_data(i).nodes);
            endif
          endfor

          groups = setfield(groups, elem_type, grp_data);
        endif
      case "nodes"
        X = mesh.nodes(:, 1:3).';

        node_idx = nodes_in_group(X, group_defs(j), tol);

        if (numel(node_idx))
          groups.nodes(end + 1).id = group_defs(j).id;
          groups.nodes(end).name = group_defs(j).name;
          groups.nodes(end).nodes = node_idx;
          group_defs_num_nodes = numel(node_idx);
        endif
      otherwise
        error("unknown element type: \"%s\"", elem_type);
    endswitch

    if (group_defs_num_nodes == 0)
      error("no nodes found for group \"%s\"", group_defs(j).name);
    endif
  endfor

  for j=1:numel(elem_types)
    grp_data = getfield(groups, elem_types{j});

    for i=1:numel(grp_data)
      grp_data(i).nodes = unique(grp_data(i).nodes);
      grp_data(i).elements = unique(grp_data(i).elements);
    endfor

    groups = setfield(groups, elem_types{j}, grp_data);
  endfor
endfunction

function status = nodes_within_cylinder(dX, group_def, tol)
  r = norm(dX(1:2, :), "columns");
  status = r >= group_def.geometry.rmin - tol & ...
           r <= group_def.geometry.rmax + tol & ...
           dX(3, :) >= group_def.geometry.zmin - tol & ...
           dX(3, :) <= group_def.geometry.zmax + tol;
endfunction

function status = nodes_within_box(dX, group_def, tol)
  status = dX(1, :) >= group_def.geometry.xmin - tol & ...
           dX(1, :) <= group_def.geometry.xmax + tol & ...
           dX(2, :) >= group_def.geometry.ymin - tol & ...
           dX(2, :) <= group_def.geometry.ymax + tol & ...
           dX(3, :) >= group_def.geometry.zmin - tol & ...
           dX(3, :) <= group_def.geometry.zmax + tol;
endfunction

function status = all_nodes_in_group(X, group_def, tol)
  status = false;

  switch group_def.type
    case {"cylinder", "box"}
      dX = group_def.R.' * (X - group_def.X0);
  endswitch

  switch group_def.type
    case "cylinder"
      status =  all(nodes_within_cylinder(dX, group_def, tol));
    case "box"
      status = all(nodes_within_box(dX, group_def, tol));
    case "all"
      status = true;
    otherwise
      error("unknown geometry \"%s\"", group_def.geometry);
  endswitch
endfunction

function idx = nodes_in_group(X, group_def, tol)
  switch group_def.type
    case {"cylinder", "box"}
      dX = group_def.R.' * (X - group_def.X0);
  endswitch

  switch group_def.type
    case "cylinder"
      idx = find(nodes_within_cylinder(dX, group_def, tol));
    case "box"
      idx = find(nodes_within_box(dX, group_def, tol));
    case "all"
      idx = 1:columns(X);
    otherwise
      error("unknown geometry \"%s\"", group_def.geometry);
  endswitch
endfunction
