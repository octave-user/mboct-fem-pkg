## Copyright (C) 2011(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{load_case}] = fem_pre_mesh_struct_create(@var{geometry}, @var{loads}, @var{material}, @var{options})
## Create a structured hexahedral mesh of arbitrary shape, defined by callback functions.
##
## @var{geometry}.mesh_size.r @dots{} Natural nodal coordinates in r-direction
##
## @var{geometry}.mesh_size.s @dots{} Natural nodal coordinates in s-direction
##
## @var{geometry}.mesh_size.t @dots{} Natural nodal coordinates in t-direction
##
## @var{geometry}.spatial_coordinates @dots{} Callback function which returns the spatial coordinates (x, y, z) related to natural coordinates (r, s, t).
##
## @var{geometry}.material_selector @dots{} Callback function which returns the material number related to the element (r(1):r(2), s(1):s(2), t(1):t(2)).
## If the callback function returns zero, no element is created at this location.
##
## @var{geometry}.pressure_boundary_condition @dots{} Callback function which returns nodal pressure values at surface elements
##
## @var{geometry}.sewing.tolerance @dots{} If the distance between two nodes is smaller than this value, the same node number will be assigned to both nodes.
##
## @var{geometry}.boundary_conditions @dots{} Callback function which returns nodal displacement and force boundary conditions.
##
## @end deftypefn

function [mesh, load_case] = fem_pre_mesh_struct_create(geometry, loads, material, options)
  if (nargin < 3 || nargin > 4 || nargout > 2)
    print_usage();
  endif

  if (nargin < 4)
    options = struct();
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  if (~isfield(options, "identical_boundary_cond"))
    options.identical_boundary_cond = false;
  endif

  if (~isfield(options, "elem_type"))
    options.elem_type = "iso8";
  endif

  switch (options.elem_type)
    case {"iso8", "iso20"}
    otherwise
      error("unknown option elem_type=\"%s\"", options.elem_type);
  endswitch

  if (isfield(geometry.mesh_size, "r") && isfield(geometry.mesh_size, "s") && isfield(geometry.mesh_size, "t"))
    r = geometry.mesh_size.r;
    s = geometry.mesh_size.s;
    t = geometry.mesh_size.t;
  else
    error("geometry.mesh_size must contain eighter num_no_r, num_no_s and num_no_t or r, s and t");
  endif

  if (options.verbose)
    fprintf(stderr, "creating nodes ...\n");
    tic();
  endif

  switch (options.elem_type)
    case "iso8"
      x = y = z = zeros(numel(r), numel(s), numel(t));

      for i=1:numel(r)
        for j=1:numel(s)
          for k=1:numel(t)
            [x(i, j, k), y(i, j, k), z(i, j, k)] = feval(geometry.spatial_coordinates, r(i), s(j), t(k), geometry);
          endfor
        endfor
      endfor

      n = numel(r) * numel(s) * numel(t);
    case "iso20"
      x = y = z = nan(2 * numel(r) - 1, 2 * numel(s) - 1, 2 * numel(t) - 1);
      n = 0;

      idx = [zeros(1,3);
             eye(3)];

      for i=1:numel(r)
        for j=1:numel(s)
          for k=1:numel(t)
            for l=1:rows(idx)
              ii = 2 * i - 1 + idx(l, 1);
              jj = 2 * j - 1 + idx(l, 2);
              kk = 2 * k - 1 + idx(l, 3);

              if (ii > 2 * numel(r) - 1 || jj > 2 * numel(s) - 1 || kk > 2 * numel(t) - 1)
                continue;
              endif

              if (idx(l, 1))
                dr = 0.5 * (r(i + 1) - r(i));
              else
                dr = 0;
              endif

              if (idx(l, 2))
                ds = 0.5 * (s(j + 1) - s(j));
              else
                ds = 0;
              endif

              if (idx(l, 3))
                dt = 0.5 * (t(k + 1) - t(k));
              else
                dt = 0;
              endif

              [x(ii, jj, kk), ...
               y(ii, jj, kk), ...
               z(ii, jj, kk)] = feval(geometry.spatial_coordinates, ...
                                      r(i) + dr, ...
                                      s(j) + ds, ...
                                      t(k) + dt, ...
                                      geometry);
              ++n;
            endfor
          endfor
        endfor
      endfor
  endswitch

  if (options.verbose)
    toc();
    fprintf(stderr, "merging nodes ...\n");
    tic();
  endif

  mesh.nodes = zeros(n, 6);
  inode_idx = zeros(size(x), "int32");

  inode = int32(0);

  switch (options.elem_type)
    case "iso8"
      for i=1:numel(r)
        for j=1:numel(s)
          for k=1:numel(t)
            sewing_cond = sqrt((x(i, j, k) - x).^2 + (y(i, j, k) - y).^2 + (z(i, j, k) - z).^2) < geometry.sewing.tolerance;
            sewing_cond(i, j, k) = false;
            inode2 = find_sewing_node(sewing_cond, inode_idx);
            if (inode2 == 0)
              inode_idx(i, j, k) = ++inode;

              mesh.nodes(inode, 1:3) = [x(i, j, k), ...
                                        y(i, j, k), ...
                                        z(i, j, k)];
            else
              inode_idx(i, j, k) = inode2;
            endif
          endfor
        endfor
      endfor
    case "iso20"
      for i=1:2:2 * numel(r) - 1
        for j=1:2:2 * numel(s) - 1
          for k=1:2:2 * numel(t) - 1
            for l=1:rows(idx)
              ii = i + idx(l, 1);
              jj = j + idx(l, 2);
              kk = k + idx(l, 3);

              if (ii > 2 * numel(r) - 1 || jj > 2 * numel(s) - 1 || kk > 2 * numel(t) - 1)
                continue;
              endif

              sewing_cond = sqrt((x(ii, jj, kk) - x).^2 + ...
                                 (y(ii, jj, kk) - y).^2 + ...
                                 (z(ii, jj, kk) - z).^2) < geometry.sewing.tolerance;
              sewing_cond(ii, jj, kk) = false;
              inode2 = find_sewing_node(sewing_cond, inode_idx);

              if (inode2 == 0)
                inode_idx(ii, jj, kk) = ++inode;

                mesh.nodes(inode, 1:3) = [x(ii, jj, kk), ...
                                          y(ii, jj, kk), ...
                                          z(ii, jj, kk)];
              else
                inode_idx(ii, jj, kk) = inode2;
              endif
            endfor
          endfor
        endfor
      endfor
  endswitch

  mesh.nodes = mesh.nodes(1:inode, :);

  if (options.verbose)
    toc();
    fprintf(stderr, "creating elements ...\n");
    tic();
  endif

  switch (options.elem_type)
    case "iso8"
      mesh.elements.iso8 = zeros((numel(r) - 1) * (numel(s) - 1) * (numel(t) - 1), 8, "int32");
    case "iso20"
      mesh.elements.iso20 = zeros((numel(r) - 1) * (numel(s) - 1) * (numel(t) - 1), 20, "int32");
  endswitch

  mesh.elements.rbe3 = struct("nodes", [], "weight", [])([]);
  mesh.material_data = material;

  switch (options.elem_type)
    case "iso8"
      mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
    case "iso20"
      mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
  endswitch

  ielement = int32(0);

  for i=2:numel(r)
    for j=2:numel(s)
      for k=2:numel(t)
        imat_id = feval(geometry.material_selector, r(i - 1:i), s(j - 1:j), t(k - 1:k), geometry);

        if (imat_id ~= 0)
          if (imat_id < 1 || imat_id > numel(material))
            error("invalid material number %d returned for %d, %d, %d", imat_id, i, j, k);
          endif

          switch (options.elem_type)
            case "iso8"
              mesh.materials.iso8(++ielement) = imat_id;
              mesh.elements.iso8(ielement, :) = [inode_idx(i,     j,     k), ...
                                                 inode_idx(i - 1, j,     k), ...
                                                 inode_idx(i - 1, j - 1, k), ...
                                                 inode_idx(i,     j - 1, k), ...
                                                 inode_idx(i,     j,     k - 1), ...
                                                 inode_idx(i - 1, j,     k - 1), ...
                                                 inode_idx(i - 1, j - 1, k - 1), ...
                                                 inode_idx(i,     j - 1, k - 1)];
            case "iso20"
              mesh.materials.iso20(++ielement) = imat_id;
              mesh.elements.iso20(ielement, :) = [inode_idx(2 * i - 1, 2 * j - 1, 2 * k - 1), ...
                                                  inode_idx(2 * i - 3, 2 * j - 1, 2 * k - 1), ...
                                                  inode_idx(2 * i - 3, 2 * j - 3, 2 * k - 1), ...
                                                  inode_idx(2 * i - 1, 2 * j - 3, 2 * k - 1), ...
                                                  inode_idx(2 * i - 1, 2 * j - 1, 2 * k - 3), ...
                                                  inode_idx(2 * i - 3, 2 * j - 1, 2 * k - 3), ...
                                                  inode_idx(2 * i - 3, 2 * j - 3, 2 * k - 3), ...
                                                  inode_idx(2 * i - 1, 2 * j - 3, 2 * k - 3), ...
                                                  inode_idx(2 * i - 2, 2 * j - 1, 2 * k - 1), ...
                                                  inode_idx(2 * i - 3, 2 * j - 2, 2 * k - 1), ...
                                                  inode_idx(2 * i - 2, 2 * j - 3, 2 * k - 1), ...
                                                  inode_idx(2 * i - 1, 2 * j - 2, 2 * k - 1), ...
                                                  inode_idx(2 * i - 2, 2 * j - 1, 2 * k - 3), ...
                                                  inode_idx(2 * i - 3, 2 * j - 2, 2 * k - 3), ...
                                                  inode_idx(2 * i - 2, 2 * j - 3, 2 * k - 3), ...
                                                  inode_idx(2 * i - 1, 2 * j - 2, 2 * k - 3), ...
                                                  inode_idx(2 * i - 1, 2 * j - 1, 2 * k - 2), ...
                                                  inode_idx(2 * i - 3, 2 * j - 1, 2 * k - 2), ...
                                                  inode_idx(2 * i - 3, 2 * j - 3, 2 * k - 2), ...
                                                  inode_idx(2 * i - 1, 2 * j - 3, 2 * k - 2)];
          endswitch
        endif
      endfor
    endfor
  endfor

  switch (options.elem_type)
    case "iso8"
      mesh.elements.iso8 = mesh.elements.iso8(1:ielement, :);
      mesh.materials.iso8 = mesh.materials.iso8(1:ielement);
    case "iso20"
      mesh.elements.iso20 = mesh.elements.iso20(1:ielement, :);
      mesh.materials.iso20 = mesh.materials.iso20(1:ielement);
  endswitch

  nodes_in_use = false(rows(mesh.nodes), 1);

  for i=1:rows(mesh.nodes)
    if (numel(find(getfield(mesh.elements, options.elem_type) == i)) > 0)
      nodes_in_use(i) = true;
    endif
  endfor

  if (options.verbose)
    toc();
    fprintf(stderr, "creating load cases ...\n");
    tic();
  endif

  if (options.identical_boundary_cond)
    load_case = create_load_case(mesh, geometry, inode_idx, r, s, t, loads(1), nodes_in_use, options);
    load_case = fem_pre_load_case_merge(load_case, repmat(setfield(load_case, "locked_dof", []), 1, numel(loads) - 1));
  else
    empty_cell = cell(1, numel(loads));
    load_case = struct("locked_dof", empty_cell, "loaded_nodes", empty_cell, "loads", empty_cell);
    for l=1:numel(loads)
      load_case(l) = create_load_case(mesh, geometry, inode_idx, r, s, t, loads(l), nodes_in_use, options);
    endfor
  endif

  if (options.verbose)
    toc();
    fprintf(stderr, "creating pressure loads ...\n");
    tic();
  endif

  switch (options.elem_type)
    case "iso8"
      pressure_perm_idx = int32([3, 4;
                                 2, 1]);
    case "iso20"
      pressure_perm_idx = int32([3, 7, 4;
                                 6, 9, 8;
                                 2, 5, 1]);
  endswitch

  for l=1:numel(loads)
    if (isfield(geometry, "pressure_boundary_condition"))
      ielement_press = int32(0);

      switch (options.elem_type)
        case "iso8"
          load_case(l).pressure.iso4.elements = zeros(numel(r) * (numel(s) - 1) * (numel(t) - 1), 4, "int32");
          load_case(l).pressure.iso4.p = zeros(rows(load_case(l).pressure.iso4.elements), 4);
        case "iso20"
          load_case(l).pressure.quad8.elements = zeros(numel(r) * (numel(s) - 1) * (numel(t) - 1), 8, "int32");
          load_case(l).pressure.quad8.p = zeros(rows(load_case(l).pressure.quad8.elements), 8);
      endswitch

      if (isfield(loads(l), "idx_r") && numel(loads(l).idx_r) > 0)
        idx_r = loads(l).idx_r;
      else
        idx_r = 1:numel(r);
      endif

      if (isfield(loads(l), "idx_s") && numel(loads(l).idx_s) > 0)
        idx_s = loads(l).idx_s;
      else
        idx_s = 1:numel(s);
      endif

      if (isfield(loads(l), "idx_t") && numel(loads(l).idx_t) > 0)
        idx_t = loads(l).idx_t;
      else
        idx_t = 1:numel(t);
      endif

      switch (options.elem_type)
        case "iso8"
          for i=idx_r
            for j=idx_s
              for k=idx_t
                if (i < 1 || j < 2 || k < 2 || i > numel(r) || j > numel(s) || k > numel(t))
                  continue;
                endif

                p_elem = feval(geometry.pressure_boundary_condition, ...
                               r(i), ...
                               s(j - 1:j), ...
                               t(k - 1:k), ...
                               geometry, loads(l), ...
                               pressure_perm_idx);

                if (numel(p_elem) ~= 0)
                  ++ielement_press;
                  load_case(l).pressure.iso4.p(ielement_press, :) = p_elem;
                  load_case(l).pressure.iso4.elements(ielement_press, :) = [inode_idx(i,     j,     k), ...
                                                                            inode_idx(i,     j, k - 1), ...
                                                                            inode_idx(i, j - 1, k - 1), ...
                                                                            inode_idx(i, j - 1,     k)];
                endif
              endfor
            endfor
          endfor
          load_case(l).pressure.iso4.elements = load_case(l).pressure.iso4.elements(1:ielement_press, :);
          load_case(l).pressure.iso4.p = load_case(l).pressure.iso4.p(1:ielement_press, :);
        case "iso20"
          for i=1:numel(idx_r)
            for j=1:numel(idx_s)
              for k=1:numel(idx_t)
                ii = idx_r(i);
                jj = idx_s(j);
                kk = idx_t(k);

                if (ii < 1 || jj < 2 || kk < 2 || ii > numel(r) || jj > numel(s) || kk > numel(t))
                  continue;
                endif

                p_elem = feval(geometry.pressure_boundary_condition, ...
                               r(ii), ...
                               linspace(s(jj - 1), s(jj), 3), ...
                               linspace(t(kk - 1), t(kk), 3), ...
                               geometry, loads(l), ...
                               pressure_perm_idx);

                if (numel(p_elem))
                  ++ielement_press;
                  load_case(l).pressure.quad8.p(ielement_press, :) = p_elem(1:8);
                  load_case(l).pressure.quad8.elements(ielement_press, :) = [inode_idx(2 * ii - 1, 2 * jj - 1, 2 * kk - 1), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 1, 2 * kk - 3), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 3, 2 * kk - 3), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 3, 2 * kk - 1), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 1, 2 * kk - 2), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 2, 2 * kk - 3), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 3, 2 * kk - 2), ...
                                                                             inode_idx(2 * ii - 1, 2 * jj - 2, 2 * kk - 1)];
                endif
              endfor
            endfor
          endfor
          load_case(l).pressure.quad8.elements = load_case(l).pressure.quad8.elements(1:ielement_press, :);
          load_case(l).pressure.quad8.p = load_case(l).pressure.quad8.p(1:ielement_press, :);
      endswitch
    endif
  endfor

  mesh.structured.r = r;
  mesh.structured.s = s;
  mesh.structured.t = t;
  mesh.structured.inode_idx = inode_idx;

  switch (options.elem_type)
    case "iso20"
      mesh = replace_degenerated(mesh);
  endswitch
endfunction

function inode = find_sewing_node(sewing_cond, inode_idx)
  inode = int32(0);

  for k2=1:size(sewing_cond, 3)
    [i2, j2] = find(sewing_cond(:, :, k2));
    if (numel(i2) > 0)
      inode = inode_idx(i2(1), j2(1), k2);
      return;
    endif
  endfor
endfunction

function load_case = create_load_case(mesh, geometry, inode_idx, r, s, t, load, nodes_in_use, options)
  load_case.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
  load_case.loads = zeros(rows(mesh.nodes), columns(mesh.nodes));
  load_case.loaded_nodes = zeros(rows(mesh.nodes), 1, "int32");

  iload = int32(0);

  if (isfield(load, "idx_r") && numel(load.idx_r) > 0)
    idx_r = load.idx_r;
  else
    idx_r = 1:numel(r);
  endif

  if (isfield(load, "idx_s") && numel(load.idx_s) > 0)
    idx_s = load.idx_s;
  else
    idx_s = 1:numel(s);
  endif

  if (isfield(load, "idx_t") && numel(load.idx_t) > 0)
    idx_t = load.idx_t;
  else
    idx_t = 1:numel(t);
  endif

  switch (options.elem_type)
    case "iso8"
      for i=idx_r
        for j=idx_s
          for k=idx_t
            [Fijk, locked_dof_ijk] = feval(geometry.boundary_condition, r(i), s(j), t(k), geometry, load);

            if (numel(locked_dof_ijk) ~= 0)
              load_case.locked_dof(inode_idx(i, j, k), 1:numel(locked_dof_ijk)) = locked_dof_ijk;
            endif

            if (~nodes_in_use(inode_idx(i, j, k)))
              load_case.locked_dof(inode_idx(i, j, k), :) = true; % eliminate nodes not connected to elements
            endif

            if (numel(Fijk) ~= 0)
              load_case.loads(++iload, 1:numel(Fijk)) = Fijk;
              load_case.loaded_nodes(iload) = inode_idx(i, j, k);
            endif
          endfor
        endfor
      endfor
    case "iso20"
      idx = [zeros(1, 3);
             eye(3)];

      for i=1:numel(idx_r)
        for j=1:numel(idx_s)
          for k=1:numel(idx_t)
            for l=1:rows(idx)
              ii = 2 * idx_r(i) - 1 + idx(l, 1);
              jj = 2 * idx_s(j) - 1 + idx(l, 2);
              kk = 2 * idx_t(k) - 1 + idx(l, 3);

              if (ii > 2 * numel(r) - 1 || jj > 2 * numel(s) - 1 || kk > 2 * numel(t) - 1)
                continue;
              endif

              if (idx(l, 1))
                dr = 0.5 * (r(idx_r(i) + 1) - r(idx_r(i)));
              else
                dr = 0;
              endif

              if (idx(l, 2))
                ds = 0.5 * (s(idx_s(j) + 1) - s(idx_s(j)));
              else
                ds = 0;
              endif

              if (idx(l, 3))
                dt = 0.5 * (t(idx_t(k) + 1) - t(idx_t(k)));
              else
                dt = 0;
              endif

              [Fijk, locked_dof_ijk] = feval(geometry.boundary_condition, ...
                                             r(idx_r(i)) + dr, ...
                                             s(idx_s(j)) + ds, ...
                                             t(idx_t(k)) + dt, ...
                                             geometry, ...
                                             load);

              if (numel(locked_dof_ijk) ~= 0)
                load_case.locked_dof(inode_idx(ii, jj, kk), 1:numel(locked_dof_ijk)) = locked_dof_ijk;
              endif

              if (~nodes_in_use(inode_idx(ii, jj, kk)))
                load_case.locked_dof(inode_idx(ii, jj, kk), :) = true; % eliminate nodes not connected to elements
              endif

              if (numel(Fijk) ~= 0)
                load_case.loads(++iload, 1:numel(Fijk)) = Fijk;
                load_case.loaded_nodes(iload) = inode_idx(ii, jj, kk);
              endif
            endfor
          endfor
        endfor
      endfor
  endswitch

  load_case.loads = load_case.loads(1:iload, :);
  load_case.loaded_nodes = load_case.loaded_nodes(1:iload);
endfunction

function mesh_rep = replace_degenerated(mesh)
  mesh_rep = mesh;
  mesh_rep.elements.iso20 = zeros(size(mesh.elements.iso20), "int32");
  mesh_rep.elements.penta15 = zeros(rows(mesh.elements.iso20), 15, "int32");
  mesh_rep.materials.iso20 = zeros(size(mesh.materials.iso20), "int32");
  mesh_rep.materials.penta15 = zeros(size(mesh.materials.iso20), "int32");

  inum_iso20 = int32(0);
  inum_penta15 = int32(0);
  idxn = 1:columns(mesh.elements.iso20);

  for i=1:rows(mesh.elements.iso20)
    elnodes = unique(mesh.elements.iso20(i, :));
    idxn = [];

    if (numel(elnodes) < columns(mesh.elements.iso20))
      idx = zeros(1, columns(mesh.elements.iso20));
      for j=1:columns(mesh.elements.iso20)
        neq = (mesh.elements.iso20(i, :) == mesh.elements.iso20(i, j));
        neq(j) = false;
        idxtmp = find(neq);
        if (numel(idxtmp))
          idx(j) = idxtmp(1);
        endif
      endfor

      idxs = sort(find(idx));

      switch (idxs)
        case [3, 4, 7, 8, 11, 15, 19, 20]
          idxn = [2,3,6,1,4,5,10,14,18,12,16,17,9,11,13];
        case [1, 2, 5, 6, 9, 13, 17, 18]
          idxn = [2,3,7,1,4,8,10,19,14,12,20,16,9,11,15];
        case [5, 6, 7, 8, 13, 14, 15, 16]
          idxn = [2, 3, 6, 1, 4, 5, 10, 19, 18, 12, 20, 17, 9, 11, 13];
        case [1, 2, 3, 4, 9, 10, 11, 12]
          idxn = [7, 6, 2, 8, 5, 1, 14, 18, 19, 16, 17, 20, 15, 13, 9];
      endswitch
    endif

    if (numel(idxn))
      mesh_rep.elements.penta15(++inum_penta15, :) = mesh.elements.iso20(i, idxn);
      mesh_rep.materials.penta15(inum_penta15) = mesh.materials.iso20(i);
    else
      mesh_rep.elements.iso20(++inum_iso20, :) = mesh.elements.iso20(i, :);
      mesh_rep.materials.iso20(inum_iso20) = mesh.materials.iso20(i);
    endif
  endfor

  mesh_rep.elements.iso20 = mesh_rep.elements.iso20(1:inum_iso20, :);
  mesh_rep.materials.iso20 = mesh_rep.materials.iso20(1:inum_iso20);
  mesh_rep.elements.penta15 = mesh_rep.elements.penta15(1:inum_penta15, :);
  mesh_rep.materials.penta15 = mesh_rep.materials.penta15(1:inum_penta15);
endfunction

%!test
%! ## TEST1
%! close all;
%! scale_stat = 2e-3;
%! scale_eig = 2e-3;
%! tol = eps;
%! number_of_modes = 3;
%! f_run_post_proc = false;
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t, varagin)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! function [F, locked] = boundary_cond(r, s, t, geo, load, varargin)
%!   locked = [];
%!   F = [];
%!   [x, y, z, R, Phi] = cylinder_geo(geo.user_data.cylinder, r, s, t);
%!   if (~load.flags.use_pressure_boundary_cond && r == 0)
%!     A = (geo.user_data.cylinder.z2 - geo.user_data.cylinder.z1) * 2 * geo.user_data.cylinder.r1 * pi ...
%!         / ((geo.mesh_size.num_no_t - 1) * (geo.mesh_size.num_no_s - 1));
%!     if (t == 0 || t == 1)
%!       alpha = 0.5;
%!     else
%!       alpha = 1;
%!     endif
%!     if (s == 1)
%!       alpha = 0;
%!     endif
%!     F = alpha * load.p * A * [cos(Phi); sin(Phi); 0];
%!   endif
%!   if ((t == 1 || t == 0) && r == 1)
%!     locked = true(3, 1);
%!   endif
%! endfunction
%!
%! function p = pressure_boundary_cond(r, s, t, geometry, load, perm_idx, varargin)
%!   if (load.flags.use_pressure_boundary_cond && r == 0)
%!     p(perm_idx) = load.p;
%!   else
%!     p = [];
%!   endif
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, 5);
%! geometry.mesh_size.s = linspace(0, 1, 16);
%! geometry.mesh_size.t = linspace(0, 1, 10);
%! geometry.user_data.cylinder.r1 = 2e-3;
%! geometry.user_data.cylinder.r2 = 7e-3;
%! geometry.user_data.cylinder.z1 = -8e-3;
%! geometry.user_data.cylinder.z2 = 8e-3;
%! geometry.user_data.cylinder.Phi1 = 0;
%! geometry.user_data.cylinder.Phi2 = 2*pi;
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.cylinder.r2;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition = @(r, s, t, geo, load, varargin) feval("boundary_cond", r, s, t, geo, load, varargin);
%! geometry.pressure_boundary_condition = @(r, s, t, geo, load, perm_idx, varargin) feval("pressure_boundary_cond", r, s, t, geo, load, perm_idx, varargin);
%! material.E = 210000e6;
%! material.nu = 0.26;
%! material.rho = 7850;
%! load.p = 1e6;
%! load.flags.use_pressure_boundary_cond = true;
%! options.elem_type = "iso20";
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load, material, options);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                               load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%! opts.scale_def = 0.5 * geometry.user_data.cylinder.r2 / max(max(max(abs(sol_eig.def))));
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   opts.rotation_angle = [30, 30, 0] * pi / 180;
%!   fem_post_sol_external(mesh, sol_eig, opts);
%!   fn = dir([opts.print_to_file, "_*.jpg"]);
%!   for i=1:numel(fn)
%!     [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("mode %d %.0fHz", i, sol_eig.f(i)));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     fn = dir([opts.print_to_file, "_*.jpg"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! figure_list();
%! endif

%!test
%! ## TEST2
%! close all;
%! scale_stat = 2e-3;
%! scale_eig = 2e-3;
%! tol = eps;
%! number_of_modes = 3;
%! f_run_post_proc = false;
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t, varargin)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! function [F, locked] = boundary_cond(r, s, t, geo, load, varargin)
%!   locked = [];
%!   F = [];
%!   [x, y, z, R, Phi] = cylinder_geo(geo.user_data.cylinder, r, s, t, varargin);
%!   if (~load.flags.use_pressure_boundary_cond && r == 0)
%!     A = (geo.user_data.cylinder.z2 - geo.user_data.cylinder.z1) * 2 * geo.user_data.cylinder.r1 * pi ...
%!         / ((geo.mesh_size.num_no_t - 1) * (geo.mesh_size.num_no_s - 1));
%!     if (t == 0 || t == 1)
%!       alpha = 0.5;
%!     else
%!       alpha = 1;
%!     endif
%!     if (s == 1)
%!       alpha = 0;
%!     endif
%!     F = alpha * load.p * A * [cos(Phi); sin(Phi); 0];
%!   endif
%!   if ((t == 1 || t == 0) && r == 1)
%!     locked = true(3, 1);
%!   endif
%! endfunction
%!
%! function p = pressure_boundary_cond(r, s, t, geometry, load, perm_idx, varargin)
%!   if (load.flags.use_pressure_boundary_cond && r == 0)
%!     p(perm_idx) = load.p;
%!   else
%!     p = [];
%!   endif
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, 5);
%! geometry.mesh_size.s = linspace(0, 1, 16);
%! geometry.mesh_size.t = linspace(0, 1, 10);
%! geometry.user_data.cylinder.r1 = 2e-3;
%! geometry.user_data.cylinder.r2 = 7e-3;
%! geometry.user_data.cylinder.z1 = -8e-3;
%! geometry.user_data.cylinder.z2 = 8e-3;
%! geometry.user_data.cylinder.Phi1 = 0;
%! geometry.user_data.cylinder.Phi2 = 2*pi;
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.cylinder.r2;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition = @(r, s, t, geo, load, varargin) feval("boundary_cond", r, s, t, geo, load, varargin);
%! geometry.pressure_boundary_condition = @(r, s, t, geo, load, perm_idx, varargin) feval("pressure_boundary_cond", r, s, t, geo, load, perm_idx, varargin);
%! material.E = 210000e6;
%! material.nu = 0.26;
%! material.rho = 7850;
%! load.p = 1e6;
%! load.flags.use_pressure_boundary_cond = true;
%! options.elem_type = "iso8";
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load, material, options);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                               load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%! opts.scale_def = 0.5 * geometry.user_data.cylinder.r2 / max(max(max(abs(sol_eig.def))));
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   opts.rotation_angle = [30, 30, 0] * pi / 180;
%!   fem_post_sol_external(mesh, sol_eig, opts);
%!   fn = dir([opts.print_to_file, "_*.jpg"]);
%!   for i=1:numel(fn)
%!     [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("mode %d %.0fHz", i, sol_eig.f(i)));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     fn = dir([opts.print_to_file, "_*.jpg"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! figure_list();
%! endif

%!test
%! ## TEST3
%! close all;
%! scale_eig = 0.15;
%! number_of_modes = 14;
%! f_run_post_proc = false;
%! ## Code_Aster SDLS109 V2.3.109
%! Rm = 0.369;
%! t = 0.048;
%! L = 0.05;
%! E = 185000e6;
%! nu = 0.3;
%! rho = 7800;
%! h = t / 8;
%! fref = [zeros(1, 6), 210.55, 210.55, 587.92, 587.92, 205.89, 205.89, 588.88, 588.88];
%!
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t, varargin)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, max([2, ceil(t/h)]));
%! geometry.mesh_size.s = linspace(0, 1, max([2, ceil(2 * pi * Rm / h)]));
%! geometry.mesh_size.t = linspace(0, 1, max([2, ceil(L / h)]));
%! geometry.user_data.cylinder.r1 = Rm - 0.5 * t;
%! geometry.user_data.cylinder.r2 = Rm + 0.5 * t;
%! geometry.user_data.cylinder.z1 = -0.5 * L;
%! geometry.user_data.cylinder.z2 = 0.5 * L;
%! geometry.user_data.cylinder.Phi1 = 0;
%! geometry.user_data.cylinder.Phi2 = 2 * pi;
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.cylinder.r2;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition =  @(r, s, t, varargin) {[],[]}{:};
%! material.E = E;
%! material.nu = nu;
%! material.rho = rho;
%! load_data = struct();
%! options.elem_type = "iso8";
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_info, ...
%!  mesh_info] = fem_ass_matrix(mesh, ...
%!                             dof_map, ...
%!                             [FEM_MAT_MASS, ...
%!                              FEM_MAT_STIFFNESS], ...
%!                             load_case);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes, 100);
%! opts.scale_def = 0.5 * Rm / max(max(max(abs(sol_eig.def))));
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   opts.rotation_angle = [0, 0, 0];
%!   opts.output_step_idx = [7, 9, 11];
%!   fem_post_sol_external(mesh, sol_eig, opts);
%!   fn = dir([opts.print_to_file, "_*.jpg"]);
%!   for i=1:numel(fn)
%!     [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("mode %d %.0fHz", i, sol_eig.f(opts.output_step_idx(i))));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     fn = dir([opts.print_to_file, "_*.jpg"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! figure_list();
%! endif
%! tol = 0.5e-2;
%! assert_simple(sol_eig.f(1:numel(fref)), sort(fref), tol * max(fref));

%!test
%! ## TEST 4
%! close all;
%! scale_eig = 0.15;
%! number_of_modes = 14;
%! f_run_post_proc = false;
%! ## Code_Aster SDLS109 V2.3.109
%! Rm = 0.369;
%! t = 0.048;
%! L = 0.05;
%! E = 185000e6;
%! nu = 0.3;
%! rho = 7800;
%! h = t / 6;
%! fref = [zeros(1, 6), 210.55, 210.55, 587.92, 587.92, 205.89, 205.89, 588.88, 588.88];
%!
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t, varargin)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, max([2, ceil(t/h)]));
%! geometry.mesh_size.s = linspace(0, 1, max([2, ceil(2 * pi * Rm / h)]));
%! geometry.mesh_size.t = linspace(0, 1, max([2, ceil(L / h)]));
%! geometry.user_data.cylinder.r1 = Rm - 0.5 * t;
%! geometry.user_data.cylinder.r2 = Rm + 0.5 * t;
%! geometry.user_data.cylinder.z1 = -0.5 * L;
%! geometry.user_data.cylinder.z2 = 0.5 * L;
%! geometry.user_data.cylinder.Phi1 = 0;
%! geometry.user_data.cylinder.Phi2 = 2 * pi;
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.cylinder.r2;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition =  @(r, s, t, varargin) {[],[]}{:};
%! material.E = E;
%! material.nu = nu;
%! material.rho = rho;
%! load_data = struct();
%! options.elem_type = "iso20";
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_info, ...
%!  mesh_info] = fem_ass_matrix(mesh, ...
%!                             dof_map, ...
%!                             [FEM_MAT_MASS, ...
%!                              FEM_MAT_STIFFNESS], ...
%!                             load_case);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes, 100);
%! opts.scale_def = 0.5 * Rm / max(max(max(abs(sol_eig.def))));
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   opts.rotation_angle = [0, 0, 0];
%!   opts.output_step_idx = [7, 9, 11];
%!   fem_post_sol_external(mesh, sol_eig, opts);
%!   fn = dir([opts.print_to_file, "_*.jpg"]);
%!   for i=1:numel(fn)
%!     [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("mode %d %.0fHz", i, sol_eig.f(opts.output_step_idx(i))));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     fn = dir([opts.print_to_file, "_*.jpg"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! figure_list();
%! endif
%! tol = 0.5e-2;
%! assert_simple(sol_eig.f(1:numel(fref)), sort(fref), tol * max(fref));

%!test
%! ## TEST 5
%! ## Code_Aster SHLV100 V2.07.100
%! close all;
%! a = 0.1;
%! b = 0.2;
%! c = 0.01;
%! t = b - a;
%! h = t / 30;
%! E = 26;
%! nu = 0.3;
%! rho = 35;
%! p = 1;
%! omega = 0.2;
%! ur_a = 7.3398e-3;
%! ur_b = 4.681610e-3;
%! taur_a = -1;
%! taur_b = 0;
%! taut_a = 1.6685;
%! taut_b = 0.66738;
%! tauz_a = 0.20055;
%! tauz_b = 0.20031;
%! f_run_post_proc = false;
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t, varargin)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! function p = pressure_callback(r, s, t, geometry, load_data, perm_idx, varargin)
%!   p = [];
%!   if (r == 0)
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%! function [F, locked] = boundary_cond_callback(r, s, t, geometry, load_data, varargin)
%!   F = [];
%!   locked = false(1, 3);
%!   if (s == 0)
%!     locked(2) = true;
%!   endif
%!   if (s == 1)
%!     locked(1) = true;
%!   endif
%!   if (t == 0 || t == 1)
%!     locked(3) = true;
%!   endif
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, max([2, ceil(t / h)]));
%! geometry.mesh_size.s = linspace(0, 1, max([2, ceil(pi / 2 * mean([a, b]) / h)]));
%! geometry.mesh_size.t = linspace(0, 1, max([2, ceil(c / h)]));
%! geometry.user_data.cylinder.r1 = a;
%! geometry.user_data.cylinder.r2 = b;
%! geometry.user_data.cylinder.z1 = 0;
%! geometry.user_data.cylinder.z2 = c;
%! geometry.user_data.cylinder.Phi1 = 0;
%! geometry.user_data.cylinder.Phi2 = pi / 2;
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.cylinder.r2;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("boundary_cond_callback", r, s, t, geometry, load_data, varargin);
%! geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("pressure_callback", r, s, t, geometry, load_data, perm_idx, varargin);
%! material.E = E;
%! material.nu = nu;
%! material.rho = rho;
%! load_data.pressure = p;
%! options.elem_type = "iso8";
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%!
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_info, ...
%!  mesh_info] = fem_ass_matrix(mesh, ...
%!                             dof_map, ...
%!                             [FEM_MAT_MASS, ...
%!                              FEM_MAT_STIFFNESS, ...
%!                              FEM_VEC_LOAD_CONSISTENT], ...
%!                             load_case);
%!
%! ## Effectively solve (-omega^2 * M + K) * U = R
%! mat_ass.K += -omega^2 * mat_ass.M;
%!
%! sol_harm = fem_sol_static(mesh, dof_map, mat_ass);
%!
%! sol_harm.stress = fem_ass_matrix(mesh, dof_map, [FEM_VEC_STRESS_CAUCH], load_case, sol_harm);
%!
%! tol_u = 0.003;
%! tol_tau = 0.08;
%! taun = nan(rows(mesh.nodes), 6);
%! idxtens = int32([1, 4, 6; 4, 2, 5; 6, 5, 3]);
%!
%! for i=1:rows(sol_harm.stress.taum.iso8)
%!   for j=1:8
%!     taun(mesh.elements.iso8(i, j), :) = sol_harm.stress.taum.iso8(i, j, :);
%!   endfor
%! endfor
%! opts.scale_def = 0.5 * a / max(max(abs(sol_harm.def)));
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   opts.rotation_angle = [0, 0, 0];
%!   fem_post_sol_external(mesh, sol_harm, opts);
%!   [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!   figure("visible", "off");
%!   imshow(img, map);
%!   title("Gmsh - deformed mesh / continuous stress tensor");
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     unlink([opts.print_to_file, "_001.jpg"]);
%!   endif
%! end_unwind_protect
%! figure_list();
%! endif
%! for i=[1, size(mesh.structured.inode_idx, 1)]
%!   switch (i)
%!   case 1
%!     ur = ur_a;
%!     taur = taur_a;
%!     taut = taut_a;
%!     tauz = tauz_a;
%!   case size(mesh.structured.inode_idx, 1)
%!     ur = ur_b;
%!     taur = taur_b;
%!     taut = taut_b;
%!     tauz = tauz_b;
%!   endswitch
%!
%!   for j=1:size(mesh.structured.inode_idx, 2)
%!     for k=1:size(mesh.structured.inode_idx, 3)
%!       inode = mesh.structured.inode_idx(i, j, k);
%!       e1 = [mesh.nodes(inode, 1:2).'; 0];
%!       e3 = [0; 0; 1];
%!       e2 = cross(e3, e1);
%!       R = [e1, e2, e3];
%!       R *= diag(1 ./ norm(R, "cols"));
%!       U = R.' * sol_harm.def(inode, 1:3).';
%!       tau = R.' * taun(inode, :)(idxtens) * R;
%!       assert_simple(U, [ur; 0; 0], tol_u * abs(ur));
%!       assert_simple(tau, diag([taur, taut, tauz]), tol_tau * abs(p));
%!     endfor
%!   endfor
%! endfor

%!test
%! ## TEST6
%! ## Code_Aster SHLV100 V2.07.100
%! close all;
%! a = 0.1;
%! b = 0.2;
%! c = 0.01;
%! t = b - a;
%! h = t / 30;
%! E = 26;
%! nu = 0.3;
%! rho = 35;
%! p = 1;
%! omega = 0.2;
%! ur_a = 7.3398e-3;
%! ur_b = 4.681610e-3;
%! taur_a = -1;
%! taur_b = 0;
%! taut_a = 1.6685;
%! taut_b = 0.66738;
%! tauz_a = 0.20055;
%! tauz_b = 0.20031;
%! f_run_post_proc = false;
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t, varargin)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! function p = pressure_callback(r, s, t, geometry, load_data, perm_idx, varargin)
%!   p = [];
%!   if (r == 0)
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%! function [F, locked] = boundary_cond_callback(r, s, t, geometry, load_data, varargin)
%!   F = [];
%!   locked = false(1, 3);
%!   if (s == 0)
%!     locked(2) = true;
%!   endif
%!   if (s == 1)
%!     locked(1) = true;
%!   endif
%!   if (t == 0 || t == 1)
%!     locked(3) = true;
%!   endif
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, max([2, ceil(t / h)]));
%! geometry.mesh_size.s = linspace(0, 1, max([2, ceil(pi / 2 * mean([a, b]) / h)]));
%! geometry.mesh_size.t = linspace(0, 1, max([2, ceil(c / h)]));
%! geometry.user_data.cylinder.r1 = a;
%! geometry.user_data.cylinder.r2 = b;
%! geometry.user_data.cylinder.z1 = 0;
%! geometry.user_data.cylinder.z2 = c;
%! geometry.user_data.cylinder.Phi1 = 0;
%! geometry.user_data.cylinder.Phi2 = pi / 2;
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.cylinder.r2;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("boundary_cond_callback", r, s, t, geometry, load_data, varargin);
%! geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("pressure_callback", r, s, t, geometry, load_data, perm_idx, varargin);
%! material.E = E;
%! material.nu = nu;
%! material.rho = rho;
%! load_data.pressure = p;
%! options.elem_type = "iso20";
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%!
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_info, ...
%!  mesh_info] = fem_ass_matrix(mesh, ...
%!                             dof_map, ...
%!                             [FEM_MAT_MASS, ...
%!                              FEM_MAT_STIFFNESS, ...
%!                              FEM_VEC_LOAD_CONSISTENT], ...
%!                             load_case);
%!
%! ## Effectively solve (-omega^2 * M + K) * U = R
%! mat_ass.K += -omega^2 * mat_ass.M;
%!
%! sol_harm = fem_sol_static(mesh, dof_map, mat_ass);
%!
%! sol_harm.stress = fem_ass_matrix(mesh, dof_map, [FEM_VEC_STRESS_CAUCH], load_case, sol_harm);
%!
%! tol_u = 0.003;
%! tol_tau = 0.08;
%! taun = nan(rows(mesh.nodes), 6);
%! idxtens = int32([1, 4, 6; 4, 2, 5; 6, 5, 3]);
%!
%! for i=1:rows(sol_harm.stress.taum.iso20)
%!   for j=1:20
%!     taun(mesh.elements.iso20(i, j), :) = sol_harm.stress.taum.iso20(i, j, :);
%!   endfor
%! endfor
%! opts.scale_def = 0.5 * a / max(max(abs(sol_harm.def)));
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! opts.skin_only = true;
%! opts.show_element = true;
%! if (f_run_post_proc)
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   opts.rotation_angle = [0, 0, 0];
%!   fem_post_sol_external(mesh, sol_harm, opts);
%!   [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!   figure("visible", "off");
%!   imshow(img, map);
%!   title("Gmsh - deformed mesh / continuous stress tensor");
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     unlink([opts.print_to_file, "_001.jpg"]);
%!   endif
%! end_unwind_protect
%! figure_list();
%! endif
%! for i=[1, size(mesh.structured.inode_idx, 1)]
%!   switch (i)
%!   case 1
%!     ur = ur_a;
%!     taur = taur_a;
%!     taut = taut_a;
%!     tauz = tauz_a;
%!   case size(mesh.structured.inode_idx, 1)
%!     ur = ur_b;
%!     taur = taur_b;
%!     taut = taut_b;
%!     tauz = tauz_b;
%!   endswitch
%!
%!   for j=1:size(mesh.structured.inode_idx, 2)
%!     for k=1:size(mesh.structured.inode_idx, 3)
%!       inode = mesh.structured.inode_idx(i, j, k);
%!       if (inode)
%!         e1 = [mesh.nodes(inode, 1:2).'; 0];
%!         e3 = [0; 0; 1];
%!         e2 = cross(e3, e1);
%!         R = [e1, e2, e3];
%!         R *= diag(1 ./ norm(R, "cols"));
%!         U = R.' * sol_harm.def(inode, 1:3).';
%!         tau = R.' * taun(inode, :)(idxtens) * R;
%!         assert_simple(U, [ur; 0; 0], tol_u * abs(ur));
%!         assert_simple(tau, diag([taur, taut, tauz]), tol_tau * abs(p));
%!       endif
%!     endfor
%!   endfor
%! endfor

%!test
%! ## TEST7
%! ## shaft with notch/shaft with step
%! close all;
%! geo_types = {"notch", "step"};
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.L0 = 0.5 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! geo.h0 = 0.4e-3;
%! geo.h = 2e-3;
%! geo.xg = [-0.5 * geo.L, -0.5 * geo.L0, -0.5 * geo.w, 0.5 * geo.w, 0.5 * geo.L0, 0.5 * geo.L];
%! geo.hg = [       geo.h,        geo.h0,       geo.h0,      geo.h0,        geo.h];
%! load_data.F = 120;
%! f_run_post_proc = false;
%! function [x, y, z] = notched_shaft_geo(geo, r, s, t, varargin)
%!   Phi = s;
%!   R = t * interp1(geo.grid.r, geo.R, r);
%!   x = interp1(geo.grid.r, geo.x, r);
%!   y = R * cos(Phi);
%!   z = R * sin(Phi);
%! endfunction
%!
%! function [F, locked] = notched_shaft_bound(r, s, t, geo, load_data, varargin)
%!   F = [];
%!   locked = false(1, 3);
%!   if (r == 1)
%!     locked(1) = true;
%!   endif
%!   if (s == 0)
%!     locked(3) = true;
%!   endif
%!   if (s == pi / 2)
%!     locked(2) = true;
%!   endif
%! endfunction
%!
%! function p = notched_shaft_pressure(r, s, t, geo, load_data, perm_idx, varargin)
%!   p = [];
%!   if (r == numel(geo.x))
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%! for j=1:numel(geo_types)
%!   switch(geo_types{j})
%!   case "notch"
%!     A = 0.22;
%!     B = 1.37;
%!   case "step"
%!     A = 0.62;
%!     B = 3.5;
%!   endswitch
%!
%!   Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%!   tauxx_n = load_data.F / (geo.d^2 * pi / 4);
%!   tauxx_a = tauxx_n * Kt_a;
%!   geo.x = [];
%!   for i=1:numel(geo.xg) - 1
%!     geo.x = [geo.x(1:end - 1), linspace(geo.xg(i), geo.xg(i + 1), ceil((geo.xg(i + 1) - geo.xg(i)) / geo.hg(i)))];
%!   endfor
%!   geo.R = repmat(0.5 * geo.D, 1, numel(geo.x));
%!   idx = find(abs(geo.x) < 0.5 * geo.w);
%!   geo.R(idx) = 0.5 * geo.D - (sqrt(geo.r^2 - geo.x(idx).^2) - (geo.r - geo.t));
%!   switch (geo_types{j})
%!   case "step"
%!     geo.R(find(geo.x > 0)) = 0.5 * geo.d;
%!   endswitch
%!
%!   geo.grid.r = geometry.mesh_size.r = 1:numel(geo.R);
%!   geo.grid.s = geometry.mesh_size.s = linspace(pi / 2, 0, ceil(0.5 * geo.d * pi / 2 / geo.h0));
%!   geo.grid.t = geometry.mesh_size.t = linspace(0, 1, ceil(0.5 * geo.d / geo.h0));
%!   load_data.pressure = load_data.F / (geo.R(end)^2 * pi);
%!   geometry.sewing.tolerance = sqrt(eps) * 0.5 * geo.D;
%!   geometry.spatial_coordinates = @(r, s, t, varargin) feval("notched_shaft_geo", geo, r, s, t, varargin);
%!   geometry.material_selector = @(r, s, t, varargin) int32(1);
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("notched_shaft_bound", r, s, t, geo, load_data, varargin);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("notched_shaft_pressure", r, s, t, geo, load_data, perm_idx, varargin);
%!   options.elem_type = "iso8";
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%!   [dof_map] = fem_ass_dof_map(mesh, load_case);
%!
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_STIFFNESS, ...
%!                                FEM_VEC_LOAD_CONSISTENT], ...
%!                               load_case);
%!   options.number_of_threads = int32(2);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, options);
%!
%!   sol_stat.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.iso8)) / tauxx_n;
%!   fprintf(stdout, "geometry: %s\n", geo_types{j});
%!   fprintf(stdout, "\tKt=%.2f\n", Kt);
%!   fprintf(stdout, "\tKt_a=%.2f\n", Kt_a);
%!   assert_simple(Kt, Kt_a, 0.26 * Kt_a);
%!   opts.scale_def = 0.25 * geo.L / max(max(abs(sol_stat.def)));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
%!   if (f_run_post_proc)
%!   unwind_protect
%!     opts.print_to_file = tempname();
%!     opts.rotation_angle = [-pi/2, 0, 0];
%!     fem_post_sol_external(mesh, sol_stat, opts);
%!     [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("deformed mesh %s / Van Mises stress", geo_types{j}));
%!   unwind_protect_cleanup
%!     if (numel(opts.print_to_file))
%!       unlink([opts.print_to_file, "_001.jpg"]);
%!     endif
%!   end_unwind_protect
%!   endif
%! endfor
%! if (f_run_post_proc)
%! figure_list();
%! endif

%!test
%! ## TEST8
%! ## shaft with notch/shaft with step
%! close all;
%! geo_types = {"notch", "step"};
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.L0 = 0.5 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! geo.h0 = 0.8e-3;
%! geo.h = 4e-3;
%! geo.xg = [-0.5 * geo.L, -0.5 * geo.L0, -0.5 * geo.w, 0.5 * geo.w, 0.5 * geo.L0, 0.5 * geo.L];
%! geo.hg = 0.2 * [       geo.h,        geo.h0,       geo.h0,      geo.h0,        geo.h];
%! load_data.F = 120;
%! f_run_post_proc = false;
%! function [x, y, z] = notched_shaft_geo(geo, r, s, t, varargin)
%!   Phi = s;
%!   x = interp1(geo.grid.r, geo.x, r, "pchip");
%!   R = t * interp1(geo.x, geo.R, x, "pchip");
%!   y = R * cos(Phi);
%!   z = R * sin(Phi);
%! endfunction
%!
%! function [F, locked] = notched_shaft_bound(r, s, t, geo, load_data, varargin)
%!   F = [];
%!   locked = false(1, 3);
%!   if (r == 1)
%!     locked(1) = true;
%!   endif
%!   if (s == 0)
%!     locked(3) = true;
%!   endif
%!   if (s == pi / 2)
%!     locked(2) = true;
%!   endif
%! endfunction
%!
%! function p = notched_shaft_pressure(r, s, t, geo, load_data, perm_idx, varargin)
%!   p = [];
%!   if (r == numel(geo.x))
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%! for j=1:numel(geo_types)
%!   switch(geo_types{j})
%!   case "notch"
%!     A = 0.22;
%!     B = 1.37;
%!   case "step"
%!     A = 0.62;
%!     B = 3.5;
%!   endswitch
%!
%!   Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%!   tauxx_n = load_data.F / (geo.d^2 * pi / 4);
%!   tauxx_a = tauxx_n * Kt_a;
%!   geo.x = [];
%!   for i=1:numel(geo.xg) - 1
%!     geo.x = [geo.x(1:end - 1), linspace(geo.xg(i), geo.xg(i + 1), ceil((geo.xg(i + 1) - geo.xg(i)) / geo.hg(i)))];
%!   endfor
%!   geo.R = repmat(0.5 * geo.D, 1, numel(geo.x));
%!   idx = find(abs(geo.x) < 0.5 * geo.w);
%!   geo.R(idx) = 0.5 * geo.D - (sqrt(geo.r^2 - geo.x(idx).^2) - (geo.r - geo.t));
%!   switch (geo_types{j})
%!   case "step"
%!     geo.R(find(geo.x > 0)) = 0.5 * geo.d;
%!   endswitch
%!
%!   geo.grid.r = geometry.mesh_size.r = 1:numel(geo.R);
%!   geo.grid.s = geometry.mesh_size.s = linspace(pi / 2, 0, ceil(0.5 * geo.d * pi / 2 / geo.h0));
%!   geo.grid.t = geometry.mesh_size.t = linspace(0, 1, ceil(0.5 * geo.d / geo.h0));
%!   load_data.pressure = load_data.F / (geo.R(end)^2 * pi);
%!   geometry.sewing.tolerance = sqrt(eps) * 0.5 * geo.D;
%!   geometry.spatial_coordinates = @(r, s, t, varargin) feval("notched_shaft_geo", geo, r, s, t, varargin);
%!   geometry.material_selector = @(r, s, t, varargin) int32(1);
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("notched_shaft_bound", r, s, t, geo, load_data, varargin);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("notched_shaft_pressure", r, s, t, geo, load_data, perm_idx, varargin);
%!   options.elem_type = "iso20";
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%!   [dof_map] = fem_ass_dof_map(mesh, load_case);
%!
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_STIFFNESS, ...
%!                                FEM_VEC_LOAD_CONSISTENT], ...
%!                               load_case);
%!   options.number_of_threads = int32(2);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, options);
%!
%!   sol_stat.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.iso20)) / tauxx_n;
%!   fprintf(stdout, "geometry: %s\n", geo_types{j});
%!   fprintf(stdout, "\tKt=%.2f\n", Kt);
%!   fprintf(stdout, "\tKt_a=%.2f\n", Kt_a);
%!   assert_simple(Kt, Kt_a, 0.26 * Kt_a);
%!   opts.scale_def = 0.25 * geo.L / max(max(abs(sol_stat.def)));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
%!   if (f_run_post_proc)
%!   unwind_protect
%!     opts.print_to_file = tempname();
%!     opts.rotation_angle = [-pi/2, 0, 0];
%!     fem_post_sol_external(mesh, sol_stat, opts);
%!     [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("deformed mesh %s / Van Mises stress", geo_types{j}));
%!   unwind_protect_cleanup
%!     if (numel(opts.print_to_file))
%!       unlink([opts.print_to_file, "_001.jpg"]);
%!     endif
%!   end_unwind_protect
%!   endif
%! endfor
%! if (f_run_post_proc)
%!   figure_list();
%! endif

%!test
%! ## TEST9
%! ## VIBRATIONS OF COMPLETE SPHERICAL SHELLS WITH IMPERFECTIONS
%! ## Thomas A. Duffey
%! ## Jason E. Pepin
%! ## Amy N. Robertson
%! ## Michael L. Steinzig
%! ## Internatial Modal Analysis Conference (IMAC-XXIII)
%! ## Orlando, Florida
%! ## January 31-February 3, 2005
%! ## Los Alamos
%! ## NATIONAL LABORATORY
%!
%! close all;
%! material.E = 28e6 * 6895;
%! material.nu = 0.28;
%! material.rho = 0.000751 * 4.4482 / (25.4e-3^4);
%! geo.D = 2 * 4.4688 * 25.4e-3;
%! geo.t = 0.0625 * 25.4e-3;
%! N = 39;
%! f_run_post_proc = false;
%! function [x, y, z] = sphere_geo(geo, r, s, t, varagrin)
%!   R = 0.5 * geo.D - r * geo.t;
%!   Zeta = s;
%!   Psi = t;
%!   x = R * cos(Zeta) * cos(Psi);
%!   y = R * cos(Zeta) * sin(Psi);
%!   z = R * sin(Zeta);
%! endfunction
%!
%! function [F, locked] = sphere_bound(r, s, t, geo, load_data, varargin)
%!   F = [];
%!   locked = false(1, 3);
%! endfunction
%!
%! function p = sphere_pressure(r, s, t, geo, load_data, perm_idx, varargin)
%!   p = [];
%!   if (r == 1)
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%!   fref = [5078, 6005, 6378, 6729];
%!
%!   geometry.mesh_size.r = linspace(-0.5, 0.5, 2);
%!   geometry.mesh_size.s = linspace(-90 * pi / 180, 90 * pi / 180, 36);
%!   geometry.mesh_size.t = linspace(-180 * pi / 180, 180 * pi / 180, 36);
%!   load_data.pressure = 1;
%!   geometry.sewing.tolerance = sqrt(eps) * 0.5 * geo.D;
%!   geometry.spatial_coordinates = @(r, s, t, varargin) feval("sphere_geo", geo, r, s, t);
%!   geometry.material_selector = @(r, s, t, varargin) int32(1);
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("sphere_bound", r, s, t, geo, load_data);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("sphere_pressure", r, s, t, geo, load_data, perm_idx);
%!   options.elem_type = "iso20";
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%!   [dof_map] = fem_ass_dof_map(mesh, load_case);
%!
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_MASS, ...
%!				  FEM_MAT_STIFFNESS, ...
%!                                FEM_VEC_LOAD_CONSISTENT], ...
%!                               load_case);
%!   options.number_of_threads = int32(2);
%!   shift = sqrt(eps) * max(abs(diag(mat_ass.K))) / max(abs(diag(mat_ass.M)));
%!   tol = eps^0.4;
%!   alg = "shift-invert";
%!   solver = "pastix";
%!   num_threads = int32(2);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, shift, tol, alg, solver, num_threads);
%!   sol_eig.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_eig);
%!   opts.scale_def = 0.25 * geo.D / max(max(max(abs(sol_eig.def))));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
%!   if (f_run_post_proc)
%!   unwind_protect
%!     opts.print_to_file = tempname();
%!     opts.rotation_angle = [-pi/2, 0, 0];
%!     fem_post_sol_external(mesh, sol_eig, opts);
%!     for i=7:N
%!       [img, map, alpha] = imread(sprintf("%s_%03d.jpg", opts.print_to_file, i));
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title(sprintf("mode %d: %.0fHz", i, sol_eig.f(i)));
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(opts.print_to_file))
%!       unlink([opts.print_to_file, "_001.jpg"]);
%!     endif
%!   end_unwind_protect
%!   figure_list();
%!   endif
%! assert_simple(sol_eig.f([7, 12, 19, 39]), fref, 5e-4 * max(fref));

%!test
%! ## TEST10
%! ## distorted sphere
%! close all;
%! material.E = 28e6 * 6895;
%! material.nu = 0.28;
%! material.rho = 0.000751 * 4.4482 / (25.4e-3^4);
%! geo.a = 2 * 4.4688 * 25.4e-3;
%! geo.b = 1 * 4.4688 * 25.4e-3;
%! geo.c = 1.5 * 4.4688 * 25.4e-3;
%! geo.t = 0.625 * 25.4e-3;
%! geo.Psii = [-180, -135, -90, -45, 0,  45, 90, 135, 180] * pi / 180;
%! geo.Zetai = [-90, -45, 0, 45, 90] * pi / 180;
%! geo.Psir = [1, 1, 1, 1.05, 1, 1.05,  1, 1;
%!             1, 1, 1, 1.10, 1, 1.10,  1, 1;
%!             1, 1, 1, 1.05, 1.2, 1.05,  1, 1];
%! geo.Psir = [geo.Psir, geo.Psir(:, 1)];
%! geo.Psir = [ones(1, columns(geo.Psir));
%!             geo.Psir;
%!             ones(1, columns(geo.Psir))];
%! geo.Psii = [geo.Psii(1:end - 1) - 2 * pi, geo.Psii, geo.Psii(2:end) + 2 * pi];
%! geo.Psir = [geo.Psir(:, 1:end - 1), geo.Psir, geo.Psir(:, 2:end)];
%! N = 10;
%! f_run_post_proc = false;
%! function X = sphere_shape(geo, Psi, Zeta, varargin)
%!   Psir = sqrt(1 - Zeta^2) * interp2(geo.Psii, geo.Zetai, geo.Psir, Psi, Zeta,  "pchip");
%!   if (~(isreal(Psir) && isfinite(Psir)))
%!    error("interpolation failed");
%!   endif
%!   X = zeros(3, 1);
%!   X(1) = Psir * geo.a * cos(Psi);
%!   X(2) = Psir * geo.b * sin(Psi);
%!   X(3) = geo.c * Zeta;
%! endfunction
%!
%! function [x, y, z] = sphere_geo(geo, r, s, t, varargin)
%!   Zeta = sin(t);
%!   Psi = s;
%!   dZeta = dPsi = 2 * pi * sqrt(eps);
%!   if (Zeta + dZeta > 1)
%!     dZeta = -dZeta;
%!   endif
%!   if (Psi + dPsi > pi)
%!     dPsi = -dPsi;
%!   endif
%!   P0 = sphere_shape(geo, Psi, Zeta);
%!   P1 = sphere_shape(geo, Psi, Zeta + dZeta);
%!   P2 = sphere_shape(geo, Psi + dPsi, Zeta);
%!   n1 = (P1 - P0) * sign(dZeta);
%!   n2 = (P2 - P0) * sign(dPsi);
%!   n3 = cross(n1, n2);
%!   if (norm(n3) == 0)
%!     if (t == -pi/2)
%!       n3 = [0; 0; 1];
%!     elseif (t == pi/2)
%!       n3 = [0; 0; -1];
%!     else
%!       error("surface normal vector is singular");
%!     endif
%!   endif
%!   n3 /= norm(n3);
%!   P = P0 - n3 * geo.t * r;
%!   if (~isfinite(P))
%!     error("P is not finite");
%!   endif
%!   x = P(1);
%!   y = P(2);
%!   z = P(3);
%! endfunction
%!
%! function [F, locked] = sphere_bound(r, s, t, geo, load_data, varargin)
%!   F = [];
%!   locked = false(1, 3);
%! endfunction
%!
%! function p = sphere_pressure(r, s, t, geo, load_data, perm_idx, varargin)
%!   p = [];
%!   if (r == 1)
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%! figure("visible", "off");
%! hold on;
%! Psi = linspace(-pi, pi, 180);
%! Zeta = linspace(-1, 1, 9);
%! for j=1:numel(Zeta)
%!   X = zeros(3, numel(Psi));
%!   for i=1:numel(Psi)
%!     X(:, i) = sphere_shape(geo, Psi(i), Zeta(j));
%!   endfor
%!   set(plot3(X(1, :), X(2, :), X(3, :)), "color", [1, 0, 0]);
%! endfor
%! Zeta = linspace(-1, 1, 50);
%! Psi = linspace(-pi, pi, 18);
%! for j=1:numel(Psi)
%!   X = zeros(3, numel(Zeta));
%!   for i=1:numel(Zeta)
%!     X(:, i) = sphere_shape(geo, Psi(j), Zeta(i));
%!   endfor
%!   set(plot3(X(1, :), X(2, :), X(3, :)), "color", [0, 0, 1]);
%! endfor
%! daspect(ones(1, 3));
%! title("wireframe");
%! xlabel("x");
%! ylabel("y");
%! zlabel("z");
%!   geometry.mesh_size.r = linspace(0, 1, 2);
%!   geometry.mesh_size.s = linspace(-pi, pi, 18);
%!   geometry.mesh_size.t = linspace(-pi/2, pi/2, 18);
%!   load_data.pressure = 1;
%!   geometry.sewing.tolerance = sqrt(eps) * geo.a;
%!   geometry.spatial_coordinates = @(r, s, t, varargin) feval("sphere_geo", geo, r, s, t, varargin);
%!   geometry.material_selector = @(r, s, t, varargin) int32(1);
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("sphere_bound", r, s, t, geo, load_data, varargin);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("sphere_pressure", r, s, t, geo, load_data, perm_idx, varargin);
%!   options.elem_type = "iso20";
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%!   [dof_map] = fem_ass_dof_map(mesh, load_case);
%!
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_MASS, ...
%!				  FEM_MAT_STIFFNESS, ...
%!                                FEM_VEC_LOAD_CONSISTENT], ...
%!                               load_case);
%!   shift = sqrt(eps) * max(abs(diag(mat_ass.K))) / max(abs(diag(mat_ass.M)));
%!   tol = eps^0.4;
%!   alg = "shift-invert";
%!   solver = "pastix";
%!   num_threads = int32(2);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, shift, tol, alg, solver, num_threads);
%!   sol_eig.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_eig);
%!   opts.scale_def = 0.5 * geo.a / max(max(max(abs(sol_eig.def))));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
%! if (f_run_post_proc)
%!   unwind_protect
%!     opts.print_to_file = tempname();
%!     opts.rotation_angle = [-pi/2, 0, 0];
%!     fem_post_sol_external(mesh, sol_eig, opts);
%!     for i=7:N
%!       [img, map, alpha] = imread(sprintf("%s_%03d.jpg", opts.print_to_file, i));
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title(sprintf("mode %d: %.0fHz", i, sol_eig.f(i)));
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(opts.print_to_file))
%!       unlink([opts.print_to_file, "_001.jpg"]);
%!     endif
%!   end_unwind_protect
%! figure_list();
%! endif

%!test
%! ## TEST11
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.number_of_modes = 10;
%! param.Udyn = eye(3) / SI_unit_meter;
%! param.fmin = 0 / (SI_unit_second^-1);
%! param.fmax = 15000 / (SI_unit_second^-1);
%! param.num_freq = 1000;
%! geometry(1).user_data.helspr.L = 25.8e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.d = 1.3e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.n = 5.3;
%! geometry(1).user_data.helspr.ni = 3;
%! geometry(1).user_data.helspr.ng = 0.75;
%! geometry(1).user_data.color = "y";
%! geometry(2).user_data.helspr.L = 27.7e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.d = 1.3e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.n = 5.7;
%! geometry(2).user_data.helspr.ni = 2.7;
%! geometry(2).user_data.helspr.ng = 0.75;
%! geometry(2).user_data.color = "g";
%! geometry(3).user_data.helspr.L = 28.63e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.d = 1.25e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.n = 6;
%! geometry(3).user_data.helspr.ni = 2.7;
%! geometry(3).user_data.helspr.ng = 0.75;
%! geometry(3).user_data.color = "b";
%! geometry = geometry(1);
%! material.E = 206000e6 / SI_unit_pascal;
%! material.G = 81500e6 / SI_unit_pascal;
%! material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! damp.D = [1e-2; 0.03e-2];
%! damp.f = [20, 15000] / (SI_unit_second^-1);
%! param.Fz = -7.2 / SI_unit_newton;
%! [material.alpha, material.beta] = fem_pre_mat_rayleigh_damping(damp.D, damp.f);
%! param.omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%!
%! function [x, y, z, R, Phi] = helspr_geo(geo, r, s, t, varargin)
%!   Phi = 2 * pi * (r * geo.n + geo.ni);
%!   r1 = 0.5 * geo.d * s;
%!   Theta = 2 * pi * t;
%!   x1 = [0.5 * geo.D * cos(Phi);
%!         0.5 * geo.D * sin(Phi);
%!         (geo.L - geo.d * (2 * (geo.ni - geo.ng) + 1)) * r + geo.d * (geo.ni - geo.ng + 0.5)];
%!   x2 = [0;
%!         0;
%!         x1(3)];
%!   e2 = x2 - x1;
%!   e1 = [-0.5 * geo.D * sin(Phi) * 2 * pi * geo.n;
%!          0.5 * geo.D * cos(Phi) * 2 * pi * geo.n;
%!          (geo.L - geo.d * (2 * (geo.ni - geo.ng) + 1))];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R1 = [e1, e2, e3];
%!   R1 *= diag(1 ./ norm(R1, "cols"));
%!   assert_simple(R1.' * R1, eye(3), eps^0.9);
%!   assert_simple(R1 * R1.', eye(3), eps^0.9);
%!   x3 = R1 * [0; r1 * cos(Theta); r1 * sin(Theta)] + x1;
%!   x = x3(1);
%!   y = x3(2);
%!   z = x3(3);
%! endfunction
%!
%! function [F, locked] = boundary_cond(r, s, t, geo, loads, varargin)
%!   locked = [];
%!   F = [];
%! endfunction
%!
%! function p = pressure_boundary_cond(r, s, t, geometry, load, perm_idx, varargin)
%!   p = [];
%! endfunction
%!
%! function [Freact, sol_stat, sol_eig, sol_eig_def] = transmissibility(geometry, material, param)
%!   material.nu = material.E / (2 * material.G) - 1;
%!   geometry.user_data.helspr.D = geometry.user_data.helspr.Di + geometry.user_data.helspr.d;
%!   kz = material.G * geometry.user_data.helspr.d^4 / (8 * geometry.user_data.helspr.n * geometry.user_data.helspr.D^3);
%!   Ustat = [0; 0; param.Fz / kz];
%!   h = geometry.user_data.helspr.d * pi / 4;
%!   geometry.user_data.helspr.nPhi = max([2, round(sqrt((geometry.user_data.helspr.D * pi * geometry.user_data.helspr.n)^2 + geometry.user_data.helspr.L^2) / h)]) + 1;
%!   geometry.user_data.helspr.nr = max([1, round(0.5 * geometry.user_data.helspr.d / h)]) + 1;
%!   geometry.user_data.helspr.nTheta = max([3, round(geometry.user_data.helspr.d * pi / h)]) + 1;
%!   geometry.mesh_size.r = linspace(0, 1, geometry.user_data.helspr.nPhi);
%!   geometry.mesh_size.s = linspace(0, 1, geometry.user_data.helspr.nr);
%!   geometry.mesh_size.t = linspace(0, 1, geometry.user_data.helspr.nTheta);
%!   geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.helspr.D;
%!   geometry.spatial_coordinates = @(r, s, t, geo, varargin) feval("helspr_geo", geometry.user_data.helspr, r, s, t, varargin);
%!   geometry.material_selector = @(r, s, t, geo, varargin) int32(1);
%!   geometry.boundary_condition = @(r, s, t, geo, load, varargin) feval("boundary_cond", r, s, t, geo, load, varargin);
%!   geometry.pressure_boundary_condition = @(r, s, t, geo, load, perm_idx, varargin) feval("pressure_boundary_cond", r, s, t, geo, load, perm_idx, varargin);
%!   options.elem_type = "iso20";
%!   loads = struct();
%!   [mesh, load_case_dof] = fem_pre_mesh_struct_create(geometry, loads, material, options);
%!   idx_node_bottom = unique(mesh.structured.inode_idx(1, :, :)(:));
%!   idx_node_top = unique(mesh.structured.inode_idx(end, :, :)(:));
%!   idx_node_bottom = idx_node_bottom(idx_node_bottom > 0);
%!   idx_node_top = idx_node_top(idx_node_top > 0);
%!   idx_node_joint = [idx_node_bottom; idx_node_top];
%!   empty_cell = cell(1, numel(idx_node_joint));
%!   mesh.elements.joints = struct("nodes", mat2cell(idx_node_joint, ones(numel(idx_node_joint), 1, "int32"), 1), "C", repmat({[eye(3), zeros(3, 3)]}, numel(idx_node_joint), 1));
%!   [dof_map] = fem_ass_dof_map(mesh, load_case_dof);
%!   load_case_stat = fem_pre_load_case_create_empty(1);
%!   for i=1:numel(load_case_stat)
%!     load_case_stat(i).joints = struct("U", repmat({zeros(3, 1)}, numel(idx_node_joint), 1));
%!     for j=1:numel(idx_node_top)
%!       load_case_stat(i).joints(numel(idx_node_bottom) + j).U = Ustat(:, i);
%!     endfor
%!   endfor
%!   load_case_dyn = fem_pre_load_case_create_empty(columns(param.Udyn));
%!   for i=1:numel(load_case_dyn)
%!     load_case_dyn(i).joints = struct("U", repmat({zeros(3, 1)}, numel(idx_node_joint), 1));
%!     for j=1:numel(idx_node_top)
%!       load_case_dyn(i).joints(numel(idx_node_bottom) + j).U = param.Udyn(:, i);
%!     endfor
%!   endfor
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS, ...
%!                                         FEM_MAT_STIFFNESS, ...
%!                                         FEM_VEC_LOAD_CONSISTENT], ...
%!                                        load_case_stat);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, param.number_of_modes);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case_stat, ...
%!                                    sol_stat);
%!   for i=1:numel(load_case_dyn)
%!     load_case_dyn(i).tau0 = sol_stat.stress.tau;
%!   endfor
%!   mesh_def = mesh;
%!   mesh_def.nodes += sol_stat.def;
%!   [mat_ass_def.M, ...
%!    mat_ass_def.D, ...
%!    mat_ass_def.K, ...
%!    mat_ass_def.KTAU0, ...
%!    mat_ass_def.R, ...
%!    mat_ass_def.mat_info, ...
%!    mat_ass_def.mesh_info] = fem_ass_matrix(mesh_def, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_TAU0, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_dyn);
%!   mat_ass_def.K += mat_ass_def.KTAU0;
%!   [sol_eig_def] = fem_sol_modal(mesh_def, dof_map, mat_ass_def, param.number_of_modes);
%!   Freact = complex(zeros(3, columns(mat_ass_def.R), numel(param.omega)));
%!   opt_factor.number_of_threads = int32(4);
%!   opt_factor.verbose = int32(0);
%!   for i=1:numel(param.omega)
%!     fprintf(stderr, "%3d:%.2fHz\n", i, param.omega(i) / (2 * pi));
%!     A = -param.omega(i)^2 * mat_ass_def.M + 1j * param.omega(i) * mat_ass_def.D + mat_ass_def.K;
%!     Uij = fem_sol_factor(A, opt_factor) \ mat_ass_def.R;
%!     for j=1:size(Freact, 2)
%!       Freact(:, j, i) = sum(Uij(:, j)(dof_map.edof.joints(1:numel(idx_node_bottom), :)), 1)(:) * mat_ass_def.mat_info.beta(3);
%!     endfor
%!   endfor
%! endfunction
%! Freact = cell(1, numel(geometry));
%! for i=1:numel(geometry)
%!   [Freact{i}] = transmissibility(geometry(i), material, param);
%! endfor
%! figure("visible", "off");
%! hold on;
%! for i=1:numel(geometry)
%!   plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * sqrt(Freact{i}(1, 1, :).^2 + Freact{i}(2, 2, :).^2 + Freact{i}(3, 3, :).^2)), sprintf("-;%d;%s", i, geometry(i).user_data.color));
%! endfor
%! xlabel("f [Hz]");
%! ylabel("kdyn [dB/(1N/m)]");
%! title("overall transmissibility");
%! grid minor on;
%! for j=1:columns(param.Udyn)
%!   figure("visible", "off");
%!   hold on;
%!   for i=1:numel(geometry)
%!     plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * abs(Freact{i}(i, i, :))), sprintf("-;%d;%s", i, geometry(i).user_data.color));
%!   endfor
%!   xlabel("f [Hz]");
%!   ylabel(sprintf("kdyn%s [dB/(1N/m)]", {"x", "y", "z"}{j}));
%!   title(sprintf("transmissibility %s-direction", {"x", "y", "z"}{j}));
%!   grid minor on;
%! endfor

%!test
%! ## TEST12
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.N = 200;
%! param.fmin = 0;
%! param.fmax = 15000 / (SI_unit_second^-1);
%! param.num_freq = 1000;
%! damp.D = [1e-2; 0.03e-2];
%! damp.f = [20, 15000] / (SI_unit_second^-1);
%! helspr1.L = 25.8e-3 / SI_unit_meter;
%! helspr1.Di = 12.12e-3 / SI_unit_meter;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.n = 5.3;
%! helspr1.ni = 3;
%! helspr1.ng = 0.75;
%! helspr1.m = 400;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! [helspr1.material.alpha, helspr1.material.beta] = fem_pre_mat_rayleigh_damping(damp.D, damp.f);
%! helspr1.material.nu = helspr1.material.E / (2 * helspr1.material.G) - 1;
%! helspr1.material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! helspr1.D = helspr1.Di + helspr1.d;
%! section1.A = helspr1.d^2 * pi / 4;
%! section1.Ay = 9 / 10 * section1.A;
%! section1.Az = section1.Ay;
%! section1.Iy = helspr1.d^4 * pi / 64;
%! section1.Iz = section1.Iy;
%! section1.It = section1.Iy + section1.Iz;
%! helspr1.Phi = linspace(0, 2 * pi * helspr1.n, ceil(helspr1.m * helspr1.n)).' + 2 * pi * helspr1.ni;
%! helspr1.x = 0.5 * helspr1.D * cos(helspr1.Phi);
%! helspr1.y = 0.5 * helspr1.D * sin(helspr1.Phi);
%! helspr1.z = (helspr1.L - helspr1.d * (2 * (helspr1.ni - helspr1.ng) + 1)) * linspace(0, 1, numel(helspr1.Phi))(:) + helspr1.d * (helspr1.ni - helspr1.ng + 0.5);
%! helspr1.e2 = [0, 0, 1];
%! elnodes = int32([(1:numel(helspr1.Phi) - 1).', (2:numel(helspr1.Phi)).']);
%! mesh.nodes = [[helspr1.x, helspr1.y, helspr1.z], zeros(numel(helspr1.Phi), 3)];
%! mesh.elements.beam2 = struct("nodes", mat2cell(elnodes, ones(numel(helspr1.Phi) - 1, 1, "int32"), 2), ...
%!                              "section", mat2cell(repmat(section1, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32")), ...
%!                              "e2", mat2cell(repmat(helspr1.e2, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32"), 3));
%! mesh.material_data(1).E = helspr1.material.E;
%! mesh.material_data(1).nu = helspr1.material.nu;
%! mesh.material_data(1).beta = helspr1.material.beta;
%! mesh.material_data(1).alpha = helspr1.material.alpha;
%! mesh.material_data(1).rho = helspr1.material.rho;
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.beam2 = repmat(int32(1), numel(helspr1.Phi) - 1, 1);
%! empty_cell = cell(1, 2);
%! mesh.elements.joints = struct("nodes", empty_cell, "C", empty_cell);
%! mesh.elements.joints(1).nodes = int32(1);
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.joints(2).nodes = rows(mesh.nodes);
%! mesh.elements.joints(2).C = eye(6);
%! Udyn = eye(6);
%! omega = linspace(0, 2 * pi * 15000, 10000);
%! load_case_dyn = fem_pre_load_case_create_empty(columns(Udyn));
%! for i=1:numel(load_case_dyn)
%!   load_case_dyn(i).joints = struct("U", repmat({zeros(6, 1)}, numel(mesh.elements.joints), 1));
%!   load_case_dyn(i).joints(2).U = Udyn(:, i);
%! endfor
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_DAMPING, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case_dyn);
%! opt_linsol.solver = "umfpack";
%! opt_linsol.pre_scaling = true;
%! opt_linsol.number_of_threads = int32(2);
%! opt_linsol.refine_max_iter = int32(50);
%! opt_linsol.epsilon_refinement = 1e-10;
%! opt_linsol.verbose = int32(0);
%! opt_eig.p = 5 * param.N;
%! opt_eig.maxit = int32(100);
%! omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%! [sol_eig, Phi, ~] = fem_sol_modal(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eig);
%! opt_post.elem_types = {"beam2"};
%! Freact = complex(zeros(columns(dof_map.edof.joints), columns(mat_ass.R), rows(dof_map.edof.joints), numel(omega)));
%! opt_factor.number_of_threads = int32(4);
%! opt_factor.verbose = int32(0);
%! for i=1:numel(omega)
%!   fprintf(stderr, "%3d:%.2fHz\n", i, omega(i) / (2 * pi));
%!   A = -omega(i)^2 * mat_ass.M + 1j * omega(i) * mat_ass.D + mat_ass.K;
%!   Uij = fem_sol_factor(A, opt_linsol) \ mat_ass.R;
%!   for k=1:size(Freact, 3)
%!     for j=1:size(Freact, 2)
%!       Freact(:, j, k, i) = Uij(dof_map.edof.joints(k, :), j) * mat_ass.mat_info.beta(2);
%!     endfor
%!   endfor
%! endfor
