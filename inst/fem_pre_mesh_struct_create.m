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

