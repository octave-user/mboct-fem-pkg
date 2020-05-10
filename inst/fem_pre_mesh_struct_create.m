## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
    options.verbose = false;
    options.identical_boundary_cond = false;
  endif

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

  x = y = z = zeros(numel(r), numel(s), numel(t));

  for i=1:numel(r)
    for j=1:numel(s)
      for k=1:numel(t)
        [x(i, j, k), y(i, j, k), z(i, j, k)] = feval(geometry.spatial_coordinates, r(i), s(j), t(k), geometry);
      endfor
    endfor
  endfor

  if (options.verbose)
    toc();
    fprintf(stderr, "merging nodes ...\n");
    tic();
  endif

  inode_idx = zeros(numel(r), numel(s), numel(t), "int32");
  mesh.nodes = zeros(numel(r) * numel(s) * numel(t), 6);

  inode = int32(0);

  for i=1:numel(r)
    for j=1:numel(s)
      for k=1:numel(t)
        sewing_cond = sqrt((x(i, j, k) - x).^2 + (y(i, j, k) - y).^2 + (z(i, j, k) - z).^2) < geometry.sewing.tolerance;
        sewing_cond(i, j, k) = false;
        inode2 = find_sewing_node(i, j, k, sewing_cond, inode_idx);
        if inode2 == 0
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

  mesh.nodes = mesh.nodes(1:inode, :);

  if (options.verbose)
    toc();
    fprintf(stderr, "creating elements ...\n");
    tic();
  endif

  mesh.elements.iso8 = zeros((numel(r) - 1) * (numel(s) - 1) * (numel(t) - 1), 8, "int32");
  mesh.elements.rbe3 = struct("nodes", [], "weight", [])([]);
  mesh.material_data = struct("C", [], "rho", [])([]);

  for i=1:numel(material)
    mesh.material_data(i).C = fem_pre_mat_isotropic(material(i).E, material(i).nu);
    mesh.material_data(i).rho = material(i).rho;
  endfor

  mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");

  ielement = int32(0);

  for i=2:numel(r)
    for j=2:numel(s)
      for k=2:numel(t)
        imat_id = feval(geometry.material_selector, r(i - 1:i), s(j - 1:j), t(k - 1:k), geometry);

        if (imat_id ~= 0)
          if (imat_id < 1 || imat_id > numel(material))
            error("invalid material number %d returned for %d, %d, %d", imat_id, i, j, k);
          endif

          mesh.materials.iso8(++ielement) = imat_id;
          mesh.elements.iso8(ielement, :) = [inode_idx(i,     j,     k), ...
                                             inode_idx(i - 1, j,     k), ...
                                             inode_idx(i - 1, j - 1, k), ...
                                             inode_idx(i,     j - 1, k), ...
                                             inode_idx(i,     j,     k - 1), ...
                                             inode_idx(i - 1, j,     k - 1), ...
                                             inode_idx(i - 1, j - 1, k - 1), ...
                                             inode_idx(i,     j - 1, k - 1)];
        endif
      endfor
    endfor
  endfor

  mesh.elements.iso8 = mesh.elements.iso8(1:ielement, :);
  mesh.materials.iso8 = mesh.materials.iso8(1:ielement);

  nodes_in_use = false(rows(mesh.nodes), 1);

  for i=1:rows(mesh.nodes)
    if (numel(find(mesh.elements.iso8 == i)) > 0)
      nodes_in_use(i) = true;
    endif
  endfor

  if (options.verbose)
    toc();
    fprintf(stderr, "creating load cases ...\n");
    tic();
  endif

  if (options.identical_boundary_cond)
    load_case = create_load_case(mesh, geometry, inode_idx, r, s, t, loads(1), nodes_in_use);
    load_case = fem_pre_load_case_merge(load_case, repmat(setfield(load_case, "locked_dof", []), 1, numel(loads) - 1));
  else
    empty_cell = cell(1, numel(loads));
    load_case = struct("locked_dof", empty_cell, "loaded_nodes", empty_cell, "loads", empty_cell);
    for l=1:numel(loads)
      load_case(l) = create_load_case(mesh, geometry, inode_idx, r, s, t, loads(l), nodes_in_use);
    endfor
  endif

  if (options.verbose)
    toc();
    fprintf(stderr, "creating pressure loads ...\n");
    tic();
  endif

  pressure_perm_idx = int32([3, 4;
                             2, 1]);

  for l=1:numel(loads)
    if (isfield(geometry, "pressure_boundary_condition"))
      ielement_press = int32(0);
      load_case(l).pressure.iso4.elements = zeros(numel(r) * (numel(s) - 1) * (numel(t) - 1), 4, "int32");
      load_case(l).pressure.iso4.p = zeros(rows(load_case(l).pressure.iso4.elements), 4);

      if (isfield(loads(l), "idx_r") && numel(loads(l).idx_r) > 0)
        idx_r = loads(l).idx_r;
      else
        idx_r = int32(1:numel(r));
      endif

      if (isfield(loads(l), "idx_s") && numel(loads(l).idx_s) > 0)
        idx_s = loads(l).idx_s;
      else
        idx_s = int32(1:numel(s));
      endif

      if (isfield(loads(l), "idx_t") && numel(loads(l).idx_t) > 0)
        idx_t = loads(l).idx_t;
      else
        idx_t = int32(1:numel(t));
      endif

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
    endif
  endfor

  mesh.structured.r = r;
  mesh.structured.s = s;
  mesh.structured.t = t;
  mesh.structured.inode_idx = inode_idx;
endfunction

function inode = find_sewing_node(i, j, k, sewing_cond, inode_idx)
  inode = int32(0);

  for k2=1:size(sewing_cond, 3)
    [i2, j2] = find(sewing_cond(:, :, k2));
    if numel(i2) > 0
      inode = inode_idx(i2(1), j2(1), k2);
      return;
    endif
  endfor
endfunction

function load_case = create_load_case(mesh, geometry, inode_idx, r, s, t, load, nodes_in_use)
  load_case.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
  load_case.loads = zeros(rows(mesh.nodes), columns(mesh.nodes));
  load_case.loaded_nodes = zeros(rows(mesh.nodes), 1, "int32");

  iload = int32(0);

  if (isfield(load, "idx_r") && numel(load.idx_r) > 0)
    idx_r = load.idx_r;
  else
    idx_r = int32(1:numel(r));
  endif

  if (isfield(load, "idx_s") && numel(load.idx_s) > 0)
    idx_s = load.idx_s;
  else
    idx_s = int32(1:numel(s));
  endif

  if (isfield(load, "idx_t") && numel(load.idx_t) > 0)
    idx_t = load.idx_t;
  else
    idx_t = int32(1:numel(t));
  endif

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

  load_case.loads = load_case.loads(1:iload, :);
  load_case.loaded_nodes = load_case.loaded_nodes(1:iload);
endfunction

%!demo
%! close all;
%! scale_stat = 2e-3;
%! scale_eig = 2e-3;
%! tol = eps;
%! number_of_modes = 3;
%!
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! function [F, locked] = boundary_cond(r, s, t, geo, load)
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
%! function p = pressure_boundary_cond(r, s, t, geometry, load)
%!   if (load.flags.use_pressure_boundary_cond && r == 0)
%!     p = repmat(load.p, 4, 1);
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
%! geometry.spatial_coordinates = @(r, s, t) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t);
%! geometry.material_selector = @(r, s, t) 1;
%! geometry.boundary_condition = @(r, s, t, geo, load) feval("boundary_cond", r, s, t, geo, load);
%! geometry.pressure_boundary_condition = @(r, s, t, geo, load) feval("pressure_boundary_cond", r, s, t, geo, load);
%! material.E = 210000e6;
%! material.nu = 0.26;
%! material.rho = 7850;
%! load.p = 1e6;
%! load.flags.use_pressure_boundary_cond = true;
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load, material);
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
%!
%! figure_list();

%!demo
%! close all;
%! scale_eig = 0.15;
%! number_of_modes = 14;
%!
%! ## Code_Aster SDLS109 V2.3.109
%! Rm = 0.369;
%! t = 0.048;
%! L = 0.05;
%! E = 185000e6;
%! nu = 0.3;
%! rho = 7850;
%! h = t / 6;
%! fref = [zeros(1, 6), 210.55, 210.55, 587.92, 587.92, 205.89, 205.89, 588.88, 588.88];
%!
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t)
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
%! geometry.spatial_coordinates = @(r, s, t) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t);
%! geometry.material_selector = @(r, s, t) 1;
%! geometry.boundary_condition =  @(r, s, t) {[],[]}{:};
%! material.E = E;
%! material.nu = nu;
%! material.rho = rho;
%! load_data = struct();
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material);
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
%! tol = 0.5e-2;
%! assert(sol_eig.f(1:numel(fref)), sort(fref), tol * max(fref));

%!demo
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
%!
%! function [x, y, z, R, Phi] = cylinder_geo(geo, r, s, t)
%!   R = (geo.r2 - geo.r1) * r + geo.r1;
%!   Phi = (geo.Phi2 - geo.Phi1) * s + geo.Phi1;
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = (geo.z2 - geo.z1) * t + geo.z1;
%! endfunction
%!
%! function p = pressure_callback(r, s, t, geometry, load_data, perm_idx)
%!   p = [];
%!   if (r == 0)
%!     p(perm_idx) = load_data.pressure;
%!   endif
%! endfunction
%!
%! function [F, locked] = boundary_cond_callback(r, s, t, geometry, load_data)
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
%! geometry.spatial_coordinates = @(r, s, t) feval("cylinder_geo", geometry.user_data.cylinder, r, s, t);
%! geometry.material_selector = @(r, s, t) 1;
%! geometry.boundary_condition =  @(r, s, t, geometry, load_data) feval("boundary_cond_callback", r, s, t, geometry, load_data);
%! geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx) feval("pressure_callback", r, s, t, geometry, load_data, perm_idx);
%! material.E = E;
%! material.nu = nu;
%! material.rho = rho;
%! load_data.pressure = p;
%!
%! [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material);
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
%!
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
%!       assert(U, [ur; 0; 0], tol_u * abs(ur));
%!       assert(tau, diag([taur, taut, tauz]), tol_tau * abs(p));
%!     endfor
%!   endfor
%! endfor

%!demo
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
%!
%! function [x, y, z] = notched_shaft_geo(geo, r, s, t)
%!   Phi = s;
%!   R = t * geo.R(r);
%!   x = geo.x(r);
%!   y = R * cos(Phi);
%!   z = R * sin(Phi);
%! endfunction
%!
%! function [F, locked] = notched_shaft_bound(r, s, t, geo, load_data)
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
%! function p = notched_shaft_pressure(r, s, t, geo, load_data, perm_idx)
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
%!   geometry.mesh_size.r = 1:numel(geo.R);
%!   geometry.mesh_size.s = linspace(pi / 2, 0, ceil(0.5 * geo.d * pi / 2 / geo.h0));
%!   geometry.mesh_size.t = linspace(0, 1, ceil(0.5 * geo.d / geo.h0));
%!   load_data.pressure = load_data.F / (geo.R(end)^2 * pi);
%!   geometry.sewing.tolerance = sqrt(eps) * 0.5 * geo.D;
%!   geometry.spatial_coordinates = @(r, s, t) feval("notched_shaft_geo", geo, r, s, t);
%!   geometry.material_selector = @(r, s, t) 1;
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data) feval("notched_shaft_bound", r, s, t, geo, load_data);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx) feval("notched_shaft_pressure", r, s, t, geo, load_data, perm_idx);
%!
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material);
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
%!   options.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, options);
%!
%!   sol_stat.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.iso8)) / tauxx_n;
%!   fprintf(stdout, "geometry: %s\n", geo_types{j});
%!   fprintf(stdout, "\tKt=%.2f\n", Kt);
%!   fprintf(stdout, "\tKt_a=%.2f\n", Kt_a);
%!   assert(Kt, Kt_a, 0.26 * Kt_a);
%!   opts.scale_def = 0.25 * geo.L / max(max(abs(sol_stat.def)));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
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
%! endfor
%! figure_list();