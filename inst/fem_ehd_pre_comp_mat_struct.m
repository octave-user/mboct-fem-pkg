## Copyright (C) 2016(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{comp_mat} = fem_ehd_pre_comp_mat_struct(@var{bearing_dimensions}, @var{material}, @var{options})
## Compute a compliance matrix for use with MBDyn's module-hydrodynamic_plain_bearing2 for a simple structured mesh
##
## @var{bearing_dimensions}.bearing_diameter @dots{} Journal plain bearing diameter.
##
## @var{bearing_dimensions}.bearing_width @dots{} Journal plain bearing width.
##
## @var{bearing_dimensions}.bearing_geometry @dots{} Callback function for building mesh and boundary conditions.
##
## @var{material} @dots{} Material data for the Finite Element Mesh.
##
## @var{options}.verbose @dots{} Enable verbose output.
##
## @var{options}.output_matrices @dots{} Return all Finite Element matrices in @var{comp_mat}.
##
## @var{options}.output_solution @dots{} Output the Finite Element solution in @var{comp_mat}.
##
## @var{options}.number_of_modes @dots{} Number of mode shapes used for a mode based EHD method.
##
## @var{options}.reference_pressure @dots{} Nodal pressure values are scaled by @var{options}.reference_pressure.
##
## @var{options}.solver.number_of_threads @dots{} Number of threads to use for a multithreaded linear solver.
##
## @var{options}.solver.refine_max_iter @dots{} Maximum number of refinement iterations for the linear solver.
##
## @seealso{fem_ehd_pre_comp_mat_export, fem_ehd_pre_comp_mat_plot}
## @end deftypefn

function [comp_mat] = fem_ehd_pre_comp_mat_struct(bearing_dimensions, material, options)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  if (~isfield(bearing_dimensions, "bearing_geometry"))
    error("missing field bearing_dimensions.bearing_geometry");
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  if (~isfield(options, "output_matrices"))
    options.output_matrices = false;
  endif

  if (~isfield(options, "output_solution"))
    options.output_solution = false;
  endif

  if (~isfield(options, "number_of_modes"))
    options.number_of_modes = -1;
  endif

  if (~isfield(options, "solver"))
    options.solver = struct();
  endif
  
  if (options.verbose)
    fprintf(stderr, "creating 2D mesh ...\n");
    tic();
  endif

  [geometry, loads, comp_mat.bearing_surf] = feval(bearing_dimensions.bearing_geometry, bearing_dimensions, options);

  if (options.verbose)
    toc();
    fprintf(stderr, "creating 3D mesh ...\n");
    tic();
  endif

  options_mesh.verbose = options.verbose;
  options_mesh.identical_boundary_cond = true;

  [mesh, load_cases] = fem_pre_mesh_struct_create(geometry, loads, material, options_mesh);

  if (options.verbose)
    toc();
    fprintf(stderr, "assembling matrices ...\n");
    tic();
  endif

  dof_map = fem_ass_dof_map(mesh, load_cases(1));
  
  [mat_ass.K, ...
   mat_ass.R, ...
   mat_ass.Rlump] = fem_ass_matrix(mesh, ...
                                   dof_map, ...
                                   [FEM_MAT_STIFFNESS, ...
                                    FEM_VEC_LOAD_CONSISTENT, ...
                                    FEM_VEC_LOAD_LUMPED], ...
                                   load_cases);

  if (options.verbose)
    toc();
    fprintf(stderr, "building compliance matrix ...\n");
    tic();
  endif

  comp_mat.bearing_dimensions = bearing_dimensions;

  if (options.output_matrices)
    comp_mat.load_cases = load_cases;
    comp_mat.mat_ass = mat_ass;
    comp_mat.loads = loads;
  endif

  comp_mat.reference_pressure = options.reference_pressure;
  comp_mat.mesh = mesh;
  comp_mat.dof_map = dof_map;
  comp_mat.dX = zeros(3, 1);
  comp_mat.dR = eye(3);
  N = length(comp_mat.bearing_surf.s) * length(comp_mat.bearing_surf.t);

  if (options.number_of_modes <= 0)
    sol_stat = fem_sol_static(mesh, dof_map, mat_ass, options.solver);

    if (options.output_solution)
      comp_mat.sol_stat = sol_stat;
    endif

    comp_mat.C = zeros(N, N);

    for i=1:length(comp_mat.bearing_surf.s)
      for j=1:length(comp_mat.bearing_surf.t)
        idxrow = (i - 1) * length(comp_mat.bearing_surf.t) + j;
        inode = mesh.structured.inode_idx(1, i, j + comp_mat.bearing_surf.offset_t);
        for k=1:length(loads)
          idxcol = (loads(k).index.i - 1) * length(comp_mat.bearing_surf.t) + loads(k).index.j;
          U = sol_stat.def(inode, 1:2, k);
          n = mesh.nodes(inode, 1:2).';
          comp_mat.C(idxrow, idxcol) = U * (n / norm(n));
        endfor
      endfor
    endfor
  else
    if (options.verbose)
      fprintf(stderr, "building condensed matrices ...\n");
      tic();
    endif

    master_dof = zeros(1, N * 3, "int32");
    k = 0;

    for i=1:length(comp_mat.bearing_surf.s)
      for j=1:length(comp_mat.bearing_surf.t)
        inode = mesh.structured.inode_idx(1, i, j + comp_mat.bearing_surf.offset_t);
        node_dofs = dof_map.ndof(inode, 1:3);
        master_dof(k + (1:length(node_dofs))) = node_dofs;
        k += length(node_dofs);
      endfor
    endfor

    master_dof = unique(master_dof);

    [Kred, Tred] = fem_cms_red_guyan(mat_ass.K, master_dof);
    Tred = full(Tred);
    Kred = full(Kred);
    Rred = full(Tred.' * mat_ass.R);

    An = zeros(3, rows(mesh.nodes));

    for k=1:(length(load_cases) - length(comp_mat.bearing_surf.t)) ## last t load cases are redundant
      if (isfield(load_cases(k).pressure, "iso4"))
        for l=1:rows(load_cases(k).pressure.iso4.elements)
          for m=1:columns(load_cases(k).pressure.iso4.elements)
            p = load_cases(k).pressure.iso4.p(l, m);
            if (p ~= 0)
              inode = load_cases(k).pressure.iso4.elements(l, m);
              node_dof = dof_map.ndof(inode, 1:3);
              for i=1:length(node_dof)
                An(i, inode) += mat_ass.Rlump(node_dof(i), k) / p;
              endfor
            endif
          endfor
        endfor
      endif
    endfor

    A = zeros(1, length(mat_ass.K));

    for i=1:rows(mesh.nodes)
      node_dof = dof_map.ndof(i, 1:3);
      for j=1:length(node_dof)
        if (node_dof(j) > 0)
          A(node_dof(j)) = norm(An(:, i)); ## Why not using An(j, i) ??
        endif
      endfor
    endfor

    comp_mat.An = An;
    comp_mat.A = A;

    A = diag(A);
    Ared = Tred.' * A * Tred;

    if (options.verbose)
      toc();
      fprintf(stderr, "computing invariants of condensed stiffness matrix ...\n");
      tic();
    endif

    [Phi, lambda] = eig_stiff_mat(Kred, Ared, options.number_of_modes);

    if (options.verbose)
      toc();
      fprintf(stderr, "computing mode shapes ...\n");
      tic();
    endif

    n = zeros(length(comp_mat.bearing_surf.s) * length(comp_mat.bearing_surf.t), length(master_dof));

    for i=1:length(comp_mat.bearing_surf.s)
      for j=1:length(comp_mat.bearing_surf.t)
        inode = mesh.structured.inode_idx(1, i, j + comp_mat.bearing_surf.offset_t);
        node_dofs = dof_map.ndof(inode, 1:3);
        nij = mesh.nodes(inode, 1:2).';
        nij /= norm(nij);
        idxrow = (i - 1) * length(comp_mat.bearing_surf.t) + j;

        for k=1:length(nij)
          idxcol = find(master_dof == node_dofs(k));
          n(idxrow, idxcol) = nij(k);
        endfor
      endfor
    endfor

    kappa = zeros(columns(Phi), 1);

    for i=1:columns(Phi)
      kappa(i) = sqrt(1 / (Phi(:, i).' * Kred * Phi(:, i)));
    endfor

    Phi *= diag(kappa);

    KPhi = Phi.' * Kred * Phi;

    RPhi = Phi.' * Rred;
    Phin = n * Phi;

    if (isfield(options, "modal_error"))
      err_red = zeros(columns(Phi), 1);
      uref = 0;

      pred = zeros(rows(Rred), columns(Rred));
      ured = zeros(rows(Rred), columns(Rred));

      for i=1:columns(Rred)
        pred(:, i) = test_press_for_modal(bearing_dimensions, comp_mat.bearing_surf, loads, Rred, i);
        ured(:, i) = Kred \ pred(:, i);
        uref = max(uref, max(abs(n * ured(:, i))));
      endfor

      for j=1:columns(Phi)
        for i=1:columns(Rred)
          Kred2 = Phi(:, 1:j).' * Kred * Phi(:, 1:j);
          ured2 = Phi(:, 1:j) * (Kred2 \ (Phi(:, 1:j).' * pred(:, i)));
          err_red(j) += max(abs(n * (ured(:, i) - ured2))) / uref;
        endfor
        err_red(j) /= columns(Rred);
      endfor

      comp_mat.err_red = err_red;
    endif

    comp_mat.master_dof = master_dof;
    comp_mat.KPhi = KPhi;
    comp_mat.Phin = Phin;
    comp_mat.RPhi = RPhi;

    if (options.output_solution)
      U = Tred * Phi;

      for i=1:columns(U)
        U(:, i) /= max(abs(U(:, i)));
      endfor

      comp_mat.modal.def = fem_post_def_nodal(mesh, dof_map, U);

      comp_mat.modal_sol.def = fem_post_def_nodal(comp_mat.mesh, ...
                                                  comp_mat.dof_map, ...
                                                  Tred * Phi * (KPhi \ RPhi));
    endif

    comp_mat.lambda = lambda;

    if (options.verbose)
      toc();
    endif
  endif
endfunction

function [x, z] = surf_idx_to_surf_pos(bearing_dimensions, bearing_surf, s, t)
  x = bearing_dimensions.bearing_diameter * pi * (s - bearing_surf.s(1)) ...
      / (bearing_surf.s(end) - bearing_surf.s(1));

  z = bearing_dimensions.bearing_width * ((t - bearing_surf.t(1)) ...
                                          / (bearing_surf.t(end) - bearing_surf.t(1)) - 0.5);
endfunction

function pred = test_press_for_modal(bearing_dimensions, bearing_surf, loads, Rred, i)
  pred = zeros(rows(Rred), 1);

  [xi, zi] = surf_idx_to_surf_pos(bearing_dimensions, ...
                                  bearing_surf, ...
                                  loads(i).position.s, ...
                                  loads(i).position.t);

  for k=1:columns(Rred)
    [xk, zk] = surf_idx_to_surf_pos(bearing_dimensions, ...
                                    bearing_surf, ...
                                    loads(k).position.s, ...
                                    loads(k).position.t);

    if (abs(xk - xi) > 0.5 * pi * bearing_dimensions.bearing_diameter)
      xk -= sign(xk - xi) * bearing_dimensions.bearing_diameter * pi;
    endif

    dr = sqrt((xk - xi)^2 + (zk - zi)^2) / min(bearing_dimensions.bearing_width, ...
                                               pi * bearing_dimensions.bearing_diameter);
    if (dr < 0.25)
      pred += Rred(:, k);
    endif
  endfor
endfunction

function [Phi, lambda] = eig_stiff_mat(K, A, max_modes)
  opts.maxiter = int32(50000);
  opts.disp = int32(0);

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);

    [Phi, lambda, status] = eigs(K, A, max_modes, "sm", opts);

    if (status ~= 0)
      error("eigs failed with status %d", status);
    endif
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  [lambda, idx_lambda] = sort(diag(lambda), "ascend");

  max_modes = min(length(lambda), max_modes);

  lambda = lambda(1:max_modes);
  idx_lambda = idx_lambda(1:max_modes);
  Phi = Phi(:, idx_lambda);
endfunction

%!function [geometry, loads, bearing_surf] = bearing_callback(bearing_dimensions, options)
%!   dx = options.element_size;
%!   b = bearing_dimensions.bearing_width;
%!   d1 = bearing_dimensions.bearing_diameter;
%!   d2 = bearing_dimensions.recess_inner_diameter;
%!   d3 = bearing_dimensions.recess_outer_diameter;
%!   d4 = bearing_dimensions.outer_diameter;
%!   h0 = bearing_dimensions.bearing_offset;
%!   h1 = bearing_dimensions.total_height;
%!   h2 = bearing_dimensions.recess_offset;
%!   h3 = bearing_dimensions.recess_height;  
%!   l1 = h1 - b - h0;
%!   z1 = linspace(-(h1 - b / 2 - h0), -b / 2, max(2, round(l1 / dx) + 1));
%!   z2 = linspace(z1(end), b / 2, max(2, 2 * ceil(round(b / dx) / 2) + 1));
%!   z3 = linspace(z2(end), b / 2 + h0, max(2, round(h0 / dx) + 1));
%!   z4e = b / 2 + h0 - h3;
%!   z5e = b / 2;
%!   dz = abs(z2 - z4e);
%!   iz4e = find(dz == min(dz));  
%!   z4 = linspace(z1(end), z4e, iz4e(1));
%!   z5 = linspace(z4(end), z3(1), length(z2) - length(z4) + 1);
%!   z6 = linspace(z4(end), b / 2 + h0 - h2, length(z5));
%!   z7 = linspace(z6(end), z3(end), length(z3));
%!   z45 = [z4, z5(2:end)];
%!   x1 = linspace(0.5 * d1, 0.5 * d2, max(2, round(0.5 * (d2 - d1) / dx) + 1));
%!   x2 = linspace(x1(end), 0.5 * d3, max(2, round(0.5 * (d3 - d2) / dx) + 1));
%!   x3 = linspace(x2(end), 0.5 * d4, max(2, round(0.5 * (d4 - d3) / dx) + 1));
%!   z123 = [z1, z2(2:end), z3(2:end)];
%!   z1453 = [z1, z4(2:end), z5(2:end), z3(2:end)];
%!   z1467 = [z1, z4(2:end), z6(2:end), z7(2:end)];
%!   x = repmat([x1, x2(2:end), x3(2:end)].', 1, length(z123));
%!   z = zeros(rows(x), columns(x));
%!   for i=1:length(x1)
%!     for j=1:length(z123)
%!       z(i, j) = interp1(x1([1, end]), [z123(j), z1453(j)], x1(i), 'linear');
%!     endfor
%!   endfor
%!   for i=2:length(x2)
%!     for j=1:length(z1453)
%!       z(length(x1) - 1 + i, j) = interp1(x2([1, end]), [z1453(j), z1467(j)], x2(i), 'linear');
%!     endfor
%!   endfor
%!   for i=2:length(x3)
%!     z(length(x1) + length(x2) - 2 + i, 1:end) = z1467;
%!   endfor
%!   matid = zeros(rows(x) - 1, columns(x) - 1, "int32");
%!   for i=1:length(x1) - 1
%!     matid(i, :) = 1;
%!   endfor
%!   for i=1:length(x2)
%!     for j=1:length(z1)+length(z4) - 2
%!       matid(i + length(x1) - 2, j) = 1;
%!     endfor
%!   endfor
%!   for i=1:length(x3)
%!     for j=1:columns(x) - length(z7)
%!       matid(i + length(x1) + length(x2) - 2, j) = 1;
%!     endfor
%!   endfor
%!   geometry.grid.x = x;
%!   geometry.grid.z = z;
%!   geometry.grid.Phi = linspace(0, 2 * pi, max(2, 2 * ceil(round(d1 * pi / dx) / 2) + 1));
%!   geometry.grid.matid = matid;
%!   geometry.mesh_size.r = r = 1:rows(geometry.grid.x);
%!   geometry.mesh_size.s = s = 1:length(geometry.grid.Phi);
%!   geometry.mesh_size.t = t = 1:columns(geometry.grid.x);
%!   geometry.sewing.tolerance = sqrt(eps) * 0.5 * d1;
%!   geometry.spatial_coordinates = @bearing_callback_geometry;
%!   geometry.material_selector = @bearing_callback_material;
%!   geometry.boundary_condition = @bearing_callback_boundary_cond;
%!   geometry.pressure_boundary_condition = @bearing_callback_pressure;
%!   idx_t = length(z1):length(z1) + length(z2) - 1;
%!   bearing_surf.t = t(idx_t);
%!   bearing_surf.s = s;
%!   bearing_surf.offset_t = length(z1) - 1;
%!   bearing_surf.grid_x = 0.5 * d1 * geometry.grid.Phi;
%!   bearing_surf.grid_z = z(1, idx_t);
%!   loads.pressure = options.reference_pressure;
%!   loads.position.r = 1;
%!   loads.position.s = nan;
%!   loads.position.t = nan;
%!   loads.limits.t = [length(z1), length(z1) + length(z2) - 1];
%!   loads.index.i = int32(0);
%!   loads.index.j = int32(0);
%!   loads.idx_r = [];
%!   loads.idx_s = [];
%!   loads.idx_t = [];
%!   loads = repmat(loads, 1, length(bearing_surf.s) * length(bearing_surf.t));
%!   for i=1:length(bearing_surf.s)
%!     for j=1:length(bearing_surf.t)
%!       idx_st = (i - 1) * length(bearing_surf.t) + j;
%!       loads(idx_st).index.i = i;
%!       loads(idx_st).index.j = j;
%!       loads(idx_st).position.s = bearing_surf.s(i);
%!       loads(idx_st).position.t = bearing_surf.t(j);
%!       if idx_st > 1
%!         loads(idx_st).idx_r = int32(1);
%!         if i == 1
%!           loads(idx_st).idx_s = bearing_surf.s([i:i + 1, end - 1:end]);
%!         elseif i == length(bearing_surf.s)
%!           loads(idx_st).idx_s = bearing_surf.s([1:2, i - 1:i]);
%!         else
%!           loads(idx_st).idx_s = bearing_surf.s(i) - 1:bearing_surf.s(i) + 1;
%!         endif
%!         loads(idx_st).idx_t = bearing_surf.t(j) - 1:bearing_surf.t(j) + 1;
%!       endif
%!     endfor
%!   endfor
%!function [x, y, z, R, Phi] = bearing_callback_geometry(r, s, t, geometry)
%!   R = geometry.grid.x(r, t);
%!   Phi = geometry.grid.Phi(s);
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = geometry.grid.z(r, t);
%!function matid = bearing_callback_material(r, s, t, geometry)
%!   matid = geometry.grid.matid(r(1), t(1));
%!function [F, locked] = bearing_callback_boundary_cond(r, s, t, geometry, load)
%!  locked = [];
%!  F = [];
%!  if (t == 1 || r == rows(geometry.grid.x))
%!    locked = true(3, 1);
%!  endif
%!function p = bearing_callback_pressure(r, s, t, geometry, load_data, perm_idx)
%!   p = [];
%!   if (r == load_data.position.r)
%!     NPhi = length(geometry.grid.Phi);
%!     for i=1:length(s)
%!       if (s(i) == load_data.position.s ...
%!           || load_data.position.s == NPhi && s(i) == 1 ...
%!                                              || load_data.position.s == 1 && s(i) == NPhi)
%!         if (all(t >= load_data.limits.t(1) & t <= load_data.limits.t(end)))
%!           for j=1:length(t)
%!             if (t(j) == load_data.position.t)
%!               p = zeros(4, 1);
%!               p(perm_idx(i, j)) = load_data.pressure;
%!               return;
%!             endif
%!           endfor
%!         endif
%!       endif
%!     endfor
%!   endif
%!test
%! do_plot = false;
%! if (do_plot)
%! close all;
%! endif
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! bearing_dimensions.bearing_geometry = @bearing_callback;
%! bearing_dimensions.bearing_width = 14e-3;
%! bearing_dimensions.bearing_diameter = 16e-3;
%! bearing_dimensions.recess_inner_diameter = 18e-3;
%! bearing_dimensions.recess_outer_diameter = 32e-3;
%! bearing_dimensions.outer_diameter = 40e-3;
%! bearing_dimensions.bearing_offset = 0.5e-3;
%! bearing_dimensions.total_height = 22e-3;
%! bearing_dimensions.recess_offset = 0.6e-3;
%! bearing_dimensions.recess_height = 14e-3;
%! options.element_size = 8e-3;
%! options.reference_pressure = 1e9;
%! options.output_solution = true;
%! options.output_matrices = true;
%! options.number_of_modes = 2;
%! options.verbose = false;
%! comp_mat_file = "";
%! unwind_protect
%! comp_mat_file = tempname();
%! if (ispc())
%!   comp_mat_file(comp_mat_file == "\\") = "/";
%! endif
%! comp_mat = fem_ehd_pre_comp_mat_struct(bearing_dimensions, material, options);
%! fem_ehd_pre_comp_mat_export(comp_mat, struct(), comp_mat_file);
%! opt_plot.scale_modes = 1e-2;
%! opt_plot.plot_mesh = true;
%! opt_plot.modal_plot = true;
%! opt_plot.nodal_plot = true;
%! opt_plot.plot_load_cases = true;
%! opt_plot.plot_const_pressure = true;
%! opt_plot.contour_levels = 5;
%! if (do_plot)
%! fem_ehd_pre_comp_mat_plot(comp_mat, opt_plot);
%! endif
%! unwind_protect_cleanup
%!   if (numel(comp_mat_file))
%!     fn = dir([comp_mat_file, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! do_plot = true;
%! if (do_plot)
%!   close all;
%! endif
%! function [geometry, loads, bearing_surf] = bearing_callback(bearing_dimensions, options)
%!   dx = options.element_size;
%!   b = bearing_dimensions.bearing_width;
%!   d1 = bearing_dimensions.bearing_diameter;
%!   d2 = bearing_dimensions.recess_inner_diameter;
%!   d3 = bearing_dimensions.recess_outer_diameter;
%!   d4 = bearing_dimensions.outer_diameter;
%!   h0 = bearing_dimensions.bearing_offset;
%!   h1 = bearing_dimensions.total_height;
%!   h2 = bearing_dimensions.recess_offset;
%!   h3 = bearing_dimensions.recess_height;  
%!   l1 = h1 - b - h0;
%!   z1 = linspace(-(h1 - b / 2 - h0), -b / 2, max(2, round(l1 / dx) + 1));
%!   z2 = linspace(z1(end), b / 2, max(2, 2 * ceil(round(b / dx) / 2) + 1));
%!   z3 = linspace(z2(end), b / 2 + h0, max(2, round(h0 / dx) + 1));
%!   z4e = b / 2 + h0 - h3;
%!   z5e = b / 2;
%!   dz = abs(z2 - z4e);
%!   iz4e = find(dz == min(dz));  
%!   z4 = linspace(z1(end), z4e, iz4e(1));
%!   z5 = linspace(z4(end), z3(1), length(z2) - length(z4) + 1);
%!   z6 = linspace(z4(end), b / 2 + h0 - h2, length(z5));
%!   z7 = linspace(z6(end), z3(end), length(z3));
%!   z45 = [z4, z5(2:end)];
%!   x1 = linspace(0.5 * d1, 0.5 * d2, max(2, round(0.5 * (d2 - d1) / dx) + 1));
%!   x2 = linspace(x1(end), 0.5 * d3, max(2, round(0.5 * (d3 - d2) / dx) + 1));
%!   x3 = linspace(x2(end), 0.5 * d4, max(2, round(0.5 * (d4 - d3) / dx) + 1));
%!   z123 = [z1, z2(2:end), z3(2:end)];
%!   z1453 = [z1, z4(2:end), z5(2:end), z3(2:end)];
%!   z1467 = [z1, z4(2:end), z6(2:end), z7(2:end)];
%!   x = repmat([x1, x2(2:end), x3(2:end)].', 1, length(z123));
%!   z = zeros(rows(x), columns(x));
%!   for i=1:length(x1)
%!     for j=1:length(z123)
%!       z(i, j) = interp1(x1([1, end]), [z123(j), z1453(j)], x1(i), 'linear');
%!     endfor
%!   endfor
%!   for i=2:length(x2)
%!     for j=1:length(z1453)
%!       z(length(x1) - 1 + i, j) = interp1(x2([1, end]), [z1453(j), z1467(j)], x2(i), 'linear');
%!     endfor
%!   endfor
%!   for i=2:length(x3)
%!     z(length(x1) + length(x2) - 2 + i, 1:end) = z1467;
%!   endfor
%!   matid = zeros(rows(x) - 1, columns(x) - 1, "int32");
%!   for i=1:length(x1) - 1
%!     matid(i, :) = 1;
%!   endfor
%!   for i=1:length(x2)
%!     for j=1:length(z1)+length(z4) - 2
%!       matid(i + length(x1) - 2, j) = 1;
%!     endfor
%!   endfor
%!   for i=1:length(x3)
%!     for j=1:columns(x) - length(z7)
%!       matid(i + length(x1) + length(x2) - 2, j) = 1;
%!     endfor
%!   endfor
%!   geometry.grid.x = x;
%!   geometry.grid.z = z;
%!   geometry.grid.Phi = linspace(0, 2 * pi, max(2, 2 * ceil(round(d1 * pi / dx) / 2) + 1));
%!   geometry.grid.matid = matid;
%!   geometry.mesh_size.r = r = 1:rows(geometry.grid.x);
%!   geometry.mesh_size.s = s = 1:length(geometry.grid.Phi);
%!   geometry.mesh_size.t = t = 1:columns(geometry.grid.x);
%!   geometry.sewing.tolerance = sqrt(eps) * 0.5 * d1;
%!   geometry.spatial_coordinates = @(varargin) feval("bearing_callback_geometry", varargin{:});
%!   geometry.material_selector = @(varargin) feval("bearing_callback_material", varargin{:});
%!   geometry.boundary_condition = @(varargin) feval("bearing_callback_boundary_cond", varargin{:});
%!   geometry.pressure_boundary_condition = @(varargin) feval("bearing_callback_pressure", varargin{:});
%!   idx_t = length(z1):length(z1) + length(z2) - 1;
%!   bearing_surf.t = t(idx_t);
%!   bearing_surf.s = s;
%!   bearing_surf.offset_t = length(z1) - 1;
%!   bearing_surf.grid_x = 0.5 * d1 * geometry.grid.Phi;
%!   bearing_surf.grid_z = z(1, idx_t);
%!   loads.pressure = options.reference_pressure;
%!   loads.position.r = 1;
%!   loads.position.s = nan;
%!   loads.position.t = nan;
%!   loads.limits.t = [length(z1), length(z1) + length(z2) - 1];
%!   loads.index.i = int32(0);
%!   loads.index.j = int32(0);
%!   loads.idx_r = [];
%!   loads.idx_s = [];
%!   loads.idx_t = [];
%!   loads = repmat(loads, 1, length(bearing_surf.s) * length(bearing_surf.t));
%!   for i=1:length(bearing_surf.s)
%!     for j=1:length(bearing_surf.t)
%!       idx_st = (i - 1) * length(bearing_surf.t) + j;
%!       loads(idx_st).index.i = i;
%!       loads(idx_st).index.j = j;
%!       loads(idx_st).position.s = bearing_surf.s(i);
%!       loads(idx_st).position.t = bearing_surf.t(j);
%!       if (idx_st > 1)
%!         loads(idx_st).idx_r = int32(1);
%!         if i == 1
%!           loads(idx_st).idx_s = bearing_surf.s([i:i + 1, end - 1:end]);
%!         elseif i == length(bearing_surf.s)
%!           loads(idx_st).idx_s = bearing_surf.s([1:2, i - 1:i]);
%!         else
%!           loads(idx_st).idx_s = bearing_surf.s(i) - 1:bearing_surf.s(i) + 1;
%!         endif
%!         loads(idx_st).idx_t = bearing_surf.t(j) - 1:bearing_surf.t(j) + 1;
%!       endif
%!     endfor
%!   endfor
%! endfunction
%! function [x, y, z, R, Phi] = bearing_callback_geometry(r, s, t, geometry)
%!   R = geometry.grid.x(r, t);
%!   Phi = geometry.grid.Phi(s);
%!   x = R .* cos(Phi);
%!   y = R .* sin(Phi);
%!   z = geometry.grid.z(r, t);
%! endfunction
%! function matid = bearing_callback_material(r, s, t, geometry)
%!   matid = geometry.grid.matid(r(1), t(1));
%! endfunction
%! function [F, locked] = bearing_callback_boundary_cond(r, s, t, geometry, load)
%!  locked = [];
%!  F = [];
%!  if (t == 1 || r == rows(geometry.grid.x))
%!    locked = true(3, 1);
%!  endif
%! endfunction
%! function p = bearing_callback_pressure(r, s, t, geometry, load_data, perm_idx)
%!   p = [];
%!   if (r == load_data.position.r)
%!     NPhi = length(geometry.grid.Phi);
%!     for i=1:length(s)
%!       if (s(i) == load_data.position.s ...
%!           || load_data.position.s == NPhi && s(i) == 1 ...
%!                                              || load_data.position.s == 1 && s(i) == NPhi)
%!         if (all(t >= load_data.limits.t(1) & t <= load_data.limits.t(end)))
%!           for j=1:length(t)
%!             if (t(j) == load_data.position.t)
%!               p = zeros(4, 1);
%!               p(perm_idx(i, j)) = load_data.pressure;
%!               return;
%!             endif
%!           endfor
%!         endif
%!       endif
%!     endfor
%!   endif
%! endfunction
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! bearing_dimensions.bearing_geometry = @(varargin) feval("bearing_callback", varargin{:});
%! bearing_dimensions.bearing_width = 14e-3;
%! bearing_dimensions.bearing_diameter = 16e-3;
%! bearing_dimensions.recess_inner_diameter = 18e-3;
%! bearing_dimensions.recess_outer_diameter = 20e-3;
%! bearing_dimensions.outer_diameter = 22e-3;
%! bearing_dimensions.bearing_offset = 0.5e-3;
%! bearing_dimensions.total_height = 16e-3;
%! bearing_dimensions.recess_offset = 0.6e-3;
%! bearing_dimensions.recess_height = 12e-3;
%! options.element_size = 1e-3;
%! options.reference_pressure = 1e9;
%! options.output_solution = true;
%! options.output_matrices = true;
%! options.number_of_modes = 3;
%! options.verbose = false;
%! comp_mat_file = "";
%! unwind_protect
%! comp_mat_file = tempname();
%! if (ispc())
%!   comp_mat_file(comp_mat_file == "\\") = "/";
%! endif
%! options.output_file = comp_mat_file;
%! comp_mat = fem_ehd_pre_comp_mat_struct(bearing_dimensions, material, options);
%! fem_ehd_pre_comp_mat_export(comp_mat, struct(), comp_mat_file);
%! opt_plot.scale_modes = 2e-3;
%! opt_plot.plot_mesh = true;
%! opt_plot.modal_plot = true;
%! opt_plot.nodal_plot = false;
%! opt_plot.plot_load_cases = false;
%! opt_plot.plot_const_pressure = false;
%! if (do_plot)
%! fem_ehd_pre_comp_mat_plot(comp_mat, opt_plot);
%! fig = get(0, "children");
%! for i=1:numel(fig)
%!   ax = get(fig(i), "children");
%!   for j=1:numel(ax)
%!     view(ax(j), 30, 30);
%!   endfor
%! endfor
%! figure_list();
%! endif
%! unwind_protect_cleanup
%!   if (numel(comp_mat_file))
%!     fn = dir([comp_mat_file, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
