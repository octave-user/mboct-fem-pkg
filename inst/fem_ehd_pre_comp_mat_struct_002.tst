## fem_ehd_pre_comp_mat_struct.m:02
%!test
%! try
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
