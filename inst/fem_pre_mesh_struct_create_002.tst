## fem_pre_mesh_struct_create.m:02
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
