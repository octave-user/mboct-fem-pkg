## fem_pre_mesh_struct_create.m:04
%!test
%! try
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
%! options.elem_type = "iso20r";
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
