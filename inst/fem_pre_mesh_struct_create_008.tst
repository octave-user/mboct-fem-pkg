## fem_pre_mesh_struct_create.m:08
%!test
%! try
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
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
