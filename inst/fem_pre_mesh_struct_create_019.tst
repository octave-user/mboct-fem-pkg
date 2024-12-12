## fem_pre_mesh_struct_create.m:09
%!test
%! try
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
%!   options.elem_type = "iso20r";
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%!   [mesh] = fem_pre_mesh_solid_to_surf(mesh);
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
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
%!   shift = sqrt(eps) * max(abs(diag(mat_ass.K))) / max(abs(diag(mat_ass.M)));
%!   tol = eps^0.4;
%!   alg = "shift-invert";
%!   solver = "pastix";
%!   num_threads = mbdyn_solver_num_threads_default();
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
