## fem_pre_mesh_struct_create.m:06
%!test
%! try
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
