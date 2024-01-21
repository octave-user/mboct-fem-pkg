## fem_tests.m:122
%!test
%! ## DEMO3
%! close all;
%! function [x, y, z, R, Phi] = cooks_membrane_geo(r, s, t, geometry, varargin)
%!   x = r * geometry.b;
%!   y = r * geometry.c + (r * (geometry.a1 + geometry.c - geometry.a0) + geometry.a0 - r * geometry.c) * s;
%!   z = t * geometry.h;
%! endfunction
%! function [F, locked] = cooks_membrane_bound_cond(r, s, t, geo, load, varargin)
%!   if (r == 1)
%!     F = load.F;
%!   else
%!     F = [];
%!   endif
%!   if (r == 0)
%!     locked = true(3, 1);
%!   else
%!     locked = [false; false; true];
%!   endif
%! endfunction
%! geometry.user_data.a0 = 44e-3;
%! geometry.user_data.a1 = 16e-3;
%! geometry.user_data.b = 44e-3;
%! geometry.user_data.c = 44e-3;
%! geometry.user_data.h = 1e-3;
%! material.E = 208.5044e6;
%! material.nu = 0.3;
%! material.rho = 1000;
%! geometry.sewing.tolerance = 0;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cooks_membrane_geo", r, s, t, geometry.user_data, varargin);
%! geometry.material_selector = @(r, s, t, varargin) int32(1);
%! geometry.boundary_condition = @(r, s, t, geo, load, varargin) feval("cooks_membrane_bound_cond", r, s, t, geo, load, varargin);
%! N = 2:22;
%! wA = zeros(3, length(N));
%! for i=1:length(N)
%!   geometry.mesh_size.r = linspace(0, 1, N(i));
%!   geometry.mesh_size.s = linspace(0, 1, N(i));
%!   geometry.mesh_size.t = [0, 1];
%!   loads.F = [0; 250; 0] / (length(geometry.mesh_size.s) * length(geometry.mesh_size.t));
%!   [mesh, load_case] = fem_pre_mesh_struct_create(geometry, loads, material);
%!   [dof_map] = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!   wA(:, i) = mean(sol_stat.def(mesh.structured.inode_idx(end, end, :), 1:3).', 2);
%! endfor
%! ## reference solution:
%! wA_ref =[-19.038e-3;
%!          24.098e-3;
%!          0];
%! for i=1:2
%!   figure("visible", "off");
%!   hold on;
%!   plot(N, 1e3 * wA(i, :), "-;wA;1");
%!   plot(N([1,end]), 1e3 * wA_ref(i)([1,end]), "-;wA_r_e_f;0");
%!   xlabel("N [1]");
%!   ylabel("wA [mm]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("w%s", {"x","y"}{i}));
%! endfor
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat, 0.3);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("Cook's membrane deformed mesh");
%! figure_list();
