## fem_tests.m:50
%!test
%! try
%! ## TEST 50
%! do_plot = false;
%! close all;
%! elem_types = {"iso8", "iso20"};
%! for k=1:numel(elem_types)
%! options.elem_type = elem_types{k};
%! switch (options.elem_type)
%! case "iso8"
%!   M = 2;
%! case "iso20"
%!   M = 1;
%! endswitch
%! geometry.user_data.l = 10e-3;
%! geometry.user_data.w = 0.02e-3;
%! geometry.user_data.h = 0.02e-3;
%! dx = 1e-3 / M;
%!
%! function [x, y, z, R, Phi] = cube_geo(geo, r, s, t, varargin)
%!   x = geo.l * r;
%!   y = geo.w * s;
%!   z = geo.h * t;
%! endfunction
%!
%! function p = pressure_callback(r, s, t, geometry, load_data, perm_idx, varargin)
%!   p = [];
%! endfunction
%!
%! function [F, locked] = boundary_cond_callback(r, s, t, geometry, load_data, varargin)
%!   F = [];
%!   locked = [];
%! endfunction
%!
%! geometry.mesh_size.r = linspace(0, 1, max([2, ceil(geometry.user_data.l / dx) + 1]));
%! geometry.mesh_size.s = linspace(0, 1, max([2, ceil(geometry.user_data.w / dx) + 1]));
%! geometry.mesh_size.t = linspace(0, 1, max([2, ceil(geometry.user_data.h / dx) + 1]));
%! geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.l;
%! geometry.spatial_coordinates = @(r, s, t, varargin) feval("cube_geo", geometry.user_data, r, s, t);
%! geometry.material_selector = @(r, s, t, varargin) 1;
%! geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("boundary_cond_callback", r, s, t, geometry, load_data);
%! geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("pressure_callback", r, s, t, geometry, load_data, perm_idx);
%! k = 50;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! material.k = diag([k, k, k]);
%! material.cp = 465;
%! load_data.pressure = 0;
%! [mesh] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%! load_case.locked_dof = false(rows(mesh.nodes), 1);
%! load_case.domain = FEM_DO_THERMAL;
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! empty_cell = cell(1, 2);
%! sol = struct("t", empty_cell, "theta", empty_cell);
%! for j=1:2
%!   switch (j)
%!     case 1
%!       R = eye(3);
%!     case 2
%!       e1 = [0.5; 0.2; 0.1];
%!       e2 = [0; 1; 0];
%!       e3 = cross(e1, e2);
%!       e2 = cross(e3, e1);
%!       R = [e1 / norm(e1), e2 / norm(e2), e3 / norm(e3)];
%!   endswitch
%!   mesh.nodes = [mesh.nodes(:, 1:3) * R.', mesh.nodes(:, 4:6) * R.'];
%!   [mat_ass.Kk, mat_ass.C] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_THERMAL_COND, ...
%!                                             FEM_MAT_HEAT_CAPACITY], ...
%!                                            load_case);
%!   qref = material.rho * material.cp * geometry.user_data.l * geometry.user_data.w * geometry.user_data.h;
%!   assert_simple(sum(sum(mat_ass.C)), qref, eps^0.5 * abs(qref));
%!   theta_x0 = 20;
%!   theta_xl = 100;
%!   theta0 = (R(:,1).' * mesh.nodes(:, 1:3).' - R(:, 1).' * mesh.nodes(1, 1:3).') / geometry.user_data.l * (theta_xl - theta_x0) + theta_x0;
%!   dt = 0.1;
%!   alpha = 0.5;
%!   sol(j).t = 0:dt:100;
%!   sol(j).theta = zeros(dof_map.totdof, numel(sol(j).t));
%!   sol(j).theta(:, 1) = theta0;
%!   A = (1 / dt) * mat_ass.C + alpha * mat_ass.Kk;
%!   opts.number_of_threads = mbdyn_solver_num_threads_default();
%!   opts.solver = "chol";
%!   Afact = fem_sol_factor(A, opts);
%!   for i=2:numel(sol(j).t)
%!     sol(j).theta(:, i) = Afact \ (mat_ass.C * (sol(j).theta(:, i - 1)) / dt - mat_ass.Kk * (sol(j).theta(:, i - 1) * (1 - alpha)));
%!   endfor
%!   assert_simple(sol(j).theta(:, end), repmat(0.5 * (theta_x0 + theta_xl), rows(sol(j).theta), 1), eps^0.5 * abs(theta_xl));
%! endfor
%! assert_simple(sol(2).theta, sol(1).theta, eps^0.5 * abs(theta_xl));
%! if (do_plot)
%!   figure("visible", "off");
%!   hold on;
%!   for j=1:2
%!     plot(sol(j).t, max(sol(j).theta, [], 1), sprintf("-;max(theta%d);1", j));
%!     plot(sol(j).t, min(sol(j).theta, [], 1), sprintf("-;min(theta%d);3", j));
%!   endfor
%!   xlabel("t [s]");
%!   ylabel("theta [degC]");
%!   grid on;
%!   grid minor on;
%!   title("transient thermal problem of a bar");
%! endif
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
