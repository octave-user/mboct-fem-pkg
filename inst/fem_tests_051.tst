## fem_tests.m:51
%!test
%! try
%! ## TEST 51
%! ## The 1-D Heat Equation
%! ## 18.303 Linear Partial Differential Equations
%! ## Matthew J. Hancock
%! ## Fall 2006
%! ## https://ocw.mit.edu/courses/mathematics/18-303-linear-partial-differential-equations-fall-2006/lecture-notes/heateqni.pdf
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%!
%! elem_types = {"iso8", "iso20"};
%! function [x, y, z, R, Phi] = cube_geo(geo, r, s, t, varargin)
%!   x = geo.l * r;
%!   y = geo.w * s;
%!   z = geo.h * t;
%! endfunction
%! function p = pressure_callback(r, s, t, geometry, load_data, perm_idx, varargin)
%!   p = [];
%! endfunction
%! function [F, locked] = boundary_cond_callback(r, s, t, geometry, load_data, varargin)
%!   F = [];
%!   locked = [];
%! endfunction
%! for k=1:numel(elem_types)
%!   options.elem_type = elem_types{k};
%!   switch (options.elem_type)
%!     case "iso8"
%!       M = 2;
%!     case "iso20"
%!       M = 1;
%!   endswitch
%!   geometry.user_data.l = 20e-3;
%!   geometry.user_data.w = 1e-3;
%!   geometry.user_data.h = 1e-3;
%!   dx = 1e-3 / M;
%!   geometry.mesh_size.r = linspace(0, 1, max([2, ceil(geometry.user_data.l / dx) + 1]));
%!   geometry.mesh_size.s = linspace(0, 1, max([2, ceil(geometry.user_data.w / dx) + 1]));
%!   geometry.mesh_size.t = linspace(0, 1, max([2, ceil(geometry.user_data.h / dx) + 1]));
%!   geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.l;
%!   geometry.spatial_coordinates = @(r, s, t, varargin) feval("cube_geo", geometry.user_data, r, s, t);
%!   geometry.material_selector = @(r, s, t, varargin) 1;
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("boundary_cond_callback", r, s, t, geometry, load_data);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("pressure_callback", r, s, t, geometry, load_data, perm_idx);
%!   K0 = 5000;
%!   material.E = 210000e6;
%!   material.nu = 0.3;
%!   material.rho = 7850;
%!   material.k = diag([K0, K0, K0]);
%!   material.cp = 465;
%!   u0 = 100;
%!   ub = 50;
%!   load_data.pressure = 0;
%!   [mesh] = fem_pre_mesh_struct_create(geometry, load_data, material, options);
%!   [mesh] = fem_pre_mesh_solid_to_surf(mesh);
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   [dof_map] = fem_ass_dof_map(mesh, load_case);
%!   empty_cell = cell(1, 2);
%!   sol = struct("t", empty_cell, "theta", empty_cell);
%!   for j=1:2
%!     switch (j)
%!       case 1
%!         R = eye(3);
%!       case 2
%!         e1 = [0.5; 0.2; 0.1];
%!         e2 = [0; 1; 0];
%!         e3 = cross(e1, e2);
%!         e2 = cross(e3, e1);
%!         R = [e1 / norm(e1), e2 / norm(e2), e3 / norm(e3)];
%!     endswitch
%!     cond_b = (mesh.nodes(:, 1) == 0 | mesh.nodes(:, 1) == geometry.user_data.l);
%!     idx_b = find(cond_b);
%!     idx_i = find(~cond_b);
%!     x = mesh.nodes(:, 1);
%!     mesh.nodes = [mesh.nodes(:, 1:3) * R.', mesh.nodes(:, 4:6) * R.'];
%!     [mat_ass.Kk, mat_ass.C] = fem_ass_matrix(mesh, ...
%!                                              dof_map, ...
%!                                              [FEM_MAT_THERMAL_COND, ...
%!                                               FEM_MAT_HEAT_CAPACITY], ...
%!                                              load_case);
%!     Kk11 = mat_ass.Kk(idx_i, idx_i);
%!     Kk12 = mat_ass.Kk(idx_i, idx_b);
%!     C11 = mat_ass.C(idx_i, idx_i);
%!     theta_b = repmat(ub, numel(idx_b), 1);
%!     theta0 = repmat(u0, dof_map.totdof, 1);
%!     qref = material.rho * material.cp * geometry.user_data.l * geometry.user_data.w * geometry.user_data.h;
%!     assert_simple(sum(sum(mat_ass.C)), qref, eps^0.5 * abs(qref));
%!     kappa = K0 / (material.rho * material.cp);
%!     T_ = geometry.user_data.l^2 / kappa;
%!     dt = dx^2 / (2 * kappa);
%!     alpha = 0.6;
%!     sol(j).t = 0:dt:0.1*T_;
%!     sol(j).theta = zeros(dof_map.totdof, numel(sol(j).t));
%!     sol(j).theta(:, 1) = theta0;
%!     A = (1 / dt) * C11 + alpha * Kk11;
%!     opts.number_of_threads = mbdyn_solver_num_threads_default();
%!     opts.solver = "pardiso";
%!     Afact = fem_sol_factor(A, opts);
%!     for i=2:numel(sol(j).t)
%!       sol(j).theta(idx_i, i) = Afact \ (C11 * (sol(j).theta(idx_i, i - 1) / dt) - Kk11 * (sol(j).theta(idx_i, i - 1) * (1 - alpha)) - Kk12 * theta_b);
%!       sol(j).theta(idx_b, i) = theta_b;
%!     endfor
%!     u0_ = 1;
%!     n = 1:10000;
%!     Bn_ = -2 * u0_ ./ (n * pi) .* ((-1).^n - 1);
%!     x_ = x / geometry.user_data.l;
%!     u_ = zeros(numel(x), numel(sol(j).t));
%!     for n=1:numel(Bn_)
%!       t_ = sol(j).t / T_;
%!       u_ += Bn_(n) * sin(n * pi * x_) .* exp(-n^2 * pi^2 * t_);
%!     endfor
%!     sol(j).theta_ref = u_ * (u0 - ub) + ub;
%!   endfor
%!   [x, idx_theta] = sort(x);
%!   if (do_plot)
%!   for j=1:numel(sol)
%!     colors = rainbow(numel(sol(j).t));
%!     figure("visible", "off");
%!     hold on;
%!     legend off;
%!     for i=1:10:numel(sol(j).t)
%!       set(plot(x, sol(j).theta(idx_theta, i), sprintf("-;t=%.2f;", sol(j).t(i))), "color", colors(i,:));
%!       hnd = plot(x, sol(j).theta_ref(idx_theta, i), sprintf("--;tref=%.2f;", sol(j).t(i)));
%!       set(hnd, "color", colors(i,:));
%!       set(hnd, "linewidth", 2);
%!     endfor
%!     xlabel("t [s]");
%!     ylabel("theta [degC]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("transient thermal problem of a bar: %s/%d", elem_types{k}, j));
%!   endfor
%!   endif
%!   tol = 1e-2;
%!   assert_simple(sol(j).theta(:, 10:end), sol(j).theta_ref(:, 10:end), tol * abs(u0 - ub));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
