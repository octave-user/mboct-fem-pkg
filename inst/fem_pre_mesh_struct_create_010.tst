## fem_pre_mesh_struct_create.m:10
%!test
%! try
%! ## TEST10
%! ## distorted sphere
%! close all;
%! material.E = 28e6 * 6895;
%! material.nu = 0.28;
%! material.rho = 0.000751 * 4.4482 / (25.4e-3^4);
%! geo.a = 2 * 4.4688 * 25.4e-3;
%! geo.b = 1 * 4.4688 * 25.4e-3;
%! geo.c = 1.5 * 4.4688 * 25.4e-3;
%! geo.t = 0.625 * 25.4e-3;
%! geo.Psii = [-180, -135, -90, -45, 0,  45, 90, 135, 180] * pi / 180;
%! geo.Zetai = [-90, -45, 0, 45, 90] * pi / 180;
%! geo.Psir = [1, 1, 1, 1.05, 1, 1.05,  1, 1;
%!             1, 1, 1, 1.10, 1, 1.10,  1, 1;
%!             1, 1, 1, 1.05, 1.2, 1.05,  1, 1];
%! geo.Psir = [geo.Psir, geo.Psir(:, 1)];
%! geo.Psir = [ones(1, columns(geo.Psir));
%!             geo.Psir;
%!             ones(1, columns(geo.Psir))];
%! geo.Psii = [geo.Psii(1:end - 1) - 2 * pi, geo.Psii, geo.Psii(2:end) + 2 * pi];
%! geo.Psir = [geo.Psir(:, 1:end - 1), geo.Psir, geo.Psir(:, 2:end)];
%! N = 10;
%! f_run_post_proc = false;
%! function X = sphere_shape(geo, Psi, Zeta, varargin)
%!   Psir = sqrt(1 - Zeta^2) * interp2(geo.Psii, geo.Zetai, geo.Psir, Psi, Zeta,  "pchip");
%!   if (~(isreal(Psir) && isfinite(Psir)))
%!    error("interpolation failed");
%!   endif
%!   X = zeros(3, 1);
%!   X(1) = Psir * geo.a * cos(Psi);
%!   X(2) = Psir * geo.b * sin(Psi);
%!   X(3) = geo.c * Zeta;
%! endfunction
%!
%! function [x, y, z] = sphere_geo(geo, r, s, t, varargin)
%!   Zeta = sin(t);
%!   Psi = s;
%!   dZeta = dPsi = 2 * pi * sqrt(eps);
%!   if (Zeta + dZeta > 1)
%!     dZeta = -dZeta;
%!   endif
%!   if (Psi + dPsi > pi)
%!     dPsi = -dPsi;
%!   endif
%!   P0 = sphere_shape(geo, Psi, Zeta);
%!   P1 = sphere_shape(geo, Psi, Zeta + dZeta);
%!   P2 = sphere_shape(geo, Psi + dPsi, Zeta);
%!   n1 = (P1 - P0) * sign(dZeta);
%!   n2 = (P2 - P0) * sign(dPsi);
%!   n3 = cross(n1, n2);
%!   if (norm(n3) == 0)
%!     if (t == -pi/2)
%!       n3 = [0; 0; 1];
%!     elseif (t == pi/2)
%!       n3 = [0; 0; -1];
%!     else
%!       error("surface normal vector is singular");
%!     endif
%!   endif
%!   n3 /= norm(n3);
%!   P = P0 - n3 * geo.t * r;
%!   if (~isfinite(P))
%!     error("P is not finite");
%!   endif
%!   x = P(1);
%!   y = P(2);
%!   z = P(3);
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
%! figure("visible", "off");
%! hold on;
%! Psi = linspace(-pi, pi, 180);
%! Zeta = linspace(-1, 1, 9);
%! for j=1:numel(Zeta)
%!   X = zeros(3, numel(Psi));
%!   for i=1:numel(Psi)
%!     X(:, i) = sphere_shape(geo, Psi(i), Zeta(j));
%!   endfor
%!   set(plot3(X(1, :), X(2, :), X(3, :)), "color", [1, 0, 0]);
%! endfor
%! Zeta = linspace(-1, 1, 50);
%! Psi = linspace(-pi, pi, 18);
%! for j=1:numel(Psi)
%!   X = zeros(3, numel(Zeta));
%!   for i=1:numel(Zeta)
%!     X(:, i) = sphere_shape(geo, Psi(j), Zeta(i));
%!   endfor
%!   set(plot3(X(1, :), X(2, :), X(3, :)), "color", [0, 0, 1]);
%! endfor
%! daspect(ones(1, 3));
%! title("wireframe");
%! xlabel("x");
%! ylabel("y");
%! zlabel("z");
%!   geometry.mesh_size.r = linspace(0, 1, 2);
%!   geometry.mesh_size.s = linspace(-pi, pi, 18);
%!   geometry.mesh_size.t = linspace(-pi/2, pi/2, 18);
%!   load_data.pressure = 1;
%!   geometry.sewing.tolerance = sqrt(eps) * geo.a;
%!   geometry.spatial_coordinates = @(r, s, t, varargin) feval("sphere_geo", geo, r, s, t, varargin);
%!   geometry.material_selector = @(r, s, t, varargin) int32(1);
%!   geometry.boundary_condition =  @(r, s, t, geometry, load_data, varargin) feval("sphere_bound", r, s, t, geo, load_data, varargin);
%!   geometry.pressure_boundary_condition = @(r, s, t, geometry, load_data, perm_idx, varargin) feval("sphere_pressure", r, s, t, geo, load_data, perm_idx, varargin);
%!   options.elem_type = "iso20";
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
%!   shift = sqrt(eps) * max(abs(diag(mat_ass.K))) / max(abs(diag(mat_ass.M)));
%!   tol = eps^0.4;
%!   alg = "shift-invert";
%!   solver = "pardiso";
%!   num_threads = mbdyn_solver_num_threads_default();
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, shift, tol, alg, solver, num_threads);
%!   sol_eig.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_eig);
%!   opts.scale_def = 0.5 * geo.a / max(max(max(abs(sol_eig.def))));
%!   opts.print_and_exit = true;
%!   opts.print_to_file = "";
%!   opts.skin_only = true;
%!   opts.show_element = true;
%! if (f_run_post_proc)
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
%! figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
