## fem_tests.m:42
%!test
%! ## TEST 42
%! close all;
%! scale_stat = 1;
%! scale_eig = 250e-3;
%! E = 210000e6;
%! nu = 0.3;
%! material.C = fem_pre_mat_isotropic(E, nu);
%! material.rho = 7850;
%! Fy = 15000;
%! h = 10e-3 / 2;
%! geometry.l = 1000e-3;
%! geometry.w = 10e-3;
%! geometry.h = 50e-3;
%! A = geometry.w * geometry.h;
%! Wz = geometry.h * geometry.w^2 / 6;
%! Iz = geometry.h * geometry.w^3 / 12;
%! tauxx_max = -Fy * geometry.l / Wz;
%! tauxy_mean = Fy / A;
%! Uy = Fy * geometry.l^3 / (3 * E * Iz);
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! number_of_modes = 10;
%! number_of_modes_disp = 3;
%! plot_def = false;
%! f = [ 0; Fy; 0 ];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! alg = {"shift-invert", "symmetric-inverse", "unsymmetric"};
%! rho = 100;
%! tol = 1e-6;
%! err = zeros(number_of_modes, numel(alg));
%! for a=1:numel(alg)
%!   [sol_eig(a), ~, err(:, a)] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes, rho, tol, alg{a});
%! endfor
%! z = linspace(0,geometry.l,100);
%! I = [ geometry.w * geometry.h, geometry.h * geometry.w^3 / 12, geometry.w * geometry.h^3 / 12 ];
%! y(1,:) = f(1) * geometry.l / ( E * I(1) ) * ( 1 - z / geometry.l );
%! for i=2:3
%!   y(i,:) = f(i) * geometry.l^3 / ( 6 * E * I(i) ) * ( 2 - 3 * z / geometry.l + ( z / geometry.l ).^3 );
%! endfor
%! B = E * I(2:3);
%! my = material.rho * I(1);
%! P0 = ( 0.3 * geometry.l ) / geometry.l^3 * 3 * B;
%! omega1 = sqrt(B / (my * geometry.l^4));
%! omega_ref = omega1.' * [3.516, 22.035, 61.697];
%! omega_ref = sort(reshape(omega_ref, 1, numel(omega_ref)));
%! f_ref = omega_ref / (2 * pi);
%! for a=1:numel(sol_eig)
%!   assert_simple(sol_eig(a).f(1:5), f_ref(1:5), 4e-2 * max(f_ref(1:5)));
%! endfor
%! if (plot_def)
%!   figure("visible","off");
%!   hold on;
%!   fem_post_sol_plot(mesh);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title('undeformed mesh');
%!   figure("visible","off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat, scale_stat);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title('deformed mesh');
%!   for i=1:min(number_of_modes_disp,length(sol_eig(1).f))
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig(1), scale_eig / max(norm(sol_eig(1).def(:, :, i), "rows")),i);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("%d. eigenmode: %gHz",i,sol_eig(1).f(i)));
%!   endfor
%!   figure_list();
%! endif
