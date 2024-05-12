## fem_tests.m:30
%!test
%! try
%! close all;
%! ## TEST 30
%! ###############################################
%! ## Stress and strain of 10 node tetrahedrons
%! ###############################################
%! a = 10e-3;
%! b = 20e-3;
%! c = 7e-3;
%! u = [2e-3;  5e-3;  3e-3];
%! l = [a; b; c];
%! x0 = [5e-3, -8e-3, -17e-3];
%! v0 = [13e-3, -16e-3, 20e-3];
%! xi =         [      0,       0,       0;
%!                     a,       0,       0;
%!                     0,       b,       0;
%!                     0,       0,       c;
%!                     0.5 * a,       0,       0;
%!                     0.5 * a, 0.5 * b,       0;
%!                     0, 0.5 * b,       0;
%!                     0,       0, 0.5 * c;
%!                     0.5 * a,       0, 0.5 * c;
%!                     0, 0.5 * b, 0.5 * c];
%! rand("seed", 0);
%! Phi = rand(3, 5) * pi / 180;
%! for r=1:columns(Phi)
%!   R = euler123_to_rotation_matrix(Phi(:, r));
%!   assert_simple(R*R.', eye(3), eps^0.8);
%!   assert_simple(R.'*R, eye(3), eps^0.8);
%!   data(r).mesh.nodes = [xi * R.' + x0, zeros(rows(xi), 3)];
%!   data(r).mesh.elements.tet10 = int32(1:10);
%!   data(r).mesh.materials.tet10 = int32(1);
%!   E = 210000e6;
%!   nu = 0.3;
%!   data(r).mesh.material_data.rho = 1;
%!   data(r).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data(r).load_case.locked_dof = false(rows(data(r).mesh.nodes), 6);
%!   [data(r).dof_map] = fem_ass_dof_map(data(r).mesh, data(r).load_case);
%!   data(r).sol_stat.def = zeros(rows(xi), 6, numel(l) * numel(u));
%!   for j=1:numel(l)
%!     for i=1:numel(u)
%!       v = u(i) * xi(:, j) / l(j);
%!       V = zeros(rows(v), 3);
%!       V(:, i) = v;
%!       V += v0;
%!       V *= R.';
%!       data(r).sol_stat.def(:, 1:3, (j - 1) * numel(u) + i) = V;
%!     endfor
%!   endfor
%!   [data(r).sol_stat.stress] = fem_ass_matrix(data(r).mesh, ...
%!                                              data(r).dof_map, ...
%!                                              [FEM_VEC_STRESS_CAUCH], ...
%!                                              data(r).load_case, ...
%!                                              data(r).sol_stat);
%!   idx = [1, 4, 6;
%!          4, 2, 5;
%!          6, 5, 3];
%!   for i=1:numel(u)
%!     for j=1:numel(l)
%!       Epsilon_a_0 = zeros(3, 3);
%!       Epsilon_a_0(i, j) = u(i) / l(j);
%!       if (i ~= j)
%!         Epsilon_a_0(i, j) /= 2;
%!         Epsilon_a_0(j, i) = Epsilon_a_0(i, j);
%!       endif
%!       Epsilon_a = R * Epsilon_a_0 * R.';
%!       epsilon_a = [Epsilon_a(1, 1);
%!                    Epsilon_a(2, 2);
%!                    Epsilon_a(3, 3);
%!                    2 * Epsilon_a(1, 2);
%!                    2 * Epsilon_a(2, 3);
%!                    2 * Epsilon_a(3, 1)];
%!       epsilon_a_0 = [Epsilon_a_0(1, 1);
%!                      Epsilon_a_0(2, 2);
%!                      Epsilon_a_0(3, 3);
%!                      2 * Epsilon_a_0(1, 2);
%!                      2 * Epsilon_a_0(2, 3);
%!                      2 * Epsilon_a_0(3, 1)];
%!       tau_a_0 = data(r).mesh.material_data.C * epsilon_a_0;
%!       Tau_a = R * tau_a_0(idx) * R.';
%!       tau_a = [Tau_a(1, 1);
%!                Tau_a(2, 2);
%!                Tau_a(3, 3);
%!                Tau_a(1, 2);
%!                Tau_a(2, 3);
%!                Tau_a(3, 1)];
%!       for n=1:size(data(r).sol_stat.stress.tau.tet10, 2)
%!         for k=1:size(data(r).sol_stat.stress.tau.tet10, 1)
%!           tau = data(r).sol_stat.stress.tau.tet10(k, n, :, (j - 1) * numel(u) + i)(:);
%!           Tau_0 = R.' * tau(idx) * R;
%!           tau_0 = [Tau_0(1, 1);
%!                    Tau_0(2, 2);
%!                    Tau_0(3, 3);
%!                    Tau_0(1, 2);
%!                    Tau_0(2, 3);
%!                    Tau_0(3, 1)];
%!           epsilon_0 = data(r).mesh.material_data.C \ tau_0;
%!           epsilon = data(r).mesh.material_data.C \ tau;
%!           tau_eps_a = data(r).mesh.material_data.C * epsilon_a;
%!           assert_simple(tau, tau_a,  max(abs(tau_a)) * (eps)^0.75);
%!           assert_simple(tau, tau_eps_a,  max(abs(tau_eps_a)) * (eps)^0.75);
%!           assert_simple(tau_0, tau_a_0,  max(abs(tau_a_0)) * (eps)^0.75);
%!           assert_simple(epsilon_0, epsilon_a_0,  max(abs(epsilon_a_0)) * (eps)^0.75);
%!           assert_simple(epsilon, epsilon_a,  max(abs(epsilon_a)) * (eps)^0.75);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
