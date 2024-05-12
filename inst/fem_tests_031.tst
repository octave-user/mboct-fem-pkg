## fem_tests.m:31
%!test
%! try
%! ## TEST 31
%! close all;
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! R = euler123_to_rotation_matrix([30; 40; 70]*pi/180);
%! a = 15e-3 / SI_unit_m;
%! b = 25e-3 / SI_unit_m;
%! c = 18e-3 / SI_unit_m;
%! ux = 12e-3 / SI_unit_m;
%! X = [ a,  0.7 * b,  0.5 * c;  ## 1
%!       0,  0.7 * b,  0.5 * c;  ## 2
%!       0, -0.7 * b,  0.5 * c;  ## 3
%!       a, -0.7 * b,  0.5 * c;  ## 4
%!       a,  0.3 * b, -0.5 * c;  ## 5
%!       0,  0.3 * b, -0.5 * c;  ## 6
%!       0, -0.3 * b, -0.5 * c;  ## 7
%!       a, -0.3 * b, -0.5 * c]; ## 8
%! mesh.nodes = [X * R.', zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8]);
%! mesh.materials.iso8 = int32([1]);
%! E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! sol_stat.def = zeros(rows(mesh.nodes), 6);
%! sol_stat.def([1,4,6,7], 1) = 0.5 * ux;
%! sol_stat.def([2,3,5,8], 1) = -0.5 * ux;
%! for i=1:rows(sol_stat.def)
%!   sol_stat.def(i, 1:3) *= R.';
%! endfor
%! epsilonxx_a = ux / a;
%! epsilonxz_a = ux / c;
%! epsilon_a = [epsilonxx_a, epsilonxx_a, -epsilonxx_a, -epsilonxx_a;
%!              zeros(4, 4);
%!              epsilonxz_a, -epsilonxz_a, epsilonxz_a, -epsilonxz_a];
%! tau_a = mesh.material_data.C * epsilon_a;
%! idx = [1, 4, 6;
%!        4, 2, 5;
%!        6, 5, 3];
%! tau_a_0 = zeros(size(tau_a));
%! for i=1:columns(tau_a)
%!   Tau_a_0 = R * tau_a(:, i)(idx) * R.';
%!   tau_a_0(:, i) = [Tau_a_0(1, 1);
%!                    Tau_a_0(2, 2);
%!                    Tau_a_0(3, 3);
%!                    Tau_a_0(1, 2);
%!                    Tau_a_0(2, 3);
%!                    Tau_a_0(3, 1)];
%! endfor
%! [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%! tau = sol_stat.stress.tau.iso8;
%! tol = sqrt(eps) * norm(tau_a);
%! assert_simple(tau(1, 1, :)(:), tau_a_0(:, 1), tol);
%! assert_simple(tau(1, 4, :)(:), tau_a_0(:, 1), tol);
%! assert_simple(tau(1, 2, :)(:), tau_a_0(:, 2), tol);
%! assert_simple(tau(1, 3, :)(:), tau_a_0(:, 2), tol);
%! assert_simple(tau(1, 5, :)(:), tau_a_0(:, 3), tol);
%! assert_simple(tau(1, 8, :)(:), tau_a_0(:, 3), tol);
%! assert_simple(tau(1, 6, :)(:), tau_a_0(:, 4), tol);
%! assert_simple(tau(1, 7, :)(:), tau_a_0(:, 4), tol);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
