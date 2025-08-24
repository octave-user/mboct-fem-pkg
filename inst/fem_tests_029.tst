## fem_tests.m:29
%!test
%! try
%! pkg load mbdyn_util_oct;
%! ## TEST 29
%! close all;
%! a = 40e-3;
%! b = 10e-3;
%! c = 30e-3;
%! ux = -3e-3;
%! uy = 5e-3;
%! uz = 7e-3;
%! epsilonxx = ux / a;
%! epsilonyy = 2 * uy / b;
%! epsilonzz = -2 * uz / c;
%! epsilonxy = uy / a;
%! epsilonyz = 0;
%! epsilonzx = uz / a;
%! epsilon_a = [epsilonxx,  epsilonxx,  epsilonxx,  epsilonxx,  epsilonxx,  epsilonxx,  epsilonxx,  epsilonxx,
%!              epsilonyy,          0,          0,  epsilonyy,  epsilonyy,          0,          0,  epsilonyy,
%!              epsilonzz,          0,          0,  epsilonzz,  epsilonzz,          0,          0,  epsilonzz,
%!              epsilonxy,  epsilonxy, -epsilonxy, -epsilonxy,  epsilonxy,  epsilonxy, -epsilonxy, -epsilonxy,
%!              epsilonyz,  epsilonyz,  epsilonyz,  epsilonyz,  epsilonyz,  epsilonyz,  epsilonyz,  epsilonyz,
%!              -epsilonzx, -epsilonzx, -epsilonzx, -epsilonzx,  epsilonzx,  epsilonzx,  epsilonzx,  epsilonzx];
%! rho = 7850;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  #  1
%!       -0.5 * a,  0.5 * b,  0.5 * c;  #  2
%!       -0.5 * a, -0.5 * b,  0.5 * c;  #  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  #  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  #  5
%!       -0.5 * a,  0.5 * b, -0.5 * c;  #  6
%!       -0.5 * a, -0.5 * b, -0.5 * c;  #  7
%!       0.5 * a, -0.5 * b, -0.5 * c]; #  8
%! R = euler123_to_rotation_matrix([25; 30; 45]*pi/180);
%! mesh.nodes = [X * R.', zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8]);
%! mesh.materials.iso8 = int32([1]);
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = rho;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS], ...
%!                                      load_case);
%! U = zeros(rows(mesh.nodes), 3);
%! U(1, 1:3) = [ux, uy, -uz];
%! U(4, 1:3) = [ux, -uy, -uz];
%! U(5, 1:3) = [ux, uy, uz];
%! U(8, 1:3) = [ux, -uy, uz];
%! sol_stat.def = [U * R.', zeros(rows(mesh.nodes), 3)];
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! tau_a = mesh.material_data.C * epsilon_a;
%! idx = [1, 4, 6;
%!        4, 2, 5;
%!        6, 5, 3];
%! for i=1:columns(tau_a)
%!   tau_a_tens = R * tau_a(:, i)(idx) * R.';
%!   tau_a(:, i) = [tau_a_tens(1, 1);
%!                  tau_a_tens(2, 2);
%!                  tau_a_tens(3, 3);
%!                  tau_a_tens(1, 2);
%!                  tau_a_tens(2, 3);
%!                  tau_a_tens(1, 3)];
%! endfor
%! for i=1:columns(tau_a)
%!   assert_simple(tau_a(:, i), sol_stat.stress.tau.iso8(1, i, :)(:), sqrt(eps) * norm(tau_a(:, i)));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
