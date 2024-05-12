## fem_tests.m:24
%!test
%! try
%! ## TEST 24
%! close all;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! rho = 7850;
%! Fy = 1;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  #  1
%!       0,  0.5 * b,  0.5 * c;  #  2
%!       0, -0.5 * b,  0.5 * c;  #  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  #  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  #  5
%!       0,  0.5 * b, -0.5 * c;  #  6
%!       0, -0.5 * b, -0.5 * c;  #  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  #  8
%!       a,  0.5 * b,  0.5 * c;  #  9
%!       a, -0.5 * b,  0.5 * c;  # 10
%!       a,  0.5 * b, -0.5 * c;  # 11
%!       a, -0.5 * b, -0.5 * c]; # 12
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8;
%!                             9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = rho;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(find(X(:, 1) == 0), 1) = true;
%! load_case.locked_dof(find(X(:, 2) == -0.5 * b), 2) = true;
%! load_case.locked_dof(find(X(:, 3) == -0.5 * c), 3) = true;
%! load_case.pressure.iso4.elements = int32([6, 5, 1, 2; ...
%!                                           5, 11, 9, 1]);
%! load_case.pressure.iso4.p = repmat(Fy / (a * c), 2, 4);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! tauyy_a = Fy / (a * c);
%! assert_simple(all(abs(sol_stat.stress.tau.iso8(:, :, 2) / tauyy_a - 1) < sqrt(eps)));
%! assert_simple(all(abs(sol_stat.stress.tau.iso8(:, :, [1,3:end]) / tauyy_a) < sqrt(eps)));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
