## fem_tests.m:23
%!test
%! ## TEST 23
%! close all;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! rho = 7850;
%! Fx = 1;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;
%!       0,  0.5 * b,  0.5 * c;
%!       0, -0.5 * b,  0.5 * c;
%!       0.5 * a, -0.5 * b,  0.5 * c;
%!       0.5 * a,  0.5 * b, -0.5 * c;
%!       0,  0.5 * b, -0.5 * c;
%!       0, -0.5 * b, -0.5 * c;
%!       0.5 * a, -0.5 * b, -0.5 * c,
%!       a,  0.5 * b,  0.5 * c;
%!       a, -0.5 * b,  0.5 * c;
%!       a,  0.5 * b, -0.5 * c;
%!       a, -0.5 * b, -0.5 * c];
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
%! load_case.loaded_nodes = int32([9:12]).';
%! load_case.loads = [repmat(Fx / 4, 4, 1), zeros(4, 5)];
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
%! tauxx_a = Fx / (b * c);
%! assert_simple(all(abs(sol_stat.stress.tau.iso8(:, :, 1) / tauxx_a - 1) < sqrt(eps)));
%! assert_simple(all(abs(sol_stat.stress.tau.iso8(:, :, 2:end) / tauxx_a) < sqrt(eps)));
