## fem_tests.m:27
%!test
%! ## TEST 27
%! close all;
%! a = 20e-3;
%! b = 20e-3;
%! c = 20e-3;
%! uy = 7e-3;
%! rho = 7850;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  #  1
%!       -0.5 * a,  0.5 * b,  0.5 * c;  #  2
%!       -0.5 * a, -0.5 * b,  0.5 * c;  #  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  #  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  #  5
%!       -0.5 * a,  0.5 * b, -0.5 * c;  #  6
%!       -0.5 * a, -0.5 * b, -0.5 * c;  #  7
%!       0.5 * a, -0.5 * b, -0.5 * c]; #  8
%! mesh.nodes = [X, zeros(rows(X), 3)];
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
%! sol_stat.def = zeros(size(mesh.nodes));
%! sol_stat.def(1:4, 2) = uy;
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! G = E / (2 * (1 + nu));
%! gamma = uy / c;
%! tauyz_a = G * gamma;
%! assert_simple(all(abs(sol_stat.stress.tau.iso8(:,:,6) / tauyz_a - 1) < sqrt(eps) * abs(tauyz_a)));
