## fem_tests.m:04
%!test
%! try
%! ##########################################################################################
%! ## TEST 4: Test case for elimination of Lagrange multipliers
%! ##########################################################################################
%! close all;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 140e-3;
%! scale_def = 20e3;
%! tol = eps^0.7;
%! tol2 = eps^0.6;
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
%!       a, -0.5 * b, -0.5 * c,  # 12
%!       d,        0,        0]; # 13
%! Phi1 = [0, 20, 45] * pi / 180;
%! Phi2 = [0, 60, 270] * pi / 180;
%! Phi3 = [0, -15, 30] * pi / 180;
%! for j=1:length(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   assert_simple(R1.' * R1, eye(3), tol);
%!   assert_simple(R1 * R1.', eye(3), tol);
%!   data(j).T1 = [R1, zeros(3, 3);
%!                 zeros(3, 3), R1];
%!   data(j).mesh.nodes = [X * R1.', zeros(rows(X), 3)];
%!   data(j).mesh.elements.iso8 = int32([1:8;
%!                                       9, 1, 4, 10, 11, 5, 8, 12]);
%!   data(j).mesh.materials.iso8 = int32([1; 1]);
%!   data(j).mesh.elements.rbe3.nodes = int32([13, 9, 10, 11, 12]);
%!   data(j).mesh.elements.rbe3.weight = ones(1, 4);
%!   n1 = [0; 0.5; 0.2];
%!   n1 /= norm(n1);
%!   n2 = [1; 1; 1];
%!   n3 = cross(n1, n2);
%!   n2 = cross(n1, n3);
%!   n2 /= norm(n2);
%!   n3 /= norm(n3);
%!   t1 = [1; 0; 0];
%!   t2 = [0; 0; 1];
%!   data(j).mesh.elements.joints(1).nodes = int32([13]);
%!   data(j).mesh.elements.joints(1).C = [[(R1 * t1).'; (R1 * t2).'],    zeros(2, 3);
%!                                        zeros(2, 3), [(R1 * n2).'; (R1 * n3).']];
%!   for i=[2, 3, 6, 7]
%!     data(j).mesh.elements.joints(end + 1).nodes = int32(i);
%!     data(j).mesh.elements.joints(end).C = [R1.', zeros(3, 3)];
%!   endfor
%!   E = 210000e6;
%!   nu = 0.3;
%!   data(j).mesh.material_data.rho = 7850;
%!   data(j).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data(j).load_case.locked_dof = false(rows(data(j).mesh.nodes), 6);
%!   data(j).load_case.loaded_nodes = int32([1; 4; 5; 8]);
%!   data(j).load_case.loads = repmat([0, 0, 10,  0,   0, 0] * data(j).T1.', 4, 1);
%!   data(j).dof_map = fem_ass_dof_map(data(j).mesh, data(j).load_case);
%!   [data(j).mat_ass.K, ...
%!    data(j).mat_ass.M, ...
%!    data(j).mat_ass.R, ...
%!    data(j).mat_ass.dm, ...
%!    data(j).mat_ass.mat_info] = fem_ass_matrix(data(j).mesh, ...
%!                                               data(j).dof_map, ...
%!                                               [FEM_MAT_STIFFNESS, ...
%!                                                FEM_MAT_MASS, ...
%!                                                FEM_VEC_LOAD_CONSISTENT, ...
%!                                                FEM_SCA_TOT_MASS], ...
%!                                               data(j).load_case);
%!   assert_simple(data(j).mat_ass.dm, data(j).mesh.material_data.rho * a * b * c, sqrt(eps) * data(j).mesh.material_data.rho * a * b * c);
%!   assert_simple(rank(data(j).mat_ass.K), columns(data(j).mat_ass.K));
%!   [data(j).Tred, data(j).Kred, data(j).Mred, data(j).Rred] = fem_cms_constr_elim(data(j).mesh, data(j).dof_map, data(j).mat_ass);
%!   assert_simple(isdefinite(data(j).Kred), true);
%!   assert_simple(isdefinite(data(j).Mred), true);
%!   [data(j).Phi, data(j).lambda] = eig(data(j).mat_ass.K, data(j).mat_ass.M);
%!   data(j).lambda = diag(data(j).lambda);
%!   idx_lambda = find(isfinite(data(j).lambda));
%!   data(j).Phi = data(j).Phi(:, idx_lambda);
%!   [data(j).lambda, idx_lambda] = sort(data(j).lambda(idx_lambda));
%!   data(j).Phi = data(j).Phi(:, idx_lambda);
%!   for k=1:columns(data(j).Phi)
%!     data(j).Phi(:, k) /= norm(data(j).Phi(data(j).dof_map.idx_node, k));
%!   endfor
%!   [data(j).Phired, data(j).lambdared] = eig(data(j).Kred, data(j).Mred);
%!   data(j).lambdared = diag(data(j).lambdared);
%!   data(j).lambdared = sort(data(j).lambdared);
%!   data(j).U = full(data(j).mat_ass.K \ data(j).mat_ass.R);
%!   data(j).Ured = data(j).Tred * (data(j).Kred \ data(j).Rred);
%!   assert_simple(data(j).Ured, data(j).U(data(j).dof_map.idx_node), tol * norm(data(j).U(data(j).dof_map.idx_node)));
%!   assert_simple(data(j).lambda(1:length(data(j).lambdared)), data(j).lambdared, tol * max(abs(data(j).lambda(1:length(data(j).lambdared)))));
%!   data(j).sol_stat.def = fem_post_def_nodal(data(j).mesh, data(j).dof_map, full(data(j).mat_ass.K \ data(j).mat_ass.R)) * data(j).T1;
%!   for k=1:columns(data(j).Phi)
%!     data(j).modal(k).f = sqrt(data(j).lambda(k)) / (2 * pi);
%!     data(j).modal(k).def = fem_post_def_nodal(data(j).mesh, data(j).dof_map, data(j).Phi(:, k)) * data(j).T1;
%!   endfor
%! endfor
%! for j=2:length(Phi1)
%!   assert_simple(data(j).sol_stat.def, data(1).sol_stat.def, tol * max(max(max(abs(data(1).sol_stat.def)))));
%!   N = min([length(data(1).lambda), length(data(j).lambda)]);
%!   assert_simple(data(j).lambda(1:N), data(1).lambda(1:N), tol * max(abs(data(1).lambda(1:N))));
%!   assert_simple(data(j).lambdared, data(1).lambdared, tol * max(abs(data(1).lambdared)));
%!   for k=1:N
%!     assert_simple(max(max(max(abs(data(j).modal(k).def - data(1).modal(k).def)))) < tol2 * max(max(max(abs(data(1).modal(k).def)))) || max(max(max(abs(data(j).modal(k).def + data(1).modal(k).def)))) < tol2 * max(max(max(abs(data(1).modal(k).def)))));
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
