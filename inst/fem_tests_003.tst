## fem_tests.m:03
%!test
%! try
%! ## TEST 3
%! ##########################################################################################
%! ## Test case for elimination of Lagrange multipliers
%! ##########################################################################################
%! close all;
%! plot_def = false;
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 140e-3;
%! scale_def = 100e-3;
%! tol = eps^0.7;
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
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8;
%!                             9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3.nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3.weight = ones(1, 4);
%! n1 = [0.5; 0.5; 0];
%! n1 /= norm(n1);
%! n2 = [1; 1; 1];
%! n3 = cross(n1, n2);
%! n2 = cross(n1, n3);
%! n2 /= norm(n2);
%! n3 /= norm(n3);
%! mesh.elements.joints(1).nodes = int32([13]);
%! mesh.elements.joints(1).C = [[1, 0, 0;
%!                               0, 0, 1],    zeros(2, 3);
%!                              zeros(2, 3), [n2.'; n3.']];
%! load_case.joints(1).U = [0; 0; 0; 0];
%! for i=[2, 3, 6, 7]
%!   mesh.elements.joints(end + 1).nodes = int32(i);
%!   mesh.elements.joints(end).C = [eye(3), zeros(3, 3)];
%!   load_case.joints(end + 1).U = [0; 0; 0];
%! endfor
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case.loads = repmat([0, 0, 10,  0,   0, 0], 4, 1);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.R, ...
%!  mat_ass.dm, ...
%!  mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_STIFFNESS, ...
%!                                      FEM_MAT_MASS, ...
%!                                      FEM_VEC_LOAD_CONSISTENT, ...
%!                                      FEM_SCA_TOT_MASS], ...
%!                                     load_case);
%! assert_simple(rank(mat_ass.K), columns(mat_ass.K));
%! [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%! assert_simple(isdefinite(Kred), true);
%! assert_simple(isdefinite(Mred), true);
%! lambda = eig(mat_ass.K, mat_ass.M);
%! lambda = sort(lambda(find(isfinite(lambda))));
%! lambdared = eig(Kred, Mred);
%! lambdared = sort(lambdared);
%! U = full(mat_ass.K \ mat_ass.R);
%! Ured = Tred * (Kred \ Rred);
%! assert_simple(Ured, U(dof_map.idx_node), tol * norm(U(dof_map.idx_node)));
%! assert_simple(lambda(1:length(lambdared)), lambdared, tol * max(abs(lambda(1:length(lambdared)))));
%! sol_stat.def = fem_post_def_nodal(mesh, dof_map, full(mat_ass.K \ mat_ass.R));
%! assert_simple(mat_ass.dm, mesh.material_data.rho * a * b * c, sqrt(eps) * mesh.material_data.rho * a * b * c);
%! if (plot_def)
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "cols")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   title("deformed mesh");
%!   figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
