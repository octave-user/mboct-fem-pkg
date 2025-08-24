## fem_tests.m:11
%!test
%! try
%! pkg load mbdyn_util_oct;
%! ## TEST 11
%! tol = eps^0.5;
%! N = 10;
%! Phi1 = linspace(0, 90 * pi / 180, N);
%! Phi2 = linspace(0, -30 * pi / 180, N);
%! Phi3 = linspace(0, 20 * pi / 180, N);
%! x0 = [linspace(-1, 1, N);
%!       linspace(0, 3, N);
%!       linspace(-3, 0, N)].';
%! a = -0.1;
%! b = 0.8;
%! c = -0.7;
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
%! lambda = zeros(30, N);
%! for i=1:N
%!   R1 = euler123_to_rotation_matrix([Phi1(i); Phi2(i); Phi3(i)]);
%!   mesh.nodes = [(xi + x0(i, :)) * R1.', zeros(rows(xi), 3)];
%!   V = a * b * c / 6;
%!   mesh.elements.tet10 = int32(1:10);
%!   mesh.materials.tet10 = int32(1);
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 1;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case(1).loaded_nodes = int32([3; 3]);
%!   load_case(1).loads = [0, 0, -0.5;
%!                         0, 0, -0.5];
%!   load_case(1).g = zeros(3, 1);
%!   for j=1:3
%!     load_case(j + 1).g = zeros(3, 1);
%!     load_case(j + 1).g(j) = 1;
%!   endfor
%!   [dof_map] = fem_ass_dof_map(mesh, load_case_dof);
%!   [K, ...
%!    M, ...
%!    Mlumped, ...
%!    R, ...
%!    dm] = fem_ass_matrix(mesh, ...
%!                         dof_map, ...
%!                         [FEM_MAT_STIFFNESS, ...
%!                          FEM_MAT_MASS, ...
%!                          FEM_MAT_MASS_LUMPED, ...
%!                          FEM_VEC_LOAD_CONSISTENT, ...
%!                          FEM_SCA_TOT_MASS], ...
%!                         load_case);
%!   R2 = eye(3) + skew([Phi1(N - i + 1); Phi2(N - i + 1); Phi3(N - i + 1)]);
%!   def = mesh.nodes(:, 1:3) * R2.' - mesh.nodes(:, 1:3) + repmat(x0(N - i + 1, :), 10, 1);
%!   if (i == 1)
%!     K0 = K;
%!     M0 = M;
%!   else
%!     T = zeros(columns(K), columns(K));
%!     for j=1:rows(mesh.nodes)
%!       T((j - 1) * 3 + (1:3), (j - 1) * 3 + (1:3)) = R1;
%!     endfor
%!     assert_simple(T.' * K * T, K0, tol * norm(K0));
%!     assert_simple(T.' * M * T, M0, tol * norm(M0));
%!   endif
%!   U = zeros(30, 1);
%!   for j=1:3
%!     U(j:3:end) = def(:, j);
%!   endfor
%!   R1 = full(K * U);
%!   [Phi, lam] = eig(K);
%!   [lam, idx_lambda] = sort(diag(lam));
%!   Phi = Phi(:, idx_lambda);
%!   lambda(:, i) = lam;
%!   assert_simple(max(abs(R1)) < tol * norm(K * Phi(:, end)));
%!   assert_simple(0.5 * U.' * K * U < tol * Phi(:, end).' * K * Phi(:, end));
%!   assert_simple(rank(K), 24);
%!   assert_simple(isdefinite(K), false);
%!   assert_simple(isdefinite(M), true);
%!   assert_simple(isdefinite(Mlumped), true);
%!   assert_simple(length(find(lambda(:, i) < tol * max(lambda(:, i)))), 6);
%!   assert_simple(all(lambda(:, i) > -tol * max(lambda(:, i))));
%!   assert_simple(lambda(:, i), lambda(:, 1), tol * max(max(abs(lambda(:, 1)))));
%!   assert_simple(dm, mesh.material_data.rho * V, mesh.material_data.rho * V * tol);
%!   assert_simple(full(sum(diag(Mlumped))) / 3, dm, tol * dm);
%!   assert_simple(all(abs(full(diag(Mlumped)) - dm / 10) < tol * dm));
%!   for j=1:3
%!     assert_simple(full(sum(R(:, j + 1), 1)), dm, tol * dm);
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
