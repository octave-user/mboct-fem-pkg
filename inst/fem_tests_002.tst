## fem_tests.m:02
%!test
%! try
%! ## TEST 2
%! close all;
%! for iorient=1:10
%!   a = 70e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   d = 140e-3;
%!   rho = 7850;
%!   m = rho * a * b * c;
%!   xgc = 0.5 * [a; b; c];
%!   Jxx = m * (b^2 + c^2) / 12;
%!   Jyy = m * (a^2 + c^2) / 12;
%!   Jzz = m * (a^2 + b^2) / 12;
%!   scale_def = 100e-3;
%!   tol = eps^0.8;
%!   Jgc = diag([Jxx, Jyy, Jzz]);
%!   J = Jgc - m * skew(xgc) * skew(xgc);
%!   R = euler123_to_rotation_matrix((2 * rand(3, 1) - 1) * pi);
%!   X = [ 0.5 * a,  0.5 * b,  0.5 * c;
%!         0,  0.5 * b,  0.5 * c;
%!         0, -0.5 * b,  0.5 * c;
%!         0.5 * a, -0.5 * b,  0.5 * c;
%!         0.5 * a,  0.5 * b, -0.5 * c;
%!         0,  0.5 * b, -0.5 * c;
%!         0, -0.5 * b, -0.5 * c;
%!         0.5 * a, -0.5 * b, -0.5 * c,
%!         a,  0.5 * b,  0.5 * c;
%!         a, -0.5 * b,  0.5 * c;
%!         a,  0.5 * b, -0.5 * c;
%!         a, -0.5 * b, -0.5 * c,
%!         d,        0,        0;
%!         -d,        0,        0];
%!   X(:, 2) += 0.5 * b;
%!   X(:, 3) += 0.5 * c;
%!   X = X * R.';
%!   xgc = R * xgc;
%!   J = R * J * R.';
%!   mesh.nodes = [X, zeros(rows(X), 3)];
%!   mesh.elements.iso8 = int32([1:8;
%!                               9, 1, 4, 10, 11, 5, 8, 12]);
%!   mesh.materials.iso8 = int32([1; 1]);
%!   mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%!   mesh.elements.rbe3(1).weight = ones(1, 4);
%!   mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%!   mesh.elements.rbe3(2).weight = ones(1, 4);
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case(1).loaded_nodes = int32([13; 14]);
%!   load_case(1).loads = [0, 0, -0.25, 0,   0, 0;
%!                         0, 0, -0.25, 0,   0, 0];
%!   load_case(1).g = zeros(3, 1);
%!   for i=1:3
%!     load_case(i + 1).g = zeros(3, 1);
%!     load_case(i + 1).g(i) = 1;
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.Mlumped, ...
%!    mat_ass.R, ...
%!    mat_ass.dm, ...
%!    mat_ass.S, ...
%!    mat_ass.J, ...
%!    mat_ass.colloc, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS, ...
%!                                         FEM_MAT_MASS, ...
%!                                         FEM_MAT_MASS_LUMPED, ...
%!                                         FEM_VEC_LOAD_CONSISTENT, ...
%!                                         FEM_SCA_TOT_MASS, ...
%!                                         FEM_VEC_INERTIA_M1, ...
%!                                         FEM_MAT_INERTIA_J, ...
%!                                         FEM_VEC_COLL_MASS, ...
%!                                         FEM_VEC_COLL_STIFFNESS], ...
%!                                        load_case);
%!   mat_ass.D = diag(zeros(columns(mat_ass.M), 1));
%!   assert_simple(mat_ass.dm, mesh.material_data.rho * a * b * c, sqrt(eps) * mesh.material_data.rho * a * b * c);
%!   assert_simple(mat_ass.S, xgc * m, m * sqrt(eps))
%!   assert_simple(mat_ass.J, J, sqrt(eps) * norm(J));
%!   slave_dofs = [];
%!   for i=1:numel(mesh.elements.rbe3)
%!     slave_ndofs = dof_map.ndof(mesh.elements.rbe3(i).nodes(2:end), :);
%!     slave_dofs = [slave_dofs, reshape(slave_ndofs, 1, numel(slave_ndofs)), dof_map.edof.rbe3(i, :)];
%!   endfor
%!   slave_dofs = slave_dofs(find(slave_dofs > 0));
%!   master_dofs = int32(1:dof_map.totdof);
%!   master_dofs(slave_dofs) = 0;
%!   master_dofs = master_dofs(find(master_dofs > 0));
%!   assert_simple(full(sum(diag(mat_ass.Mlumped))) / 3, mat_ass.dm, tol * mat_ass.dm);
%!   for i=1:3
%!     assert_simple(full(sum(mat_ass.R(:, i + 1), 1)), m, tol * m);
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
