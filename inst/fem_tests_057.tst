## fem_tests.m:57
%!test
%! try
%! ## TEST 57
%! rndstate = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   a = 100e-3;
%!   b = 75e-3;
%!   c = 20e-3;
%!   ri = a * [0, 1/3, 2/3,  1, 2/3, 1/3,   0,   0,   0, 1/3, 2/3, 1/3,   0,   0,   0, 1/3, 1/3,   0,   0,   0];
%!   si = b * [0,   0,   0,  0, 1/3, 2/3,   1, 2/3, 1/3, 1/3,   0, 1/3, 2/3, 1/3,   0,   0,   0, 1/3,   0,   0];
%!   ti = c * [0,   0,   0,  0,   0,   0,   0,   0,   0,   0, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 2/3, 2/3, 2/3,   1];
%!   X = [ri; si; ti];
%!   g = [0; 0; -9.81];
%!   for i=1:100
%!     e1 = rand(3, 1);
%!     e2 = rand(3, 1);
%!     e3 = cross(e1, e2);
%!     e2 = cross(e3, e1);
%!     R = [e1, e2, e3];
%!     R *= diag(1 ./ norm(R, "cols"));
%!     mesh.nodes = [(R * X).', zeros(numel(ri), 3)];
%!     mesh.elements.tet20 = int32(1:20);
%!     mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!     mesh.material_data.E = 210000e6;
%!     mesh.material_data.nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     mesh.material_data.alpha = 1e-10;
%!     mesh.material_data.beta = 1e-8;
%!     load_case.locked_dof = false(size(mesh.nodes));
%!     load_case.domain = FEM_DO_STRUCTURAL;
%!     load_case.g = R * g;
%!     dof_map = fem_ass_dof_map(mesh, load_case);
%!     [mat_ass.M, ...
%!      mat_ass.D, ...
%!      mat_ass.K, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J, ...
%!      mat_ass.mat_info, ...
%!      mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_MAT_MASS, ...
%!                                           FEM_MAT_DAMPING, ...
%!                                           FEM_MAT_STIFFNESS, ...
%!                                           FEM_VEC_LOAD_CONSISTENT, ...
%!                                           FEM_SCA_TOT_MASS, ...
%!                                           FEM_VEC_INERTIA_M1, ...
%!                                           FEM_MAT_INERTIA_J], ...
%!                                          load_case);
%!     F = zeros(rows(mesh.nodes), 3);
%!     for j=1:3
%!       F(:, j) = mat_ass.R(dof_map.ndof(:, j));
%!     endfor
%!     F = sum(F, 1).';
%!     mref = mesh.material_data.rho * a * b * c / 6;
%!     tolm = eps^0.9;
%!     tollambda = eps^0.9;
%!     assert_simple(mat_ass.dm, mref, tolm * mref);
%!     assert_simple(rank(mat_ass.M), 3 * columns(mesh.elements.tet20));
%!     assert_simple(rank(mat_ass.K), 3 * columns(mesh.elements.tet20) - 6);
%!     [Phi, lambda] = eig(mat_ass.K, mat_ass.M);
%!     lambda = diag(lambda);
%!     sol.def = zeros(rows(mesh.nodes), 6, columns(Phi));
%!     for j=1:3
%!       for k=1:columns(Phi)
%!         sol.def(:, j, k) = Phi(dof_map.ndof(:, j), k);
%!       endfor
%!     endfor
%!     assert_simple(isdefinite(mat_ass.M), true);
%!     assert_simple(isdefinite(mat_ass.K), false);
%!     assert_simple(R.' * F, mref * g, tolm * norm(mref * g));
%!     if (i == 1)
%!       Rref = R;
%!       Sref = mat_ass.S;
%!       Jref = mat_ass.J;
%!       lambdaref = lambda;
%!     else
%!       assert_simple(lambda, lambdaref, tollambda * max(abs(lambdaref)));
%!       assert_simple(R.' * mat_ass.S, Rref.' * Sref, tolm * norm(Sref));
%!       assert_simple(R.' * mat_ass.J * R, Rref.' * Jref * Rref, tolm * norm(Jref));
%!     endif
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", rndstate);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
