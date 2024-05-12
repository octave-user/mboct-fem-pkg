## fem_tests.m:08
%!test
%! try
%! ## TEST 8
%! X = [ 1,  1,  1;
%!       -1,  1,  1;
%!       -1, -1,  1;
%!       1, -1,  1;
%!       1,  1, -1;
%!       -1,  1, -1;
%!       -1, -1, -1;
%!       1, -1, -1];
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32(1:8);
%! mesh.materials.iso8 = int32(1);
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 1;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof([2; 3; 6; 7], :) = true;
%! load_case.loaded_nodes = int32([1; 4; 5; 8]);
%! load_case.loads = [0, 0, -0.25;
%!                    0, 0, -0.25;
%!                    0, 0, -0.25;
%!                    0, 0, -0.25];
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [K, M, R, dm] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS, FEM_MAT_MASS, FEM_VEC_LOAD_CONSISTENT, FEM_SCA_TOT_MASS], load_case);
%! assert_simple(isdefinite(K), true);
%! assert_simple(isdefinite(M), true);
%! assert_simple(dm, mesh.material_data.rho * 8, sqrt(eps) * mesh.material_data.rho * 8);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
