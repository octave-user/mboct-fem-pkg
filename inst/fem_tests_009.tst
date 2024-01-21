## fem_tests.m:09
%!test
%! ## TEST 9
%! X = [1,  1,  1;
%!      -1,  1,  1;
%!      -1, -1,  1;
%!      1, -1,  1;
%!      1,  1, -1;
%!      -1,  1, -1;
%!      -1, -1, -1;
%!      1, -1, -1];
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32(1:8);
%! mesh.materials.iso8 = int32(1);
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 1;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [K, M, dm] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS, FEM_MAT_MASS, FEM_SCA_TOT_MASS]);
%! assert_simple(isdefinite(K), false);
%! assert_simple(isdefinite(M), true);
%! assert_simple(dm, mesh.material_data.rho * 8, sqrt(eps) * mesh.material_data.rho * 8);
