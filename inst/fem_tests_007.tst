## fem_tests.m:07
%!test
%! try
%! ## TEST 7
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 140e-3;
%! tol = eps^0.8;
%! X = [ a,  0.5 * b,  0.5 * c;
%!       0,  0.5 * b,  0.5 * c;
%!       0, -0.5 * b,  0.5 * c;
%!       a, -0.5 * b,  0.5 * c;
%!       a,  0.5 * b, -0.5 * c;
%!       0,  0.5 * b, -0.5 * c;
%!       0, -0.5 * b, -0.5 * c;
%!       a, -0.5 * b, -0.5 * c,
%!       d,        0,        0];
%! for i=1:2
%!   if i<2
%!     idx_node = 1:8;
%!   else
%!     idx_node = 1:9;
%!   endif
%!   mesh(i).nodes = [X(idx_node, :), zeros(length(idx_node), 3)];
%!   mesh(i).elements.iso8 = int32(1:8);
%!   mesh(i).materials.iso8 = int32(1);
%!   if i>=2
%!     mesh(i).elements.rbe3.nodes = int32([9, 1, 4, 5, 8]);
%!     mesh(i).elements.rbe3.weight = ones(1, 4);
%!   endif
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh(i).material_data.rho = 1;
%!   mesh(i).material_data.C = fem_pre_mat_isotropic(E, nu);
%!   load_case(i).locked_dof = false(rows(mesh(i).nodes), 6);
%!   load_case(i).locked_dof([2; 3; 6; 7], 1:3) = true;
%!   load_case(i).loaded_nodes = int32([1; 4; 5; 8]);
%!   load_case(i).loads = [0, 0, -0.25, 0,   0, 0;
%!                         0, 0, -0.25, 0,   0, 0;
%!                         0, 0, -0.25, 0,   0, 0;
%!                         0, 0, -0.25, 0,   0, 0];
%!   dof_map{i} = fem_ass_dof_map(mesh(i), load_case(i));
%!   [mat_ass(i).K, ...
%!    mat_ass(i).M, ...
%!    mat_ass(i).dm, ...
%!    mat_ass(i).R] = fem_ass_matrix(mesh(i), ...
%!                                   dof_map{i}, ...
%!                                   [FEM_MAT_STIFFNESS, ...
%!                                    FEM_MAT_MASS, ...
%!                                    FEM_SCA_TOT_MASS, ...
%!                                    FEM_VEC_LOAD_CONSISTENT], ...
%!                                   load_case(i));
%!   sol_stat(i).def = fem_post_def_nodal(mesh(i), dof_map{i}, full(mat_ass(i).K \ mat_ass(i).R));
%!   assert_simple(mat_ass(i).dm, mesh(i).material_data.rho * a * b * c, sqrt(eps) * mesh(i).material_data.rho * a * b * c);
%! endfor
%! assert_simple(sol_stat(2).def(1:8, :), sol_stat(1).def, tol * max(max(abs(sol_stat(1).def))));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
