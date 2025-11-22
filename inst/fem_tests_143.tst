%!test
%! try
%! U_pre = [1.1; 2.2; 3.3];
%! mesh.nodes = [0,   0,   0, 0, 0, 0;
%!               1,   2,   3, 0, 0, 0;];
%! mesh.elements.joints(1).nodes = 1;
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.rbe2(1).nodes = [1, 2];
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! load_case(1).rbe2(1).U = -U_pre;
%! load_case(2).rbe2(1).U = -2 * U_pre;
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! tol = eps^0.9;
%! assert_simple(sol_stat.def(2, 1:3, 1).', U_pre, tol * norm(U_pre));
%! assert_simple(sol_stat.def(2, 1:3, 2).', 2 * U_pre, tol * norm(U_pre));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
