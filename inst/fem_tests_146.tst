%!test
%! ## Boege, Aufgabensammlung Mechanik, 1995
%! ## page 8, example 55
%! try
%! l1 = 3;
%! l2 = 1.5;
%! l3 = 4;
%! F = 20e3;
%! mesh.nodes = [l3, 0, l1 + l2, 0, 0, 0;
%!                0, 0,      l1, 0, 0, 0;
%!                0, 0,       0, 0, 0, 0];
%! mesh.elements.joints(1).nodes = int32([1, 2]);
%! mesh.elements.joints(1).C = [l3, 0, l2, 0, 0, 0, -l3, 0, -l2, 0, 0, 0] / sqrt(l2^2 + l3^2);
%! mesh.elements.joints(2).nodes = int32([1, 3]);
%! mesh.elements.joints(2).C = [l3, 0, l1 + l2, 0, 0, 0, -l3, 0, -l1 - l2, 0, 0, 0] / sqrt((l1 + l2)^2 + l3^2);
%! mesh.elements.joints(3).nodes = int32(2);
%! mesh.elements.joints(3).C = [eye(3), zeros(3, 3)]([1, 3], :);
%! mesh.elements.joints(4).nodes = int32(3);
%! mesh.elements.joints(4).C = [eye(3), zeros(3, 3)]([1, 3], :);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! load_case_dof.locked_dof(:, [2, 4:6]) = true;
%! load_case_dof.dof_in_use = false(size(mesh.nodes));
%! load_case_dof.dof_in_use(:, [1, 3]) = true;
%! load_case.loaded_nodes = int32(1);
%! load_case.loads = [0, 0, -F, 0, 0, 0];
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, mat_ass.R, mat_ass.mat_info, mat_ass.mesh_info] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! F_Z = sol_stat.lambda.joints(1, 1);
%! F_D = sol_stat.lambda.joints(2, 1);
%! F_Zx = sol_stat.lambda.joints(3, 1);
%! F_Zz = sol_stat.lambda.joints(3, 2);
%! F_Dx = sol_stat.lambda.joints(4, 1);
%! F_Dz = sol_stat.lambda.joints(4, 2);
%! tol = 1e-3 * F;
%! assert_simple(F_Z, 28.48e3, tol);
%! assert_simple(F_D, -40.14e3, tol);
%! assert_simple(F_Zx, 26.67e3, tol);
%! assert_simple(F_Zz, 10e3, tol);
%! assert_simple(F_Dx, -26.67e3, tol);
%! assert_simple(F_Dz, -30e3, tol);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
