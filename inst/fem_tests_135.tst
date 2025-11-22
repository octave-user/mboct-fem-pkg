%!test
%! try
%! ## K.J. Bathe 2002, page 113-115, example 3.11
%! L = 2350e-3;
%! k = 12500;
%! kr = 15000;
%! mesh.nodes = [0,   0,   0, 0, 0, 0;  # C
%!               0,   0,   L, 0, 0, 0;  # B
%!               0,   0, 2*L, 0, 0, 0]; # A
%! mesh.elements.line2 = [1, 2;
%!                        2, 3];
%! mesh.elements.point1 = [1;2;3];
%! mesh.elements.joints(1).nodes = 1;
%! mesh.elements.joints(1).C = eye(6)([1:4,6], :);
%! mesh.elements.joints(2).nodes = 3;
%! mesh.elements.joints(2).C = eye(6)([2,6], :);
%! mesh.elements.rbe2(1).nodes = [1, 2];
%! mesh.elements.rbe2(2).nodes = [3, 2];
%! mesh.elements.springs(1).nodes = 2;
%! mesh.elements.springs(1).K = diag([k, zeros(1, 5)]);
%! mesh.elements.springs(1).F = zeros(6, 1);
%! mesh.elements.springs(2).nodes = 3;
%! mesh.elements.springs(2).K = diag([k, zeros(1, 5)]);
%! mesh.elements.springs(2).F = zeros(6, 1);
%! mesh.elements.springs(3).nodes = [1,3];
%! mesh.elements.springs(3).K = zeros(12, 12);
%! mesh.elements.springs(3).K(5, 5) = kr;
%! mesh.elements.springs(3).K(11, 11) = kr;
%! mesh.elements.springs(3).K(5, 11) = -kr;
%! mesh.elements.springs(3).K(11, 5) = -kr;
%! mesh.elements.springs(3).F = zeros(12, 1);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! load_case.loads = [0, 0, -1, 0, 0, 0];
%! load_case.loaded_nodes = [3];
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
%! load_case_eig.lambda = sol_stat.lambda;
%! [mat_ass.KTAU0] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case_eig);
%! number_of_modes = 2;
%! opt_sol.solver = "pardiso";
%! opt_sol.problem = "buckling";
%! [U, sol_eig.lambda] = fem_sol_eigs(mat_ass.K, mat_ass.KTAU0, number_of_modes, opt_sol);
%! sol_eig.def = fem_post_def_nodal(mesh, dof_map, U);
%! opt_post.elem_types = {"line2", "point1"};
%! A = [k * L + kr / L,        -2 * kr / L;
%!         -2 * kr / L, k * L + 4 * kr / L];
%! B = [ 1, -1;
%!      -1,  2];
%! [lambdaref] = eig(A, B);
%! tol = 1e-10;
%! assert_simple(sol_eig.lambda(:), lambdaref(:), tol * norm(lambdaref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
