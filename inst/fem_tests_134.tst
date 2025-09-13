%!test
%! try
%! r = 1250e-3;
%! Phi0 = 0 * pi / 180;
%! x0 = r * cos(Phi0);
%! y0 = r * sin(Phi0);
%! zu = 100e-3;
%! zl = -50e-3;
%! k = 500;
%! mesh.nodes = [x0, y0,   0, 0, 0, 0;
%!                0,  0,  zu, 0, 0, 0;
%!                0,  0,  zl, 0, 0, 0];
%! mesh.elements.line2 = [1, 2;
%!                        1, 3;
%!                        2, 3];
%! mesh.elements.point1 = [1;2;3];
%! mesh.elements.joints(1).nodes = 2;
%! mesh.elements.joints(1).C = [eye(3), zeros(3, 3)];
%! mesh.elements.joints(2).nodes = 3;
%! mesh.elements.joints(2).C = [eye(3), zeros(3, 3)](1:2,:);
%! mesh.elements.rbe2(1).nodes = [1, 2];
%! mesh.elements.rbe2(2).nodes = [1, 3];
%! mesh.elements.springs(1).nodes = 1;
%! mesh.elements.springs(1).K = diag([0, k, zeros(1, 4)]);
%! mesh.elements.springs(1).F = zeros(6, 1);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! load_case.loads = [1, 0, 0, 0, 0, 0];
%! load_case.loaded_nodes = [1];
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
%! number_of_modes = 1;
%! opt_sol.solver = "pardiso";
%! opt_sol.problem = "buckling";
%! opt_sol.p = 23;
%! [U, sol_eig.lambda] = fem_sol_eigs(mat_ass.K, mat_ass.KTAU0, number_of_modes, opt_sol);
%! sol_eig.def = fem_post_def_nodal(mesh, dof_map, U);
%! Fref = -r * k * cos(Phi0);
%! tol = 1e-10;
%! opt_post.elem_types = {"point1", "line2"};
%! assert(sol_eig.lambda(1), Fref, tol * abs(Fref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
