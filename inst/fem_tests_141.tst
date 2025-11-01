%!test
%! try
%! ## simple pendulum
%! m = 1.2;
%! l = 2.3;
%! g = 9.81;
%! mesh.nodes = [0, 0, 0, 0, 0, 0;
%!               0, 0, -l, 0, 0, 0];
%! mesh.elements.joints(1).nodes = 1;
%! mesh.elements.joints(1).C = eye(6)(1:3, :);
%! mesh.elements.joints(2).nodes = 2;
%! mesh.elements.joints(2).C = eye(6)(6, :);
%! mesh.elements.rbe2(1).nodes = [2, 1];
%! mesh.elements.bodies(1).nodes = 2;
%! mesh.elements.bodies(1).m = m;
%! mesh.elements.bodies(1).J = zeros(3, 3);
%! mesh.elements.bodies(1).lcg = [0; 0; 0];
%! mesh.elements.springs(1).nodes = 1;
%! mesh.elements.springs(1).K = [zeros(3, 3), zeros(3, 3);
%!                               zeros(3, 3), eye(3)];
%! mesh.elements.springs(1).F = zeros(6, 1);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! load_case.g = [0; 0; -g];
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
%! load_case.lambda = sol_stat.lambda;
%! mesh.elements.springs(1).K(:, :) = 0;
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.KTAU0, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case);
%! mat_ass.K += mat_ass.KTAU0;
%! opt_sol.algorithm = "unsymmetric";
%! opt_sol.solver = "umfpack";
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 2, opt_sol);
%! fref = sqrt(g / l) / (2 * pi);
%! tol = 1e-10;
%! assert(sol_eig.f(1:2), repmat(fref, 1, 2), tol * fref);
%! assert(mat_ass.M(dof_map.ndof(2, 1:3), dof_map.ndof(2, 1:3)), m * eye(3), tol * m);
%! assert(mat_ass.K(dof_map.ndof(2, 4:5), dof_map.ndof(2, 4:5)), m * g * l * eye(2), tol * m * g * l);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
