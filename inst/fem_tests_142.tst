%!test
%! ## double pendulum
%! m1 = 1.2;
%! l1 = 2.3;
%! m2 = 0.5;
%! l2 = 1.3;
%! g = 9.81;
%! N = 4;
%! mesh.nodes = [0, 0,        0, 0, 0, 0;
%!               0, 0,      -l1, 0, 0, 0;
%!               0, 0, -l1 - l2, 0, 0, 0];
%! mesh.elements.line2 = [1, 2;
%!                        2, 3];
%! mesh.elements.bodies(1).nodes = 1;
%! mesh.elements.bodies(1).m = m1;
%! mesh.elements.bodies(1).J = zeros(3, 3);
%! mesh.elements.bodies(1).lcg = [0; 0; -l1];
%! mesh.elements.bodies(2).nodes = 2;
%! mesh.elements.bodies(2).m = m2;
%! mesh.elements.bodies(2).J = zeros(3, 3);
%! mesh.elements.bodies(2).lcg = [0; 0; -l2];
%! mesh.elements.rbe2(1).nodes = [1, 2];
%! mesh.elements.rbe2(2).nodes = [2, 3];
%! mesh.elements.joints(1).nodes = 1;
%! mesh.elements.joints(1).C = eye(6)([1:3, 6], :);
%! mesh.elements.joints(2).nodes = 2;
%! mesh.elements.joints(2).C = eye(6)(6, :);
%! mesh.elements.springs(1).nodes = 1;
%! mesh.elements.springs(1).K = [zeros(3, 6);
%!                               zeros(3, 3), eye(3)];
%! mesh.elements.springs(1).F = zeros(6, 1);
%! mesh.elements.springs(2).nodes = 2;
%! mesh.elements.springs(2).K = [zeros(3, 6);
%!                               zeros(3, 3), eye(3)];
%! mesh.elements.springs(2).F = zeros(6, 1);
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
%! mesh.elements.springs(1).K(:, :) = 0;
%! mesh.elements.springs(2).K(:, :) = 0;
%! load_case.lambda = sol_stat.lambda;
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.KTAU0, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_STIFFNESS_TAU0, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! mat_ass.K += mat_ass.KTAU0;
%! opt_sol.algorithm = "shift-invert";
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, opt_sol);
%! opt_post.elem_types = {"line2", "point1"};
%! M = [(m1 + m2) * l1^2, l1 * l2 * m2;
%!      l1 * l2 * m2, l2^2 * m2];
%! K = [l1 * (m1 + m2) * g, 0;
%!      0, l2 * m2 * g];
%! lambdaref = eig(K, M);
%! fref = sqrt(lambdaref) / (2 * pi);
%! tol = 1e-10;
%! assert(sol_eig.f(1), fref(1), tol * fref(1));
%! assert(sol_eig.f(2), fref(1), tol * fref(1));
%! assert(sol_eig.f(3), fref(2), tol * fref(2));
%! assert(sol_eig.f(4), fref(2), tol * fref(2));
