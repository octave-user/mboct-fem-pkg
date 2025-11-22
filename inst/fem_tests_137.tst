%!test
%! try
%! ## Robert Gasch, Klaus Knothe, Strukturdynamik Band1, 1987, page 340-343, chapter 7.4.6
%! a = 1750e-3;
%! m = 8.8;
%! J = 2 / 3 * m * a^2;
%! g = 9.81;
%! G = m * g;
%! mesh.nodes = [    0, 0,     0, 0, 0, 0; ## 1: body
%!                   a, 0,     a, 0, 0, 0; ## 2: hinge 1
%!                  -a, 0,     a, 0, 0, 0; ## 3: hinge 2
%!               3 * a, 0, 3 * a, 0, 0, 0; ## 4: clamp 1
%!              -3 * a, 0, 3 * a, 0, 0, 0];## 5: clamp 2
%! mesh.elements.line2 = [1, 2;
%!                        1, 3;
%!                        2, 4;
%!                        3, 5];
%! mesh.elements.point1 = [1;2;3;4;5];
%! mesh.elements.joints(1).nodes = 4;
%! mesh.elements.joints(1).C = eye(6)([1:3],:);
%! mesh.elements.joints(2).nodes = 5;
%! mesh.elements.joints(2).C = eye(6)([1:3],:);
%! mesh.elements.joints(3).nodes = 4;
%! mesh.elements.joints(3).C = [0,0,0,1,0,1];
%! mesh.elements.joints(4).nodes = 5;
%! mesh.elements.joints(4).C = [0,0,0,-1,0,1];
%! mesh.elements.springs(1).nodes = 1;
%! mesh.elements.springs(1).K = diag(repmat(1, 1, 6));
%! mesh.elements.springs(1).F = zeros(6, 1);
%! mesh.elements.rbe2(1).nodes = [1, 2];
%! mesh.elements.rbe2(2).nodes = [1, 3];
%! mesh.elements.rbe2(3).nodes = [4, 2];
%! mesh.elements.rbe2(4).nodes = [5, 3];
%! mesh.elements.bodies(1).nodes = 1;
%! mesh.elements.bodies(1).m = m;
%! mesh.elements.bodies(1).J = diag(repmat(J, 1, 3));
%! mesh.elements.bodies(1).lcg = zeros(3, 1);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! load_case.g = [0; 0; -g];
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! opt_sol.solver = "umfpack";
%! opt_sol.refine_max_iter = 100;
%! opt_sol.pre_scaling = true;
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%! load_case_eig.lambda = sol_stat.lambda;
%! mesh.elements.springs(1).K(:, :) = 0;
%! [mat_ass.K, ...
%!  mat_ass.KTAU0, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case_eig);
%! number_of_modes = 4;
%! mat_ass.K += mat_ass.KTAU0;
%! opt_sol.p = 6;
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes, opt_sol);
%! opt_post.elem_types = {"point1", "line2"};
%!
%! Mref = diag([m, m, m, J, J, J, 0, 0]);
%!
%! Kref = 1/4 * [       G / a,         0,            0,         0,      2 * G,            0, -2 * sqrt(2), -2 * sqrt(2);
%!                          0, 2 * G / a,            0,    -2 * G,          0,            0,            0,            0;
%!                          0,         0,        G / a,         0,          0,            0, -2 * sqrt(2), -2 * sqrt(2);
%!                          0,    -2 * G,            0, 6 * G * a,          0,            0,            0,            0;
%!                      2 * G,         0,            0,         0, 12 * G * a,            0,            0,            0;
%!                          0,         0,            0,         0,          0,    6 * G * a,            0,            0;
%!               -2 * sqrt(2),         0, -2 * sqrt(2),         0,          0,            0,            0,            0;
%!               -2 * sqrt(2),         0, -2 * sqrt(2),         0,          0,            0,            0,            0];
%! Mref = Mref([2,4:6], [2,4:6]);
%! Kref = Kref([2,4:6], [2,4:6]);
%! lambdaref = eig(Kref, Mref);
%! lambdaref = sort(lambdaref(isfinite(lambdaref)));
%! lambdaref = sqrt(-lambdaref);
%! fref = imag(lambdaref)/(2*pi);
%! tol = 1e-8;
%! assert_simple(sol_eig.f(:), fref(1:numel(sol_eig.f)), tol * norm(fref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
