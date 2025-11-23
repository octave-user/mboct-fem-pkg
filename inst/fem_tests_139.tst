%!test
%! ## Robert Gasch, Klaus Knothe, Strukturdynamik Band 1, 1987, page 132, figure 2.31a-c
%! try
%!  b = 500e-3;
%!  h = 450e-3;
%!  m = 3.5;
%!  g = 9.81;
%!  k = 1000;
%!  N = 1;
%!  mesh.nodes = [      0, 0,       0, 0, 0, 0;  #1
%!                0.5 * b, 0, 0.5 * h, 0, 0, 0;  #2
%!                      b, 0,       h, 0, 0, 0;  #3
%!                      b, 0,       0, 0, 0, 0;  #4
%!                0.5 * b, 0, 0.5 * h, 0, 0, 0;  #5
%!                      0, 0,       h, 0, 0, 0]; #6
%!  mesh.elements.line2 = [1, 2;
%!                         2, 3;
%!                         4, 5;
%!                         5, 6;
%!                         3, 6];
%!  mesh.elements.point1 = [1;2;3;4;5;6];
%!  mesh.elements.rbe2(1).nodes = [1, 2];
%!  mesh.elements.rbe2(2).nodes = [1, 3];
%!  mesh.elements.rbe2(3).nodes = [4, 5];
%!  mesh.elements.rbe2(4).nodes = [4, 6];
%!  mesh.elements.joints(1).nodes = 1;
%!  mesh.elements.joints(1).C = eye(6)([1:4, 6], :);
%!  mesh.elements.joints(2).nodes = 4;
%!  mesh.elements.joints(2).C = eye(6)([2:4, 6], :);
%!  mesh.elements.joints(3).nodes = [2, 5];
%!  mesh.elements.joints(3).C = [eye(6), -eye(6)]([1,3], :);
%!  l63 = mesh.nodes(6, 1:3).' - mesh.nodes(3, 1:3).';
%!  mesh.elements.joints(4).nodes = [3, 6];
%!  mesh.elements.joints(4).C = [eye(3), -skew(l63), -eye(3), zeros(3, 3)]([3],:);
%!  mesh.elements.joints(5).nodes = [1, 3];
%!  mesh.elements.joints(5).C = [eye(6), -eye(6)]([4,6],:);
%!  mesh.elements.bodies(1).nodes = 3;
%!  mesh.elements.bodies(1).m = m;
%!  mesh.elements.bodies(1).lcg = [-0.5 * b; 0; 0];
%!  mesh.elements.bodies(1).J = zeros(3, 3);
%!  mesh.elements.springs(1).nodes = [1, 4];
%!  mesh.elements.springs(1).K = zeros(12, 12);
%!  mesh.elements.springs(1).K([1, 7], [1, 7]) = [k, -k; -k, k];
%!  mesh.elements.springs(1).F = zeros(12,1);
%!  mesh.materials = struct();
%!  mesh.material_data = struct("E", cell(1, 0), "nu", cell(1, 0));
%!  load_case_dof.locked_dof = false(size(mesh.nodes));
%!  load_case_dof.dof_in_use = false(size(mesh.nodes));
%!  load_case_dof.dof_in_use(3, 4:6) = true;
%!  dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!  load_case.g = [0; 0; -g];
%!  [mat_ass.M, ...
%!   mat_ass.K, ...
%!   mat_ass.R, ...
%!   mat_ass.mat_info, ...
%!   mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_MASS, ...
%!                                        FEM_MAT_STIFFNESS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case);
%!  opt_sol.solver = "umfpack";
%!  sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!  opt_post.elem_types = {"line2", "point1"};
%!  opt_sol.p = 3;
%!  sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, opt_sol);
%!  lambdaref = sqrt(k / m) * sqrt(1 / (1 + (b / h)^2));
%!  fref = lambdaref / (2 * pi);
%!  tol = 1e-10;
%!  assert_simple(sol_eig.f, fref, tol * abs(fref));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
