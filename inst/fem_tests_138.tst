%!test
%! try
%! ## Robert Gasch, Klaus Knothe, Strukturdynamik Band1, 1987, page 340-343, chapter 7.4.6
%! a = 1750e-3;
%! m = 8.8;
%! J = 2 / 3 * m * a^2;
%! g = 9.81;
%! G = m * g;
%! mesh_data(1).mesh.nodes = [    0, 0,     0, 0, 0, 0; ## 1: body
%!                                a, 0,     a, 0, 0, 0; ## 2: hinge 1
%!                               -a, 0,     a, 0, 0, 0; ## 3: hinge 2
%!                            3 * a, 0, 3 * a, 0, 0, 0; ## 4: clamp 1
%!                           -3 * a, 0, 3 * a, 0, 0, 0; ## 5: clamp 2
%!                                0, 0,   0.5, 0, 0, 0];
%! mesh_data(1).mesh.elements.line2 = [1, 2;
%!                                     1, 3;
%!                                     2, 4;
%!                                     3, 5];
%! mesh_data(1).mesh.elements.point1 = [1;2;3;4;5];
%! mesh_data(1).mesh.elements.joints(1).nodes = 4;
%! mesh_data(1).mesh.elements.joints(1).C = eye(6)([1:3],:);
%! mesh_data(1).mesh.elements.joints(2).nodes = 5;
%! mesh_data(1).mesh.elements.joints(2).C = eye(6)([1:3],:);
%! mesh_data(1).mesh.elements.joints(3).nodes = 4;
%! mesh_data(1).mesh.elements.joints(3).C = [0,0,0,1,0,1];
%! mesh_data(1).mesh.elements.joints(4).nodes = 5;
%! mesh_data(1).mesh.elements.joints(4).C = [0,0,0,-1,0,1];
%! mesh_data(1).mesh.elements.rbe2(1).nodes = [1, 2];
%! mesh_data(1).mesh.elements.rbe2(2).nodes = [1, 3];
%! mesh_data(1).mesh.elements.rbe2(3).nodes = [4, 2];
%! mesh_data(1).mesh.elements.rbe2(4).nodes = [5, 3];
%! mesh_data(1).mesh.elements.rbe3(1).nodes = [6, 1, 2, 3, 4, 5];
%! mesh_data(1).mesh.elements.rbe3(1).weight = ones(1, 5);
%! mesh_data(1).mesh.elements.bodies(1).nodes = 1;
%! mesh_data(1).mesh.elements.bodies(1).m = m;
%! mesh_data(1).mesh.elements.bodies(1).J = diag(repmat(J, 1, 3));
%! mesh_data(1).mesh.elements.bodies(1).lcg = zeros(3, 1);
%! mesh_data(1).mesh.materials = struct();
%! mesh_data(1).mesh.material_data = struct()([]);
%! mesh_data(1).load_case_dof.locked_dof = false(size(mesh_data(1).mesh.nodes));
%! mesh_data(1).dof_map = fem_ass_dof_map(mesh_data(1).mesh, mesh_data(1).load_case_dof);
%! mesh_data = repmat(mesh_data, 1, 2);
%! mesh_data(2).mesh.nodes(:, 2) += 10000e-3;
%! [mesh, dof_map] = fem_post_mesh_merge(mesh_data, struct());
%! load_case.g = [0; 0; -g];
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
%! opt_sol.solver = "pardiso";
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%! opt_post.elem_types = {"line2", "point1"};
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
