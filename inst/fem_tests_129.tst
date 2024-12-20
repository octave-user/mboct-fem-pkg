## fem_tests.m:129
%!test
%! try
%! ## TEST 129
%! m1 = 1.5;
%! m2 = 30;
%! m3 = 0.01;
%! J1 = zeros(3, 3);
%! J2 = zeros(3, 3);
%! s1 = 2000;
%! s2 = 10000;
%! s3 = 50;
%! d1 = 50;
%! d2 = 800;
%! d3 = 0.15;
%! f1 = 1000;
%! f2 = 2000;
%! f3 = 3000;
%! mesh.nodes = zeros(3, 6);
%! mesh.materials = struct();
%! mesh.material_data = struct()([]);
%! mesh.elements.bodies(1).nodes = int32(1);
%! mesh.elements.bodies(1).m = m1;
%! mesh.elements.bodies(1).J = J1;
%! mesh.elements.bodies(1).lcg = zeros(3, 1);
%! mesh.elements.bodies(2).nodes = int32(2);
%! mesh.elements.bodies(2).m = m2;
%! mesh.elements.bodies(2).J = J2;
%! mesh.elements.bodies(2).lcg = zeros(3, 1);
%! mesh.elements.bodies(3).nodes = int32(3);
%! mesh.elements.bodies(3).m = 0;
%! mesh.elements.bodies(3).J = diag([0, 0, m3]);
%! mesh.elements.bodies(3).lcg = zeros(3, 1);
%! mesh.elements.springs(1).nodes = int32(1);
%! mesh.elements.springs(1).K = diag([s1, zeros(1, 5)]);
%! mesh.elements.springs(1).F = [f1; zeros(5, 1)];
%! mesh.elements.dashpots(1).nodes = int32(1);
%! mesh.elements.dashpots(1).D = diag([d1, zeros(1, 5)]);
%! mesh.elements.springs(2).nodes = int32(2);
%! mesh.elements.springs(2).K = diag([zeros(1, 2), s2, zeros(1, 3)]);
%! mesh.elements.springs(2).F = [zeros(2, 1); f2; zeros(3, 1)];
%! mesh.elements.dashpots(2).nodes = int32(2);
%! mesh.elements.dashpots(2).D = diag([zeros(1, 2), d2, zeros(1, 3)]);
%! mesh.elements.springs(3).nodes = int32(3);
%! mesh.elements.springs(3).K = diag([zeros(1, 5), s3]);
%! mesh.elements.springs(3).F = [zeros(5, 1); f3];
%! mesh.elements.dashpots(3).nodes = int32(3);
%! mesh.elements.dashpots(3).D = diag([zeros(1, 5), d3]);
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.locked_dof(1, 2:6) = true;
%! load_case.locked_dof(2, [1:2, 4:6]) = true;
%! load_case.locked_dof(3, 1:5) = true;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_info, ...
%!  mesh_info] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_MASS, ...
%!                               FEM_MAT_STIFFNESS, ...
%!                               FEM_MAT_DAMPING, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%!  [sol, U] = fem_sol_static(mesh, dof_map, mat_ass);
%!  U_ref = [f1 / s1; f2 / s2; f3 / s3];
%!  R_ref = [f1; f2; f3];
%!  C = [ -mat_ass.M \ mat_ass.D,           -mat_ass.M \ mat_ass.K;
%!         eye(columns(mat_ass.M)),   zeros(size(mat_ass.M))];
%!  [R, lambda] = eig(C);
%!  lambda = diag(lambda);
%!  [~, idx] = sort(imag(lambda), "ascend");
%!  lambda = lambda(idx);
%!  omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!  omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!  omega3 = sqrt(s3 / m3 - (d3 / (2 * m3))^2);
%!  alpha1 = -d1 / (2 * m1);
%!  alpha2 = -d2 / (2 * m2);
%!  alpha3 = -d3 / (2 * m3);
%!  lambda1 = alpha1 + [-1j; 1j] * omega1;
%!  lambda2 = alpha2 + [-1j; 1j] * omega2;
%!  lambda3 = alpha3 + [-1j; 1j] * omega3;
%!  lambda_ref = [lambda1; lambda2; lambda3];
%!  [~, idx] = sort(imag(lambda_ref), "ascend");
%!  lambda_ref = lambda_ref(idx);
%!  tol = eps^0.9;
%!  assert_simple(max(abs(lambda ./ lambda_ref - 1)) < tol);
%!  assert_simple(max(abs(U - U_ref)) /norm(U_ref) < tol);
%!  assert_simple(max(abs(full(mat_ass.R) - R_ref)) / norm(R_ref) < tol); 
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
