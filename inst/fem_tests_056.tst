## fem_tests.m:56
%!test
%! try
%! ## TEST 56
%! l1 = 1000;
%! l2 = 800;
%! d = 10;
%! A = d^2 * pi / 4;
%! Ay = Az = 9/10 * A;
%! It = d^4 * pi / 32;
%! Iy = Iz = d^4 * pi / 64;
%! E = 210000;
%! nu = 0.3;
%! rho = 0;
%! m1 = 100;
%! gz = -9.81;
%! mesh.nodes = [0, 0, 0, 0, 0, 0;
%!               l1, 0, 0, 0, 0, 0];
%! mesh.elements.beam2.nodes = int32([1, 2]);
%! mesh.elements.beam2.material = int32(1);
%! mesh.elements.beam2.e2 = [0; 1; 0];
%! mesh.elements.beam2.section.A = A;
%! mesh.elements.beam2.section.Ay = Ay;
%! mesh.elements.beam2.section.Az = Az;
%! mesh.elements.beam2.section.It = It;
%! mesh.elements.beam2.section.Iy = Iy;
%! mesh.elements.beam2.section.Iz = Iz;
%! mesh.materials.beam2 = int32(1);
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.elements.bodies.nodes = int32(2);
%! mesh.elements.bodies.m = m1;
%! mesh.elements.bodies.J = zeros(3, 3);
%! mesh.elements.bodies.lcg = [l2; 0; 0];
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(1, 2:3) = true;
%! load_case.locked_dof(2, 1:4) = true;
%! load_case.g = [0; 0; gz];
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol = fem_sol_static(mesh, dof_map, mat_ass);
%! Phiref = -m1 * gz * l2 * l1 / (E * Iy) * [-1/6, 1/3];
%! tol = 1e-4;
%! for i=1:2
%!   assert_simple(sol.def(i, 5), Phiref(i), tol * norm(Phiref));
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
