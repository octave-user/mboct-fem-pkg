## fem_tests.m:46
%!test
%! try
%! ## TEST 46
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   a = 50e-3;
%!   b = 20e-3;
%!   c = 15e-3;
%!   x = a * [  1, -1, -1,  1,  1, -1, -1,  1, 0, -1,  0, 1,  0, -1,  0,  1, 1, -1, -1,  1] + 2 * a;
%!   y = b * [  1,  1, -1, -1,  1,  1, -1, -1, 1,  0, -1, 0,  1,  0, -1,  0, 1,  1, -1, -1] + 2 * b;
%!   z = c * [  1,  1,  1,  1, -1, -1, -1, -1, 1,  1,  1, 1, -1, -1, -1, -1, 0,  0,  0,  0] + 2 * c;
%!   E = 210000e6;
%!   nu = 0.3;
%!   for j=1:100
%!     for k=1:6
%!       switch (k)
%!       case {1, 2, 3}
%!         epsilon = [2 * rand(3, 1) - 1; zeros(3, 1)];
%!       otherwise
%!         epsilon = zeros(6, 1);
%!         epsilon(k) = 2 * rand() - 1;
%!       endswitch
%!       e1 = rand(3, 1);
%!       e2 = rand(3, 1);
%!       e3 = cross(e1, e2);
%!       e2 = cross(e3, e1);
%!       e1 /= norm(e1);
%!       e2 /= norm(e2);
%!       e3 /= norm(e3);
%!       R = [e1, e2, e3];
%!       mesh.nodes = [(R * [x; y; z]).', zeros(numel(x), 3)];
%!       mesh.elements.iso20r = int32(1:rows(mesh.nodes));
%!       mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!       mesh.material_data.rho = 7850;
%!       mesh.materials.iso20r = int32(1);
%!       load_case.locked_dof = false(size(mesh.nodes));
%!       dof_map = fem_ass_dof_map(mesh, load_case);
%!       sol_stat.def = zeros(size(mesh.nodes));
%!       gammaxx = epsilon(1);
%!       gammayy = epsilon(2);
%!       gammazz = epsilon(3);
%!       gammaxy = epsilon(4);
%!       gammayz = epsilon(5);
%!       gammazx = epsilon(6);
%!       X = mesh.nodes(:, 1);
%!       Y = mesh.nodes(:, 2);
%!       Z = mesh.nodes(:, 3);
%!       u = (X.*Z.*gammazx-Y.*Z.*gammayz+X.*Y.*gammaxy)./(2.*X) + gammaxx * X;
%!       v = ((-X.*Z.*gammazx)+Y.*Z.*gammayz+X.*Y.*gammaxy)./(2.*Y) + gammayy * Y;
%!       w = -((-X.*Z.*gammazx)-Y.*Z.*gammayz+X.*Y.*gammaxy)./(2.*Z) + gammazz * Z;
%!       sol_stat.def(:, 1) = u;
%!       sol_stat.def(:, 2) = v;
%!       sol_stat.def(:, 3) = w;
%!       tau = mesh.material_data.C * epsilon;
%!       sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        FEM_VEC_STRESS_CAUCH, ...
%!                                        load_case, ...
%!                                        sol_stat);
%!      switch (k)
%!      case {1,2,3}
%!        for i=1:3
%!          assert_simple(max(abs(sol_stat.stress.tau.iso20r(1, :, i) ./ tau(i) - 1)) < sqrt(eps));
%!          assert_simple(max(abs(sol_stat.stress.tau.iso20r(1, :, i + 3) ./ max(abs(tau)))) < sqrt(eps));
%!        endfor
%!      otherwise
%!        assert_simple(max(abs(sol_stat.stress.tau.iso20r(1, :, k) ./ tau(k) - 1)) < sqrt(eps));
%!      endswitch
%!   endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
