## fem_tests.m:35
%!test
%! try
%! ## TEST 35
%! state = rand("state");
%! unwind_protect
%! rand("seed", 0);
%! for i=1:30
%! for L = 0.1:100:1000
%! w = 10;
%! h = 50;
%! c2 = 0.291;
%! f2 = -1.25;
%! E = 210000;
%! nu = 0.3;
%! G = E / (2 * (1 + nu));
%! rho = 7850e-12;
%! A = w * h;
%! Ay = 5 / 6 * w * h;
%! Az = 5 / 6 * w * h;
%! It = c2 * h * w^3;;
%! Iy = w * h^3 / 12;
%! Iz = w^3 * h / 12;
%! U1 = zeros(3, 1);
%! Phi1 = zeros(3, 1);
%! U2 = [0; 0; f2 * L^3 / (3 * E * Iy) + f2 * L / (G * Ay)];
%! Phi2 = [0; -f2 * L^2 / (2 * E * Iy); 0];
%! R = euler123_to_rotation_matrix(2 * pi * rand(3, 1));
%! Uref = [(R * U1).', (R * Phi1).';
%!         (R * U2).', (R * Phi2).'];
%! X1 = zeros(3, 1);
%! X2 = [L; 0; 0];
%! F2 = [0; 0; f2];
%! mesh.nodes = [(R * X1).', zeros(1,3);
%!               (R * X2).',  zeros(1, 3)];
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.materials.beam2 = int32(1);
%! if (mod(i, 2))
%!   mesh.elements.beam2.nodes = int32([1, 2]);
%! else
%!   mesh.elements.beam2.nodes = int32([2, 1]);
%! endif
%! mesh.elements.beam2.section.A = A;
%! mesh.elements.beam2.section.Ay = Ay;
%! mesh.elements.beam2.section.Az = Az;
%! mesh.elements.beam2.section.It = It;
%! mesh.elements.beam2.section.Iy = Iy;
%! mesh.elements.beam2.section.Iz = Iz;
%! mesh.elements.beam2.e2 = R * [0; -1; 0];
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.locked_dof(1, :) = true;
%! load_case.loads = [(R * F2).', zeros(1, 3)];
%! load_case.loaded_nodes = int32(2);
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! assert_simple(sol_stat.def, Uref, 1e-3 * norm(Uref));
%! endfor
%! endfor
%! unwind_protect_cleanup
%! rand("state", state);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
