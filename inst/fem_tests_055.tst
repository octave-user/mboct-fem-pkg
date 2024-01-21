## fem_tests.m:55
%!test
%! ## TEST 55
%! state = rand("state");
%! unwind_protect
%!   rand("seed", 0);
%!   for i=1:30
%!     for L=2000:1000:4000
%!       N = 50;
%!       w = 10;
%!       h = 50;
%!       c2 = 0.291;
%!       f2 = -1.25;
%!       E = 210000;
%!       nu = 0.3;
%!       G = E / (2 * (1 + nu));
%!       rho = 7850e-12;
%!       A = w * h;
%!       Ay = 5 / 6 * w * h;
%!       Az = 5 / 6 * w * h;
%!       It = c2 * h * w^3;;
%!       Iy = w * h^3 / 12;
%!       Iz = w^3 * h / 12;
%!       gz = -9.81;
%!       qz = rho * A * gz;
%!       z = L - linspace(0, L, N);
%!       wz = qz * L^4 / (24 * E * Iy) * (3 - 4 * z / L + (z / L).^4);
%!       R = euler123_to_rotation_matrix(2 * pi * rand(3, 1));
%!       X = [linspace(0, L, N); zeros(2, N)];
%!       mesh.nodes = [(R * X).',  zeros(N, 3)];
%!       mesh.material_data.E = E;
%!       mesh.material_data.nu = nu;
%!       mesh.material_data.rho = rho;
%!       mesh.materials.beam2 = ones(N - 1, 1, "int32");
%!       section1.A = A;
%!       section1.Ay = Ay;
%!       section1.Az = Az;
%!       section1.It = It;
%!       section1.Iy = Iy;
%!       section1.Iz = Iz;
%!       mesh.elements.beam2 = struct("section", mat2cell(repmat(section1, N - 1, 1), ones(N - 1, 1), 1), ...
%!                                    "e2", mat2cell(repmat((R * [0; -1; 0]).', N - 1, 1), ones(N - 1, 1), 3), ...
%!                                    "nodes", mat2cell([1:N - 1; 2:N].', ones(1, N - 1), 2));
%!       load_case.locked_dof = false(size(mesh.nodes));
%!       load_case.locked_dof(1, :) = true;
%!       load_case.g = R * [0; 0; gz];
%!       [dof_map] = fem_ass_dof_map(mesh, load_case);
%!       [mat_ass.K, ...
%!        mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS, ...
%!                                     FEM_VEC_LOAD_CONSISTENT], ...
%!                                    load_case);
%!       [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!       tol = 2e-3;
%!       assert_simple(R(:, 3).' * sol_stat.def(:, 1:3).', wz, tol * max(abs(wz)));
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   rand("state", state);
%! end_unwind_protect
