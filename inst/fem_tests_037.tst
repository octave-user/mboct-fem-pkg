## fem_tests.m:37
%!test
%! ## TEST 37
%! for L = 8000e-3:1000e-3:50000e-3;
%! w = 10e-3;
%! h = 50e-3;
%! c2 = 0.291;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! N = 50;
%! A = w * h;
%! Ay = 5 / 6 * w * h;
%! Az = 5 / 6 * w * h;
%! It = c2 * h * w^3;;
%! Iy = w * h^3 / 12;
%! Iz = w^3 * h / 12;
%! B = E * [Iy];
%! mu = rho * A;
%! omega1 = sqrt(B / (mu * L^4));
%! omega_ref = omega1.' * [3.516, 22.035, 61.697]; ## valid only for lean beams
%! omega_ref = sort(omega_ref(:));
%! f_ref = omega_ref / (2 * pi);
%! R = eye(3);
%! X = [linspace(0, L, N);
%!      zeros(2, N)];
%! mesh.nodes = [(R * X).', zeros(N, 3)];
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.materials.beam2 = ones(N - 1, 1, "int32");
%! beam1.nodes = int32([]);
%! beam1.material = int32(1);
%! beam1.section.A = A;
%! beam1.section.Ay = Ay;
%! beam1.section.Az = Az;
%! beam1.section.It = It;
%! beam1.section.Iy = Iy;
%! beam1.section.Iz = Iz;
%! beam1.e2 = R * [0; 1; 0];
%! mesh.elements.beam2 = repmat(beam1, N - 1, 1);
%! for i=1:N - 1
%!   mesh.elements.beam2(i).nodes = int32(i:i+1);
%! endfor
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.locked_dof(1, :) = true;
%! load_case.locked_dof(:, [1,2,4,6]) = true;
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_MAT_MASS], ...
%!                              load_case);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, 3);
%! assert_simple(max(abs(sol_eig.f(:) ./ f_ref(:) - 1)) < 0.1e-2);
%! endfor
