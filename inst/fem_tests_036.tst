## fem_tests.m:36
%!test
%! try
%! ## TEST 36
%! for L=2000e-3:1000e-3:50000e-3;
%! w = 10e-3;
%! h = 50e-3;
%! c2 = 0.291;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! N = 20;
%! A = w * h;
%! Ay = 5 / 6 * w * h;
%! Az = 5 / 6 * w * h;
%! It = c2 * h * w^3;;
%! Iy = w * h^3 / 12;
%! Iz = w^3 * h / 12;
%! B = E * [Iz];
%! mu = rho * A;
%! alpha = 1e-10;
%! beta = 1e-8;
%! omega1 = sqrt(B / (mu * L^4)); ## valid only for lean beams
%! omega_ref = omega1.' * [3.516, 22.035, 61.697];
%! omega_ref = sort(omega_ref(:));
%! f_ref = omega_ref / (2 * pi);
%! R = eye(3);
%! X = [linspace(0, L, N);
%!      zeros(2, N)];
%! mesh.nodes = [(R * X).', zeros(N, 3)];
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.material_data.alpha = alpha;
%! mesh.material_data.beta = beta;
%! mesh.materials.beam2 = ones(N - 1, 1, "int32");
%! beam1.nodes = int32([]);
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
%! load_case.locked_dof(:, [1,3,4,5]) = true;
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.D, ...
%!  mat_ass.dm] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_STIFFNESS, ...
%!                                FEM_MAT_MASS, ...
%!                                FEM_MAT_DAMPING, ...
%!                                FEM_SCA_TOT_MASS], ...
%!                               load_case);
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, 3);
%! tolf = 0.1e-2;
%! tolm = 5 * eps;
%! assert_simple(max(abs(sol_eig.f(:) ./ f_ref(:) - 1)) < tolf);
%! assert_simple(mat_ass.dm, rho * A * L, tolm * rho * A * L);
%! assert_simple(mat_ass.D, mat_ass.M * alpha + mat_ass.K * beta, tolm);
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
