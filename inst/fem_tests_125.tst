## fem_tests.m:125
%!test
%! try
%! ## TEST 125
%! close all;
%! N = 20;
%! d = 1.6;
%! D = 12.45 + d;
%! L = 34.7 - 6 * d;
%! n = 6;
%! m = 40;
%! E = 206000;
%! G = 81500;
%! nu = E / (2 * G) - 1;
%! rho = 7850e-12;
%! omegai = 2 * pi * [500;
%!                    1000];
%! zetai = [0.1e-2;
%!          0.1e-2];
%! Ai = [ones(2, 1), omegai.^2];
%! Bi = 2 * omegai .* zetai;
%! Xi = Ai \ Bi;
%! alpha = Xi(1);
%! beta = Xi(2);
%! A = d^2 * pi / 4;
%! Ay = 9 / 10 * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64;
%! Iz = Iy;
%! It = Iy + Iz;
%! Phi = linspace(0, 2 * pi * n, ceil(m * n)).';
%! x = 0.5 * D * cos(Phi);
%! y = 0.5 * D * sin(Phi);
%! z = L * Phi / (2 * pi * n);
%! e2 = repmat([0, 0, 1], numel(Phi) - 1, 1);
%! mesh.nodes = [x, y, z, zeros(numel(Phi), 3)];
%! elnodes = int32([(1:numel(Phi) - 1).', (2:numel(Phi)).']);
%! section.A = A;
%! section.Ay = Ay;
%! section.Az = Az;
%! section.It = It;
%! section.Iy = Iy;
%! section.Iz = Iz;
%! mesh.elements.beam2 = struct("nodes", mat2cell(elnodes, ones(numel(Phi) - 1, 1, "int32"), 2), ...
%!                              "section", mat2cell(repmat(section, numel(Phi) - 1, 1), ones(numel(Phi) - 1, 1, "int32")), ...
%!                              "e2", mat2cell(e2, ones(numel(Phi) - 1, 1, "int32"), 3));
%! node1 = 1;
%! node2 = numel(Phi);
%! mesh.elements.joints = struct("nodes", {node1}, "C", {eye(6)});
%! U1 = zeros(6, 6);
%! load_case = struct("joints", cell(1, 6), "loads", cell(1, 6), "loaded_nodes", cell(1, 6));
%! for i=1:numel(load_case)
%!   load_case(i).joints = struct("U", {U1(:, i)});
%!   load_case(i).loaded_nodes = node2;
%!   load_case(i).loads = zeros(1, 6);
%!   load_case(i).loads(i) = 1;
%! endfor
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.material_data.alpha = alpha;
%! mesh.material_data.beta = beta;
%! mesh.materials.beam2 = ones(numel(Phi) - 1, 1, "int32");
%! load_case_dof.locked_dof = false(size(mesh.nodes));
%! [dof_map] = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_MASS, ...
%!                                       FEM_MAT_DAMPING, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case);
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N);
%! omega = linspace(0, 2 * pi * 3200, 1000);
%! U2 = F1 = complex(zeros(6, numel(omega), columns(mat_ass.R)));
%! opts.solver = "pardiso";
%! opts.refine_max_iter = int32(250);
%! opts.epsilon_refinement = eps^0.8;
%! opts.number_of_threads = mbdyn_solver_num_threads_default();
%! opts.verbose = int32(0);
%! zeta = (alpha ./ omega + beta * omega) / 2;
%! for i=1:numel(omega)
%!   U = fem_sol_factor(-omega(i)^2 * mat_ass.M + 1j * omega(i) * mat_ass.D + mat_ass.K, opts) \ mat_ass.R;
%!   for j=1:columns(U)
%!     U2(:, i, j) = U(dof_map.ndof(node2, :), j);
%!     F1(:, i, j) = U(dof_map.edof.joints(1, :), j) * mat_ass.mat_info.beta(1);
%!   endfor
%!   if (mod(i, 20) == 0)
%!     fprintf(stderr, "%d/%d\n", i, numel(omega));
%!   endif
%! endfor
%! for i=1:columns(mat_ass.R)
%!   figure("visible", "off");
%!   hold on;
%!   semilogy(omega/(2*pi), abs(U2(i, :, i)), "-;abs(U2/F2);k");
%!   semilogy(omega/(2*pi), abs(F1(i, :, i)), "-;abs(F1/F2);r");
%!   xlabel("f [Hz]");
%!   ylabel("abs(U1/U2) [1]");
%!   title({"Fx","Fy","Fz","Mx","My","Mz"}{i});
%!   grid on;
%!   grid minor on;
%! endfor
%! figure_list();
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
