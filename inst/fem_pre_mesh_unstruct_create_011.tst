## fem_pre_mesh_unstruct_create.m:11
%!test
%! try
%! ## TEST11
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.m1 = 2.2 / SI_unit_kilogram;
%! param.J1 = 1e-6 * [2.4427674e+03 -3.2393301e+01 -5.3623318e+02
%!                    -3.2393301e+01  3.7506484e+03  7.9745924e+01
%!                    -5.3623318e+02  7.9745924e+01  4.2943071e+03] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg1 = zeros(3, 1);
%! param.m2 = 2.2978223 / SI_unit_kilogram;
%! param.J2 = 1e-6 * [8.6104517e+03, 5.5814351e+01, -3.0103453e+02;
%!                    5.5814351e+01, 1.1905289e+04,  2.0425595e+02;
%!                    -3.0103453e+02, 2.0425595e+02,  1.2109596e+04] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg2 = 1e-3 * [12.940131;
%!                      -6.0192164e-01;
%!                      5.9683120e+01] / SI_unit_meter;
%! param.offset1 = (6.63e-3 + 2.4e-3) / SI_unit_meter;
%! param.N = 200;
%! param.maxdef = 10e-3 / SI_unit_meter;
%! param.fmin = 0;
%! param.fmax = 15000/ (SI_unit_second^-1);
%! param.num_freq_modal = 5000;
%! param.num_freq_nodal = 0;
%! param.damping = "local";
%! param.solver = "undamped";
%! switch (param.solver)
%! case "damped"
%!   param.N *= 2;
%! endswitch
%! damp.D_helspr1 = [1e-2; 0.03e-2];
%! damp.f_helspr1 = [20, 1000] / (SI_unit_second^-1);
%! damp.D_rubspr1 = [10e-2; 30e-2];
%! damp.f_rubspr1 = [30, 100] / (SI_unit_second^-1);
%! damp.D_global = [1e-3; 1e-3];
%! damp.f_global = [10, 10000] / (SI_unit_second^-1);
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.D = 12.12e-3 / SI_unit_meter + helspr1.d;
%! helspr1.L = 24.2e-3 / SI_unit_meter;
%! helspr1.n = 5.7;
%! helspr1.ni = 2.7;
%! helspr1.ng = 0.75;
%! helspr1.m = 200;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! [helspr1.material.alpha, helspr1.material.beta] = fem_pre_mat_rayleigh_damping(damp.D_helspr1, damp.f_helspr1);
%! helspr1.material.nu = helspr1.material.E / (2 * helspr1.material.G) - 1;
%! helspr1.material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! helspr1.X = [[45.20e-3,  45.20e-3, -62.5e-3;
%!               43.13e-3, -43.13e-3,   0.0e-3] / SI_unit_meter;
%!              [-helspr1.L, -helspr1.L, -helspr1.L]];
%! section1.A = helspr1.d^2 * pi / 4;
%! section1.Ay = 9 / 10 * section1.A;
%! section1.Az = section1.Ay;
%! section1.Iy = helspr1.d^4 * pi / 64;
%! section1.Iz = section1.Iy;
%! section1.It = section1.Iy + section1.Iz;
%! rubspr1.d = 11.1e-3 / SI_unit_meter;
%! rubspr1.D = 28.6e-3 / SI_unit_meter;
%! rubspr1.L = 10.5e-3 / SI_unit_meter;
%! rubspr1.material.E = 2.6e6 / SI_unit_pascal;
%! rubspr1.material.nu = 0.499;
%! [rubspr1.material.alpha, rubspr1.material.beta] = fem_pre_mat_rayleigh_damping(damp.D_rubspr1, damp.f_rubspr1);
%! rubspr1.material.rho = 910 / (SI_unit_kilogram / SI_unit_meter^3);
%! rubspr1.X = [(170e-3 * [0.5, 0.5, -0.5, -0.5] + 8e-3) / SI_unit_meter;
%!              70e-3 * [-0.5, 0.5, -0.5, 0.5] / SI_unit_meter;
%!              -[1, 1, 1, 1] * (helspr1.L + rubspr1.L + param.offset1)];
%! rubspr1.e2 = [1, 0, 0];
%! section2.A = (rubspr1.D^2 - rubspr1.d^2) * pi / 4;
%! section2.Ay = 0.8 * section2.A;
%! section2.Az = section2.Ay;
%! section2.Iy = (rubspr1.D^4 - rubspr1.d^4) * pi / 64;
%! section2.Iz = section2.Iy;
%! section2.It = section2.Iy + section2.Iz;
%! node_idx_rb1 = int32(1);
%! node_idx_rb2 = int32(2);
%! mesh.nodes(node_idx_rb1, :) = [0, 0, rubspr1.L + param.offset1 + helspr1.L, zeros(1, 3)];
%! mesh.nodes(node_idx_rb2, :) = [0, 0, rubspr1.L, zeros(1, 3)];
%! node_idx_helspr1 = int32(rows(mesh.nodes) + 1);
%! helspr1.Phi = linspace(0, 2 * pi * helspr1.n, ceil(helspr1.m * helspr1.n)).' + 2 * pi * helspr1.ni;
%! node_idx_rubspr1 = node_idx_helspr1 + numel(helspr1.Phi) * columns(helspr1.X);
%! helspr1.x = 0.5 * helspr1.D * cos(helspr1.Phi);
%! helspr1.y = 0.5 * helspr1.D * sin(helspr1.Phi);
%! helspr1.z = (helspr1.L - helspr1.d * (2 * (helspr1.ni - helspr1.ng) + 1)) * linspace(0, 1, numel(helspr1.Phi))(:) + helspr1.d * (helspr1.ni - helspr1.ng + 0.5);
%! helspr1.e2 = [0, 0, 1];
%! helspr1.Phi0 = [pi/2, -pi/2, pi] - 2 * pi * helspr1.ni;
%! idx_beam = int32(0);
%! for i=1:columns(helspr1.X)
%!   R = euler123_to_rotation_matrix([0; 0; helspr1.Phi0(i)]);
%!   idx_node = node_idx_helspr1 - 1 + (1:numel(helspr1.Phi)) + (i - 1) * numel(helspr1.Phi);
%!   elnodes = int32([(1:numel(helspr1.Phi) - 1).', (2:numel(helspr1.Phi)).']) + node_idx_helspr1 - 1 + (i - 1) * numel(helspr1.Phi);
%!   mesh.nodes(idx_node, 1:6) = [[helspr1.x, helspr1.y, helspr1.z] * R.' + helspr1.X(:, i).', zeros(numel(helspr1.Phi), 3)];
%!   mesh.elements.beam2(idx_beam + (1:numel(helspr1.Phi) - 1)) = struct("nodes", mat2cell(elnodes, ones(numel(helspr1.Phi) - 1, 1, "int32"), 2), ...
%!                                                                       "section", mat2cell(repmat(section1, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32")), ...
%!                                                                       "e2", mat2cell(repmat(helspr1.e2, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32"), 3));
%!   idx_beam += numel(helspr1.Phi) - 1;
%! endfor
%! for i=1:columns(rubspr1.X)
%!   elnodes = int32([1, 2]) + node_idx_rubspr1 - 1 + (i - 1) * 2;
%!   mesh.nodes(node_idx_rubspr1 - 1 + (1:2) + 2 * (i - 1), 1:6) = [[zeros(1, 3); zeros(1, 2), rubspr1.L] + rubspr1.X(:, i).', zeros(2, 3)];
%!   mesh.elements.beam2(idx_beam + 1) = struct("nodes", mat2cell(elnodes, ones(1, 1, "int32"), 2), ...
%!                                              "section", mat2cell(repmat(section2, 1, 1), ones(1, 1, "int32")), ...
%!                                              "e2", mat2cell(rubspr1.e2, ones(1, 1, "int32"), 3));
%!   ++idx_beam;
%! endfor
%! for i=1:columns(helspr1.X)
%!   idx_node = node_idx_helspr1 - 1 + (1:numel(helspr1.Phi)) + (i - 1) * numel(helspr1.Phi);
%!   elnodes = int32([idx_node(1), node_idx_rb2;
%!                    idx_node(end), node_idx_rb1]);
%!   mesh.elements.beam2(idx_beam + (1:rows(elnodes))) = struct("nodes", mat2cell(elnodes, ones(rows(elnodes), 1, "int32"), 2), ...
%!                                                              "section", mat2cell(repmat(section1, rows(elnodes), 1), ones(rows(elnodes), 1, "int32")), ...
%!                                                              "e2", mat2cell(repmat(helspr1.e2, rows(elnodes), 1), ones(rows(elnodes), 1, "int32"), 3));
%!   idx_beam += rows(elnodes);
%! endfor
%! for i=1:columns(rubspr1.X)
%!   idx_node = int32([1, 2]) + node_idx_rubspr1 - 1 + (i - 1) * 2;
%!   elnodes = int32([idx_node(2), node_idx_rb2]);
%!   mesh.elements.beam2(idx_beam + (1:rows(elnodes))) = struct("nodes", mat2cell(elnodes, ones(rows(elnodes), 1, "int32"), 2), ...
%!                                                              "section", mat2cell(repmat(section2, rows(elnodes), 1), ones(rows(elnodes), 1, "int32")), ...
%!                                                              "e2", mat2cell(repmat(rubspr1.e2, rows(elnodes), 1), ones(rows(elnodes), 1, "int32"), 3));
%!   idx_beam += rows(elnodes);
%! endfor
%! empty_cell = cell(1, 2);
%! mesh.elements.bodies = struct("m", empty_cell, "J", empty_cell, "lcg", empty_cell, "nodes", empty_cell);
%! mesh.elements.bodies(1).m = param.m1;
%! mesh.elements.bodies(1).J = param.J1;
%! mesh.elements.bodies(1).lcg = param.lcg1;
%! mesh.elements.bodies(1).nodes = node_idx_rb1;
%! mesh.elements.bodies(2).m = param.m2;
%! mesh.elements.bodies(2).J = param.J2;
%! mesh.elements.bodies(2).lcg = param.lcg2;
%! mesh.elements.bodies(2).nodes = node_idx_rb2;
%! empty_cell = cell(1, 2);
%! mesh.material_data = struct("E", empty_cell, "nu", empty_cell, "rho", empty_cell, "c", empty_cell);
%! mesh.material_data(1).E = helspr1.material.E;
%! mesh.material_data(1).nu = helspr1.material.nu;
%! mesh.material_data(1).beta = helspr1.material.beta;
%! mesh.material_data(1).alpha = helspr1.material.alpha;
%! mesh.material_data(1).rho = helspr1.material.rho;
%! mesh.material_data(2).E = rubspr1.material.E;
%! mesh.material_data(2).nu = rubspr1.material.nu;
%! mesh.material_data(2).beta = rubspr1.material.beta;
%! mesh.material_data(2).alpha = rubspr1.material.alpha;
%! mesh.material_data(2).rho = rubspr1.material.rho;
%! mesh.material_data(3).E = helspr1.material.E;
%! mesh.material_data(3).nu = helspr1.material.nu;
%! mesh.material_data(3).beta = helspr1.material.beta;
%! mesh.material_data(3).alpha = helspr1.material.alpha;
%! mesh.material_data(3).rho = 0;
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%! for i=1:columns(rubspr1.X)
%!   load_case_dof.locked_dof(node_idx_rubspr1 + 2 * (i - 1), 1:6) = true;
%! endfor
%! load_case_dof.domain = FEM_DO_STRUCTURAL;
%! mesh.materials.beam2 = [repmat(int32(1), columns(helspr1.X) * (numel(helspr1.Phi) - 1), 1);
%!                         repmat(int32(2), columns(rubspr1.X), 1);
%!                         repmat(int32(3), columns(helspr1.X) * 2, 1);
%!                         repmat(int32(3), columns(rubspr1.X), 1)];
%! empty_cell = cell(1, numel(2 * columns(helspr1.X) + columns(rubspr1.X)));
%! mesh.elements.joints = struct("nodes", empty_cell, "C", empty_cell);
%! idx_joint = int32(0);
%! for i=1:columns(helspr1.X)
%!   node_idx_slave = node_idx_helspr1 - 1 + numel(helspr1.Phi) * i;
%!   lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb1, 1:3);
%!   ++idx_joint;
%!   mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb1, node_idx_slave]);
%!   mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                             zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%! endfor
%! for i=1:columns(helspr1.X)
%!   node_idx_slave = node_idx_helspr1 + numel(helspr1.Phi) * (i - 1);
%!   lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!   ++idx_joint;
%!   mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!   mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                             zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%! endfor
%! for i=1:columns(rubspr1.X)
%!   node_idx_slave = node_idx_rubspr1 + 2 * i - 1;
%!   lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!   ++idx_joint;
%!   mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!   mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                             zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%! endfor
%! empty_cell = cell(1, 6);
%! load_case = struct("nodes", empty_cell, "loaded_nodes", empty_cell);
%! for i=1:numel(load_case)
%!   load_case(i).loaded_nodes = int32(node_idx_rb1);
%!   load_case(i).loads = zeros(1, 6);
%!   load_case(i).loads(i) = 1;
%! endfor
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_DAMPING, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! switch (param.damping)
%! case "global"
%!  [damp.alpha_global, damp.beta_global] = fem_pre_mat_rayleigh_damping(damp.D_global, damp.f_global);
%!  mat_ass.D = damp.alpha_global * mat_ass.M + damp.beta_global * mat_ass.K;
%! case "local"
%! otherwise
%!  error("unknown damping model: \"%s\"", param.damping);
%! endswitch
%! opt_linsol.solver = "umfpack";
%! opt_linsol.pre_scaling = true;
%! opt_linsol.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_linsol.refine_max_iter = int32(50);
%! opt_linsol.epsilon_refinement = 1e-10;
%! opt_linsol.verbose = int32(0);
%! opt_eig.p = 5 * param.N;
%! opt_eig.maxit = int32(100);
%! omega = linspace(2 * pi * param.fmin, 2 * pi * param.fmax, param.num_freq_modal);
%! switch(param.solver)
%! case "damped"
%!   [sol_eig, Phi] = fem_sol_modal_damped(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eig);
%!   [Phi, h] = fem_sol_modes_scale(mat_ass.M, mat_ass.K, sol_eig.lambda, Phi, mat_ass.R);
%!   U2 = fem_sol_harmonic_modal(h, sol_eig.lambda, Phi(dof_map.ndof(node_idx_rb2, :), :), omega);
%! case "undamped"
%!   [sol_eig, Phi] = fem_sol_modal(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eig);
%!   [dgen, kgen, rgen] = fem_sol_modes_scale2(mat_ass.M, mat_ass.D, mat_ass.K, Phi, mat_ass.R);
%!   U2 = fem_sol_harmonic_modal2(dgen, kgen, rgen, Phi(dof_map.ndof(node_idx_rb2, :), :), omega);
%! endswitch
%! for i=1:size(sol_eig.def, 3)
%!   sre = max(max(abs(real(sol_eig.def(:, 1:3, i)))));
%!   sim = max(max(abs(imag(sol_eig.def(:, 1:3, i)))));
%!   if (sim > sre)
%!     sol_eig.def(:, :, i) = imag(sol_eig.def(:, :, i)) * param.maxdef / sim;
%!   else
%!     sol_eig.def(:, :, i) = real(sol_eig.def(:, :, i)) * param.maxdef / sre;
%!   endif
%! endfor
%! if (param.num_freq_nodal > 0)
%! omegan = linspace(2 * pi * param.fmin, 2 * pi * param.fmax, param.num_freq_nodal);
%! U2n = zeros(6, numel(omegan), columns(mat_ass.R));
%! for i=1:numel(omegan)
%!  Un = fem_sol_factor(-omegan(i)^2 * mat_ass.M + 1j * omegan(i) * mat_ass.D + mat_ass.K, opt_linsol) \ mat_ass.R;
%!  for j=1:columns(Un)
%!    U2n(:, i, j) = Un(dof_map.ndof(node_idx_rb2, :), j);
%!  endfor
%! endfor
%! endif
%! figure("visible", "off");
%! hold on;
%! semilogy(sol_eig.f * SI_unit_second^-1, 100 * sol_eig.D, "-;D(f);r");
%! xlabel("f [Hz]");
%! ylabel("D [%]");
%! title("modal damping ratio");
%! grid on;
%! grid minor on;
%! for i=1:3
%!   figure("visible", "off");
%!   subplot(2, 1, 1);
%!   hold on;
%!   plot(omega / (2*pi) * (SI_unit_second^-1), 20 * log10(abs(-omega.^2 .* U2(i, :, i)) * (SI_unit_meter / SI_unit_second^2)), "-;modal;r");
%!   if (param.num_freq_nodal > 0)
%!     plot(omegan / (2*pi) * (SI_unit_second^-1), 20 * log10(abs(-omegan.^2 .* U2n(i, :, i)) * (SI_unit_meter / SI_unit_second^2)), "-;nodal;b");
%!   endif
%!   xlabel("f [Hz]");
%!   ylabel("a [dB/(1 m/s^2/N)]");
%!   ylim([-80, 10]);
%!   yticks(-100:5:10);
%!   grid on;
%!   grid minor on;
%!   title("frequency response magnitude");
%!   subplot(2, 1, 2);
%!   hold on;
%!   plot(omega / (2*pi) * (SI_unit_second^-1), 180 / pi * arg(-omega.^2 .* U2(i, :, i)), "-;modal;r");
%!   if (param.num_freq_nodal > 0)
%!     plot(omegan / (2*pi) * (SI_unit_second^-1), 180 / pi * arg(-omegan.^2 .* U2n(i, :, i)), "-;nodal;b");
%!   endif
%!   xlabel("f [Hz]");
%!   ylabel("arg(a) [deg]");
%!   ylim([-180,180]);
%!   yticks(-180:30:180);
%!   grid on;
%!   grid minor on;
%!   title("frequency response phase");
%! endfor
%! opt_post.elem_types = {"beam2"};
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
