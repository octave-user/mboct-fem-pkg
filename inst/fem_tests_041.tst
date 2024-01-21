## fem_tests.m:41
%!test
%! ## TEST 41
%! ## Robert Gasch, Klaus Knothe 1989
%! ## Strukturdynamik Band 2
%! ## Kontinua und ihre Diskretisierung
%! ## Equation 2.38b, figure 2.11, page 26, chapter 2
%! do_plot = false;
%! L = 2000e-3;
%! w = 10e-3;
%! h = 20e-3;
%! c2 = 0.291;
%! E = 70000e6;
%! nu = 0.3;
%! rho = 2700;
%! N = 100;
%! A = w * h;
%! Ay = 5 / 6 * w * h;
%! Az = 5 / 6 * w * h;
%! It = c2 * h * w^3;;
%! Iy = w * h^3 / 12;
%! Iz = w^3 * h / 12;
%! B = E * Iz;
%! mu = rho * A;
%! Delta = linspace(1e-6, 10, 1000) / L;
%! omega = sqrt(Delta.^4 * B / mu);
%! P0 = 1;
%! wstat_a = P0 * L^3 / (3 * B);
%! V_a = 3 ./ (Delta * L).^3 .* (sin(Delta * L) .* cosh(Delta * L) - cos(Delta * L) .* sinh(Delta * L)) ./ (1 + cos(Delta * L) .* cosh(Delta * L));
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
%! mesh.elements.beam2 = repmat(beam1, 1, N - 1);
%! for i=1:N - 1
%!   mesh.elements.beam2(i).nodes = int32(i:i+1);
%! endfor
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.locked_dof(1, :) = true;
%! load_case.loaded_nodes = int32(N);
%! load_case.loads = [0, P0, 0, 0, 0, 0];
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_MAT_MASS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! wdyn = zeros(1, numel(omega));
%! Ustat = mat_ass.K \ mat_ass.R;
%! wstat = Ustat(dof_map.ndof(N, 2));
%! for i=1:numel(omega)
%!   U = (-omega(i)^2 * mat_ass.M + mat_ass.K) \ mat_ass.R;
%!   wdyn(i) = U(dof_map.ndof(N, 2));
%! endfor
%! if (do_plot)
%!   figure("visible", "off");
%!   hold on;
%!   plot(Delta * L, wdyn / wstat_a, "-;V(omega);1");
%!   plot(Delta * L, V_a, "-;V_r_e_f(omega);0");
%!   xlabel("Delta * L [1]");
%!   ylabel("wdyn/wstat [1]");
%!   ylim([-1, 2]);
%!   grid on;
%!   grid minor on;
%!   title("frequency response cantilever beam y-direction");
%! endif
%! assert_simple(wstat, wstat_a, 1e-3 * wstat_a);
%! idx = find(abs(V_a) < 2);
%! assert_simple(mean(abs(wdyn(idx) / wstat_a - V_a(idx))) < 1e-2);
