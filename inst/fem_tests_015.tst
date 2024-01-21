## fem_tests.m:15
%!test
%! ## TEST 15
%! do_plot = false;
%! close all;
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! scale_def = 25e-3 / SI_unit_m;
%! tol = eps^0.5;
%! tol2 = eps^0.2;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ## 1
%!       0,  0.5 * b,  0.5 * c;  ## 2
%!       0, -0.5 * b,  0.5 * c;  ## 3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ## 4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ## 5
%!       0,  0.5 * b, -0.5 * c;  ## 6
%!       0, -0.5 * b, -0.5 * c;  ## 7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ## 8
%!       a,  0.5 * b,  0.5 * c;  ## 9
%!       a, -0.5 * b,  0.5 * c;  ## 10
%!       a,  0.5 * b, -0.5 * c;  ## 11
%!       a, -0.5 * b, -0.5 * c,  ## 12
%!       a + d,        0,        0;  ## 13
%!       -e,        0,        0]; ## 14
%! algorithms = {"eliminate", "unsymmetric", "shift-invert", "diag-shift-invert"};
%! for ialg=1:numel(algorithms)
%!   data(ialg).mesh.nodes = [X, zeros(rows(X), 3)];
%!   data(ialg).mesh.elements.iso8 = int32([1:8;
%!                                          9, 1, 4, 10, 11, 5, 8, 12]);
%!   data(ialg).mesh.materials.iso8 = int32([1; 1]);
%!   data(ialg).mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%!   data(ialg).mesh.elements.rbe3(1).weight = ones(1, 4);
%!   data(ialg).mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%!   data(ialg).mesh.elements.rbe3(2).weight = ones(1, 4);
%!   E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%!   nu = 0.3;
%!   data(ialg).mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%!   data(ialg).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data(ialg).load_case.locked_dof = false(rows(data(ialg).mesh.nodes), 6);
%!   data(ialg).cms_opt.verbose = false;
%!   data(ialg).cms_opt.modes.number = int32(6);
%!   data(ialg).cms_opt.nodes.modal.number = int32(14);
%!   data(ialg).cms_opt.nodes.interfaces.number = int32(13);
%!   data(ialg).cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!   data(ialg).cms_opt.algorithm = algorithms{ialg};
%!   data(ialg).cms_opt.invariants = false;
%!   [data(ialg).mesh_cms, ...
%!    data(ialg).mat_ass_cms, ...
%!    data(ialg).dof_map_cms, ...
%!    data(ialg).sol_eig_cms] = fem_cms_create(data(ialg).mesh, data(ialg).load_case, data(ialg).cms_opt);
%!   data(ialg).mesh.elements.joints = struct("nodes", cell(1, 1), "C", cell(1,1));
%!   data(ialg).load_case.joints = struct("U", cell(1, 1));
%!   data(ialg).mesh.elements.joints(1).nodes = data(ialg).cms_opt.nodes.modal.number;
%!   data(ialg).mesh.elements.joints(1).C = eye(6);
%!   data(ialg).load_case.joints(1).U = zeros(6, 1);
%!   data(ialg).load_case.loaded_nodes = [data(ialg).cms_opt.nodes.interfaces.number];
%!   data(ialg).load_case.loads = [0, 0, -1, 0, 0, 0] / SI_unit_N;
%!   data(ialg).dof_map = fem_ass_dof_map(data(ialg).mesh, data(ialg).load_case);
%!   [data(ialg).mat_ass.M, ...
%!    data(ialg).mat_ass.K, ...
%!    data(ialg).mat_ass.R, ...
%!    data(ialg).mat_ass.m, ...
%!    data(ialg).mat_ass.mat_info] = fem_ass_matrix(data(ialg).mesh, ...
%!                                                  data(ialg).dof_map, ...
%!                                                  [FEM_MAT_MASS, ...
%!                                                   FEM_MAT_STIFFNESS, ...
%!                                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                                   FEM_SCA_TOT_MASS], ...
%!                                                  data(ialg).load_case);
%!   mref = data(ialg).mesh.material_data.rho * a * b * c;
%!   assert_simple(data(ialg).mat_ass.m, mref, sqrt(eps) * mref);
%!   [Mred, Rred] = fem_ass_matrix(data(ialg).mesh_cms, ...
%!                                 data(ialg).dof_map_cms, ...
%!                                 [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                  rmfield(data(ialg).load_case, "joints"));
%!   U = full(data(ialg).mat_ass.K \ data(ialg).mat_ass.R);
%!   def_stat = fem_post_def_nodal(data(ialg).mesh, data(ialg).dof_map, U);
%!   ured = data(ialg).mat_ass_cms.Tred * (data(ialg).mat_ass_cms.Kred \ (data(ialg).mat_ass_cms.Tred.' * Rred(data(ialg).dof_map_cms.idx_node, :)));
%!   Ured = zeros(data(ialg).dof_map_cms.totdof, columns(Rred));
%!   Ured(data(ialg).dof_map_cms.idx_node, :) = ured;
%!   def_red_stat = fem_post_def_nodal(data(ialg).mesh, data(ialg).dof_map_cms, Ured);
%!   [data(ialg).mat_ass.Tc, data(ialg).mat_ass.Kc, data(ialg).mat_ass.Mc] = fem_cms_constr_elim(data(ialg).mesh, data(ialg).dof_map, data(ialg).mat_ass);
%!   [phi, lambda] = fem_sol_eigs(data(ialg).mat_ass.Kc, data(ialg).mat_ass.Mc, data(ialg).cms_opt.modes.number);
%!   Phi = zeros(data(ialg).dof_map.totdof, columns(phi));
%!   Phi(data(ialg).dof_map.idx_node, :) = data(ialg).mat_ass.Tc * phi;
%!   def_modal = fem_post_def_nodal(data(ialg).mesh, data(ialg).dof_map, Phi);
%!   [phi_red, lambda_red] = fem_sol_eigs(data(ialg).mat_ass_cms.Kred, data(ialg).mat_ass_cms.Mred, data(ialg).cms_opt.modes.number);
%!   Phi_red = zeros(data(ialg).dof_map_cms.totdof, columns(phi_red));
%!   Phi_red(data(ialg).dof_map_cms.idx_node, :) = data(ialg).mat_ass_cms.Tred * phi_red;
%!   def_red_modal = fem_post_def_nodal(data(ialg).mesh, data(ialg).dof_map_cms, Phi_red);
%!   for i=1:size(def_red_modal, 3)
%!     phi = reshape(def_modal(:, :, i).', size(def_modal, 1) * size(def_modal, 2), 1);
%!     phi /= max(abs(phi));
%!     min_fPhi = inf;
%!     min_flambda = inf;
%!     for j=1:columns(Phi_red)
%!       phi_red = reshape(def_red_modal(:, :, j).', size(def_red_modal, 1) * size(def_red_modal, 2), 1);
%!       phi_red /= max(abs(phi_red));
%!       if norm(phi + phi_red) < norm(phi - phi_red)
%!         phi_red *= -1;
%!       endif
%!       fPhi = norm(phi_red - phi);
%!       flambda = abs(lambda_red(j) - lambda(i));
%!       if fPhi < min_fPhi && flambda < min_flambda
%!         min_fPhi = fPhi;
%!         min_flambda = flambda;
%!         phi_red_opt = phi_red;
%!         lambda_red_opt = lambda_red(j);
%!       endif
%!     endfor
%!     assert_simple(phi_red_opt, phi, tol2 * norm(phi));
%!     assert_simple(lambda_red_opt, lambda(i), tol2 * max(abs(lambda)));
%!   endfor
%!   assert_simple(lambda_red, lambda, tol2 * max(abs(lambda)));
%!   for i=1:size(def_stat, 3)
%!     assert_simple(def_red_stat, def_stat, tol * max(norm(def_stat(:, :, i), "rows")));
%!   endfor
%!   assert_simple(isdefinite(data(ialg).mat_ass_cms.Kred), true);
%!   assert_simple(isdefinite(data(ialg).mat_ass_cms.Mred), true);
%!   for i=1:3
%!     assert_simple(sum(data(ialg).mat_ass_cms.diagM(i:3:end)), data(ialg).mat_ass.m, tol * data(ialg).mat_ass.m);
%!   endfor
%!   if (do_plot)
%!     for i=1:numel(data(ialg).sol_eig_cms.f)
%!       figure("visible","off");
%!       fem_post_sol_plot(data(ialg).mesh, data(ialg).sol_eig_cms, scale_def / max(norm(data(ialg).sol_eig_cms.def(:, 1:3, i), "rows")), i);
%!       xlabel("x [m]");
%!       ylabel("y [m]");
%!       zlabel("z [m]");
%!       grid on;
%!       grid minor on;
%!       title(sprintf("mode %d: f=%.0fHz", i, data(ialg).sol_eig_cms.f(i)));
%!     endfor
%!     figure_list();
%!   endif
%! endfor
%! for i=2:numel(data)
%!   assert_simple([data(i).sol_eig_cms.f],[data(1).sol_eig_cms.f], tol * max([data(1).sol_eig_cms.f]));
%!   assert_simple(data(i).mat_ass.M, data(1).mat_ass.M, 0);
%!   assert_simple(data(i).mat_ass.K, data(1).mat_ass.K, 0);
%! endfor
