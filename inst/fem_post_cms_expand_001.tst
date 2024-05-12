## fem_post_cms_expand.m:01
%!test
%! try
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!             0,  0.5 * b,  0.5 * c;  ##  2
%!             0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!             0,  0.5 * b, -0.5 * c;  ##  6
%!             0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!             a,  0.5 * b,  0.5 * c;  ##  9
%!             a, -0.5 * b,  0.5 * c;  ## 10
%!             a,  0.5 * b, -0.5 * c;  ## 11
%!             a, -0.5 * b, -0.5 * c,  ## 12
%!         a + d,        0,        0;  ## 13
%!            -e,        0,        0]; ## 14
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! cms_data.load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%! cms_opt.algorithm = "eliminate";
%! cms_opt.invariants = true;
%! [cms_data.mesh, ...
%!  cms_data.mat_ass, ...
%!  cms_data.dof_map, ...
%!  cms_data.sol_eig, ...
%!  cms_data.cms_opt] = fem_cms_create(mesh, cms_data.load_case, cms_opt);
%! sol_dyn.t = linspace(0, 1e-3, 100);
%! sol_dyn.bodies.q = zeros(columns(cms_data.mat_ass.Tred), numel(sol_dyn.t));
%! f = 1000 + rand(1, columns(cms_data.mat_ass.Tred)) * 1000;
%! for i=1:rows(sol_dyn.bodies.q)
%!   sol_dyn.bodies.q(i, :) = sin(2 * pi * f(i) * sol_dyn.t);
%! endfor
%! sol_dyn.bodies.X = zeros(3, numel(sol_dyn.t));
%! sol_dyn.bodies.R = zeros(3, 3, numel(sol_dyn.t));
%! for i=1:numel(sol_dyn.t)
%!   sol_dyn.bodies.R(:, :, i) = eye(3);
%! endfor
%! sol_dyn.bodies.X_ref = zeros(3, numel(sol_dyn.t));
%! sol_dyn.bodies.R_ref = sol_dyn.bodies.R;
%! scale_types = {"modal node", "modal node*", "least square", "least square*", "reference node", "reference node*"};
%! sol_tot = struct("bodies", cell(numel(scale_types), 2));
%! for j=1:2
%!   for i=1:numel(scale_types)
%!     opt_exp.scale = 1;
%!     opt_exp.scale_type = scale_types{i};
%!     opt_exp.output_stress = FEM_SCA_STRESS_VMIS;
%!     sol_tot(i, j) = fem_post_cms_expand(sol_dyn, cms_data, 1:10:numel(sol_dyn.t), opt_exp);
%!   endfor
%!   if (j == 1)
%!     cms_data.cms_opt.selected_modes = 1:2:columns(cms_data.mat_ass.Tred);
%!     sol_dyn.bodies.q = sol_dyn.bodies.q(cms_data.cms_opt.selected_modes, :);
%!   endif
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
