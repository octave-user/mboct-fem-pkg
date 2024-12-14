## fem_pre_mesh_struct_create.m:11
%!test
%! try
%! ## TEST11
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.number_of_modes = 10;
%! param.Udyn = eye(3) / SI_unit_meter;
%! param.fmin = 0 / (SI_unit_second^-1);
%! param.fmax = 15000 / (SI_unit_second^-1);
%! param.num_freq = 1000;
%! geometry(1).user_data.helspr.L = 25.8e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.d = 1.3e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.n = 5.3;
%! geometry(1).user_data.helspr.ni = 3;
%! geometry(1).user_data.helspr.ng = 0.75;
%! geometry(1).user_data.color = "y";
%! geometry(2).user_data.helspr.L = 27.7e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.d = 1.3e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.n = 5.7;
%! geometry(2).user_data.helspr.ni = 2.7;
%! geometry(2).user_data.helspr.ng = 0.75;
%! geometry(2).user_data.color = "g";
%! geometry(3).user_data.helspr.L = 28.63e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.d = 1.25e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.n = 6;
%! geometry(3).user_data.helspr.ni = 2.7;
%! geometry(3).user_data.helspr.ng = 0.75;
%! geometry(3).user_data.color = "b";
%! geometry = geometry(1);
%! material.E = 206000e6 / SI_unit_pascal;
%! material.G = 81500e6 / SI_unit_pascal;
%! material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! damp.D = [1e-2; 0.03e-2];
%! damp.f = [20, 15000] / (SI_unit_second^-1);
%! param.Fz = -7.2 / SI_unit_newton;
%! [material.alpha, material.beta] = fem_pre_mat_rayleigh_damping(damp.D, damp.f);
%! param.omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%!
%! function [x, y, z, R, Phi] = helspr_geo(geo, r, s, t, varargin)
%!   Phi = 2 * pi * (r * geo.n + geo.ni);
%!   r1 = 0.5 * geo.d * s;
%!   Theta = 2 * pi * t;
%!   x1 = [0.5 * geo.D * cos(Phi);
%!         0.5 * geo.D * sin(Phi);
%!         (geo.L - geo.d * (2 * (geo.ni - geo.ng) + 1)) * r + geo.d * (geo.ni - geo.ng + 0.5)];
%!   x2 = [0;
%!         0;
%!         x1(3)];
%!   e2 = x2 - x1;
%!   e1 = [-0.5 * geo.D * sin(Phi) * 2 * pi * geo.n;
%!          0.5 * geo.D * cos(Phi) * 2 * pi * geo.n;
%!          (geo.L - geo.d * (2 * (geo.ni - geo.ng) + 1))];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R1 = [e1, e2, e3];
%!   R1 *= diag(1 ./ norm(R1, "cols"));
%!   assert_simple(R1.' * R1, eye(3), eps^0.9);
%!   assert_simple(R1 * R1.', eye(3), eps^0.9);
%!   x3 = R1 * [0; r1 * cos(Theta); r1 * sin(Theta)] + x1;
%!   x = x3(1);
%!   y = x3(2);
%!   z = x3(3);
%! endfunction
%!
%! function [F, locked] = boundary_cond(r, s, t, geo, loads, varargin)
%!   locked = [];
%!   F = [];
%! endfunction
%!
%! function p = pressure_boundary_cond(r, s, t, geometry, load, perm_idx, varargin)
%!   p = [];
%! endfunction
%!
%! function [Freact, sol_stat, sol_eig, sol_eig_def] = transmissibility(geometry, material, param)
%!   material.nu = material.E / (2 * material.G) - 1;
%!   geometry.user_data.helspr.D = geometry.user_data.helspr.Di + geometry.user_data.helspr.d;
%!   kz = material.G * geometry.user_data.helspr.d^4 / (8 * geometry.user_data.helspr.n * geometry.user_data.helspr.D^3);
%!   Ustat = [0; 0; param.Fz / kz];
%!   h = geometry.user_data.helspr.d * pi / 4;
%!   geometry.user_data.helspr.nPhi = max([2, round(sqrt((geometry.user_data.helspr.D * pi * geometry.user_data.helspr.n)^2 + geometry.user_data.helspr.L^2) / h)]) + 1;
%!   geometry.user_data.helspr.nr = max([1, round(0.5 * geometry.user_data.helspr.d / h)]) + 1;
%!   geometry.user_data.helspr.nTheta = max([3, round(geometry.user_data.helspr.d * pi / h)]) + 1;
%!   geometry.mesh_size.r = linspace(0, 1, geometry.user_data.helspr.nPhi);
%!   geometry.mesh_size.s = linspace(0, 1, geometry.user_data.helspr.nr);
%!   geometry.mesh_size.t = linspace(0, 1, geometry.user_data.helspr.nTheta);
%!   geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.helspr.D;
%!   geometry.spatial_coordinates = @(r, s, t, geo, varargin) feval("helspr_geo", geometry.user_data.helspr, r, s, t, varargin);
%!   geometry.material_selector = @(r, s, t, geo, varargin) int32(1);
%!   geometry.boundary_condition = @(r, s, t, geo, load, varargin) feval("boundary_cond", r, s, t, geo, load, varargin);
%!   geometry.pressure_boundary_condition = @(r, s, t, geo, load, perm_idx, varargin) feval("pressure_boundary_cond", r, s, t, geo, load, perm_idx, varargin);
%!   options.elem_type = "iso20";
%!   loads = struct();
%!   [mesh, load_case_dof] = fem_pre_mesh_struct_create(geometry, loads, material, options);
%!   [mesh] = fem_pre_mesh_solid_to_surf(mesh);
%!   idx_node_bottom = unique(mesh.structured.inode_idx(1, :, :)(:));
%!   idx_node_top = unique(mesh.structured.inode_idx(end, :, :)(:));
%!   idx_node_bottom = idx_node_bottom(idx_node_bottom > 0);
%!   idx_node_top = idx_node_top(idx_node_top > 0);
%!   idx_node_joint = [idx_node_bottom; idx_node_top];
%!   empty_cell = cell(1, numel(idx_node_joint));
%!   mesh.elements.joints = struct("nodes", mat2cell(idx_node_joint, ones(numel(idx_node_joint), 1, "int32"), 1), "C", repmat({[eye(3), zeros(3, 3)]}, numel(idx_node_joint), 1));
%!   [dof_map] = fem_ass_dof_map(mesh, load_case_dof);
%!   load_case_stat = fem_pre_load_case_create_empty(1);
%!   for i=1:numel(load_case_stat)
%!     load_case_stat(i).joints = struct("U", repmat({zeros(3, 1)}, numel(idx_node_joint), 1));
%!     for j=1:numel(idx_node_top)
%!       load_case_stat(i).joints(numel(idx_node_bottom) + j).U = Ustat(:, i);
%!     endfor
%!   endfor
%!   load_case_dyn = fem_pre_load_case_create_empty(columns(param.Udyn));
%!   for i=1:numel(load_case_dyn)
%!     load_case_dyn(i).joints = struct("U", repmat({zeros(3, 1)}, numel(idx_node_joint), 1));
%!     for j=1:numel(idx_node_top)
%!       load_case_dyn(i).joints(numel(idx_node_bottom) + j).U = param.Udyn(:, i);
%!     endfor
%!   endfor
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS, ...
%!                                         FEM_MAT_STIFFNESS, ...
%!                                         FEM_VEC_LOAD_CONSISTENT], ...
%!                                        load_case_stat);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, param.number_of_modes);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case_stat, ...
%!                                    sol_stat);
%!   for i=1:numel(load_case_dyn)
%!     load_case_dyn(i).tau0 = sol_stat.stress.tau;
%!   endfor
%!   mesh_def = mesh;
%!   mesh_def.nodes += sol_stat.def;
%!   [mat_ass_def.M, ...
%!    mat_ass_def.D, ...
%!    mat_ass_def.K, ...
%!    mat_ass_def.KTAU0, ...
%!    mat_ass_def.R, ...
%!    mat_ass_def.mat_info, ...
%!    mat_ass_def.mesh_info] = fem_ass_matrix(mesh_def, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_TAU0, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_dyn);
%!   mat_ass_def.K += mat_ass_def.KTAU0;
%!   [sol_eig_def] = fem_sol_modal(mesh_def, dof_map, mat_ass_def, param.number_of_modes);
%!   Freact = complex(zeros(3, columns(mat_ass_def.R), numel(param.omega)));
%!   opt_factor.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_factor.verbose = int32(0);
%!   for i=1:numel(param.omega)
%!     fprintf(stderr, "%3d:%.2fHz\n", i, param.omega(i) / (2 * pi));
%!     A = -param.omega(i)^2 * mat_ass_def.M + 1j * param.omega(i) * mat_ass_def.D + mat_ass_def.K;
%!     Uij = fem_sol_factor(A, opt_factor) \ mat_ass_def.R;
%!     for j=1:size(Freact, 2)
%!       Freact(:, j, i) = sum(Uij(:, j)(dof_map.edof.joints(1:numel(idx_node_bottom), :)), 1)(:) * mat_ass_def.mat_info.beta(3);
%!     endfor
%!   endfor
%! endfunction
%! Freact = cell(1, numel(geometry));
%! for i=1:numel(geometry)
%!   [Freact{i}] = transmissibility(geometry(i), material, param);
%! endfor
%! figure("visible", "off");
%! hold on;
%! for i=1:numel(geometry)
%!   plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * sqrt(Freact{i}(1, 1, :).^2 + Freact{i}(2, 2, :).^2 + Freact{i}(3, 3, :).^2)), sprintf("-;%d;%s", i, geometry(i).user_data.color));
%! endfor
%! xlabel("f [Hz]");
%! ylabel("kdyn [dB/(1N/m)]");
%! title("overall transmissibility");
%! grid minor on;
%! for j=1:columns(param.Udyn)
%!   figure("visible", "off");
%!   hold on;
%!   for i=1:numel(geometry)
%!     plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * abs(Freact{i}(i, i, :))), sprintf("-;%d;%s", i, geometry(i).user_data.color));
%!   endfor
%!   xlabel("f [Hz]");
%!   ylabel(sprintf("kdyn%s [dB/(1N/m)]", {"x", "y", "z"}{j}));
%!   title(sprintf("transmissibility %s-direction", {"x", "y", "z"}{j}));
%!   grid minor on;
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
