## fem_tests.m:32
%!test
%! try
%! ## TEST 32
%! close all;
%! scale_stat = 1;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! material.alpha = 1e-10;
%! material.beta = 1e-8;
%! h = 2.5e-3;
%! geometry.l = 1000e-3;
%! geometry.w = 2.5e-3;
%! geometry.h = 25e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! do_plot = false;
%! options = struct();
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, zeros(3,1));
%! load_case.g = [0; 0; -9.81];
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_MAT_STIFFNESS, ...
%!                                           FEM_VEC_LOAD_CONSISTENT], ...
%!                                          load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! idx_itf = find(mesh.nodes(:, 1) == geometry.l);
%! mesh_cms = mesh;
%! cms_opt.nodes.interfaces.number = rows(mesh_cms.nodes) + 1;
%! cms_opt.nodes.modal.number = rows(mesh_cms.nodes) + 2;
%! mesh_cms.nodes(cms_opt.nodes.interfaces.number, 1:3) = [geometry.l + 0.05, 0.5 * geometry.w, 0.5 * geometry.h];
%! mesh_cms.nodes(cms_opt.nodes.modal.number, 1:3) = [-0.05, 0.5 * geometry.w, 0.5 * geometry.h];
%! mesh_cms.elements.rbe3(1).nodes = [cms_opt.nodes.interfaces.number, idx_itf.'];
%! mesh_cms.elements.rbe3(1).weight = ones(1, numel(idx_itf));
%! load_case_cms.locked_dof = false(rows(mesh_cms.nodes), 6);
%! load_case_cms.locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%! load_case_cms.locked_dof(find(mesh_cms.nodes(:, 1) == 0), 1:3) = true;
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(30);
%! cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%! cms_opt.algorithm = "shift-invert";
%! cms_opt.scaling = "mean K,M";
%! cms_opt.invariants = true;
%! cms_opt.refine_max_iter = int32(250);
%! cms_opt.epsilon_refinement = eps^0.8;
%! [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms] = fem_cms_create(mesh_cms, load_case_cms, cms_opt);
%! Dred = mat_ass_cms.Mred * material.alpha + mat_ass_cms.Kred * material.beta;
%! assert_simple(mat_ass_cms.Dred, Dred, sqrt(eps) * max(max(abs(Dred))));
%! mat_ass_cms.Rred = mat_ass_cms.Inv3.' * load_case.g;
%! Ured = mat_ass_cms.Kred \ mat_ass_cms.Rred;
%! sol_cms.def = fem_post_def_nodal(mesh_cms, dof_map_cms, mat_ass_cms.Tred * Ured);
%! w = geometry.w;
%! h = geometry.h;
%! l = geometry.l;
%! rho = material.rho;
%! A = w * h;
%! Iy = w * h^3 / 12;
%! qz = rho * A * load_case.g(3);
%! z = l - mesh.nodes(:, 1);
%! wz = qz * l^4 / (24 * material.E * Iy) * (3 - 4 * z / l + (z / l).^4);
%! tol = 1e-2;
%! if (do_plot)
%!   figure("visible","off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat, scale_stat);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title('deformed mesh');
%!   figure_list();
%! endif
%! assert_simple(sol_stat.def(:, 3), wz, tol * max(abs(wz)));
%! assert_simple(sol_cms.def(1:rows(mesh.nodes), 3), wz, tol * max(abs(wz)));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
