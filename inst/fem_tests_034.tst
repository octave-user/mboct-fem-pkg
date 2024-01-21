## fem_tests.m:34
%!test
%! ## TEST 34
%! ## Cantilever beam with rectangular cross section and lateral load
%! ## W.Beitz, K.-H.Grothe, 1997, Dubbel, section 2.4.6, page C17, figure 23
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! Fz = -15000;
%! h = 0.5e-3;
%! geometry.l = 150e-3;
%! geometry.w = h;
%! geometry.h = 15e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [ 0; 0; Fz ];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! mesh.nodes(:, 3) -= 0.5 * geometry.h;
%! [dof_map] = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! x = linspace(0, geometry.l, 100);
%! z = linspace(-0.5 * geometry.h, 0.5 * geometry.h, 50);
%! [xx, zz] = meshgrid(x, z);
%! xtauel = mesh.nodes(:, 1)(mesh.elements.iso8);
%! ytauel = mesh.nodes(:, 1)(mesh.elements.iso8);
%! ztauel = mesh.nodes(:, 3)(mesh.elements.iso8);
%! tauxxel = sol_stat.stress.taum.iso8(:, :, 1);
%! tauxzel = sol_stat.stress.taum.iso8(:, :, 6);
%! tauxx = griddata(xtauel(:), ztauel(:), tauxxel(:), xx, zz);
%! tauxz = griddata(xtauel(:), ztauel(:), tauxzel(:), xx, zz);
%! Iy = geometry.w * geometry.h^3 / 12;
%! tauxx_a = -Fz / Iy * (geometry.l - xx) .* zz;
%! tauxz_a = 3 / 2 * Fz / (geometry.h * geometry.w) * (1 - (zz / (0.5 * geometry.h)).^2);
%! idx_x = find((xx(:) > 0.1 * geometry.l) & (xx(:) < 0.9 * geometry.l));
%! assert_simple(tauxx(:)(idx_x), tauxx_a(:)(idx_x), 1e-2 * max(tauxx_a(:)(idx_x)));
%! assert_simple(tauxz(:)(idx_x), tauxz_a(:)(idx_x), 7e-2 * max(abs(tauxz_a(:)(idx_x))));
