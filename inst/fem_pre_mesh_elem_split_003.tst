## fem_pre_mesh_elem_split.m:03
%!test
%! do_plot = false;
%! if (do_plot)
%! close all;
%! endif
%! a = 20e-3;
%! b = 15e-3;
%! c = 10e-3;
%! rho = 7850;
%! m = 3 / 8 * a * b * c * rho;
%! tol_m = eps^0.9;
%! N = 10;
%! r = 1;
%! X = [0.5 * a, 0.5 * b, c;
%!            0, 0.5 * b, c;
%!            0,       0, c;
%!            a,       0, c;
%!      0.5 * a, 0.5 * b, 0;
%!            0, 0.5 * b, 0;
%!            0,       0, 0;
%!            a,       0, 0];
%!
%! mesh_data(1).mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh_data(1).mesh.elements.iso8 = int32([1:8]);
%! mesh_data(1).mesh.materials.iso8 = int32([1]);
%! mesh_data(1).load_case.locked_dof = false(rows(mesh_data(1).mesh.nodes), 6);
%! E = 210000e6;
%! nu = 0.3;
%! mesh_data(1).mesh.material_data.rho = rho;
%! mesh_data(1).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! options.iso8.max_angle = 95 * pi / 180;
%! mesh_data(2).mesh = fem_pre_mesh_elem_split(mesh_data(1).mesh, options);
%! mesh_data(2).load_case = mesh_data(1).load_case;
%! for i=1:numel(mesh_data)
%!   mesh_data(i).dof_map = fem_ass_dof_map(mesh_data(i).mesh, mesh_data(i).load_case);
%!   [mesh_data(i).mat_ass.M, ...
%!    mesh_data(i).mat_ass.K, ...
%!    mesh_data(i).mat_ass.m] = fem_ass_matrix(mesh_data(i).mesh, ...
%!                                             mesh_data(i).dof_map, ...
%!                                             [FEM_MAT_MASS, ...
%!                                              FEM_MAT_STIFFNESS, ...
%!                                              FEM_SCA_TOT_MASS], ...
%!                                             mesh_data(i).load_case);
%!   mesh_data(i).sol_eig = fem_sol_modal(mesh_data(i).mesh, mesh_data(i).dof_map, mesh_data(i).mat_ass, N, r);
%! endfor
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   for i=1:numel(mesh_data)
%!     for j=8:9
%!       opts.scale_def = 2.5e-3 / max(max(abs(mesh_data(i).sol_eig.def(:, 1:3, j))));
%!       opts.output_step_idx = j;
%!       if (do_plot)
%!         fem_post_sol_external(mesh_data(i).mesh, mesh_data(i).sol_eig, opts);
%!         [img, map, alpha] = imread(sprintf("%s_%03d.jpg", opts.print_to_file, 1));
%!         figure("visible", "off");
%!         imshow(img, map);
%!         title(sprintf("mode %d f=%.0fHz", j, mesh_data(i).sol_eig.f(j)));
%!       endif
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     fn = dir([opts.print_to_file, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! for i=1:numel(mesh_data)
%!   assert_simple(mesh_data(i).mat_ass.m, m, tol_m * m);
%! endfor
%! if (do_plot)
%!   figure_list();
%! endif
