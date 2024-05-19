## fem_cms_create.m:02
%!test
%! try
%! close all;
%! f_run_post_proc = false;
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
%!       0,  0.5 * b,  0.5 * c;  ##  2
%!       0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!       0,  0.5 * b, -0.5 * c;  ##  6
%!       0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!       a,  0.5 * b,  0.5 * c;  ##  9
%!       a, -0.5 * b,  0.5 * c;  ## 10
%!       a,  0.5 * b, -0.5 * c;  ## 11
%!       a, -0.5 * b, -0.5 * c,  ## 12
%!       a + d,        0,        0;  ## 13
%!       -e,        0,        0]; ## 14
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
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.verbose = false;
%! cms_opt.modes.number = 5;
%! cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%! cms_opt.algorithm = "shift-invert";
%! cms_opt.invariants = true;
%! cms_opt.refine_max_iter = int32(10);
%! cms_opt.solver = "pastix";
%! [mesh_cms, ...
%!  mat_ass_cms, ...
%!  dof_map_cms, ...
%!  sol_eig_cms, ...
%!  cms_opt] = fem_cms_create(mesh, load_case, cms_opt);
%! if (f_run_post_proc)
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   opt_post.scale_def = 0.2;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_eig_cms, opt_post);
%!   fn = dir([opt_post.print_to_file, "*.jpg"]);
%!   for i=1:numel(fn)
%!     figure("visible", "off");
%!     [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!     imshow(img, alpha);
%!     title(sprintf("mode %d f=%.0fHz", i, sol_eig_cms.f(i) / SI_unit_s));
%!   endfor
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
