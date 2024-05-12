## fem_cms_export.m:02
%!test
%! try
%! close all;
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
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.element.name = "elem_id_cube";
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.modal.name = "node_id_cube_modal";
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.nodes.interfaces.name = "node_id_cube_interface1";
%! cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%! cms_opt.algorithm = "shift-invert";
%! cms_opt.invariants = true;
%! [mesh_cms, ...
%!  mat_ass_cms, ...
%!  dof_map_cms, ...
%!  sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
