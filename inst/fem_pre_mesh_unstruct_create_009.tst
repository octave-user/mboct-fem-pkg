## fem_pre_mesh_unstruct_create.m:09
%!test
%! try
%! ## TEST9
%! close all;
%! SI_unit_meter = 1;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! geo.D = 148e-3 / SI_unit_meter;
%! geo.d = 140e-3 / SI_unit_meter;
%! geo.h = 10e-3 / SI_unit_meter;
%! param.rho = 1.33 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.c = 225 / (SI_unit_meter / SI_unit_second);
%! param.eta = 8.3492e-06 / (SI_unit_pascal * SI_unit_second);
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     geo_file = [filename, ".geo"];
%!     fd = fopen(geo_file, "w");
%!     if (fd == -1)
%!       error("failed to open file %s", geo_file);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "vout = newv;\n");
%!     fputs(fd, "Sphere(vout) = {0, 0, 0, 0.5 * D};\n");
%!     fputs(fd, "vin = newv;\n");
%!     fputs(fd, "Sphere(vin) = {0, 0, 0, 0.5 * d};\n");
%!     fputs(fd, "v = newv;\n");
%!     fputs(fd, "BooleanDifference(v) = {Volume{vout}; Delete; }{ Volume{vin}; Delete; };\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {v};\n");
%!     fputs(fd, "Physical Surface(\"surface\", 8) = {2};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v}; } } = h;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   opt.mesh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh = fem_pre_mesh_reorder(mesh);
%!   mesh.material_data.c = 220;
%!   mesh.material_data.rho = 1.33;
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.c = param.c;
%!   mesh.material_data.eta = param.eta;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.coll_Ka, ...
%!    mat_ass.coll_Ma] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                       FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                       FEM_VEC_COLL_STIFF_ACOUSTICS, ...
%!                                       FEM_VEC_COLL_MASS_ACOUSTICS], ...
%!                                      load_case);
%!   N = 80;
%!   fref = 2 * pi * 1 / (1 / SI_unit_second);
%!   opt_sol.rho = fref^2;
%!   opt_sol.tolerance = sqrt(eps);
%!   opt_sol.algorithm = "shift-invert";
%!   opt_sol.solver = "umfpack";
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.refine_max_iter = int32(10);
%!   opt_sol.pre_scaling = true;
%!   opt_so.problem = "acoustic";
%!   [Phi, lambda, err] = fem_sol_eigs(mat_ass.Ka, mat_ass.Ma, N, opt_sol);
%!   sol.p = -Phi(dof_map.ndof, :) * diag(imag(lambda));
%!   sol.f = imag(lambda) / (2 * pi);
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
