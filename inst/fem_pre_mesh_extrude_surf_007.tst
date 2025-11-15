%!test
%! try
%! ### TEST6
%! ##########################################################################
%! ## VIBRATIONS OF COMPLETE SPHERICAL SHELLS WITH IMPERFECTIONS
%! ## Thomas A. Duffey
%! ## Jason E. Pepin
%! ## Amy N. Robertson
%! ## Michael L. Steinzig
%! ## Internatial Modal Analysis Conference (IMAC-XXIII)
%! ## Orlando, Florida
%! ## January 31-February 3, 2005
%! ## Los Alamos
%! ## NATIONAL LABORATORY
%! ##########################################################################
%! R = 4.4688 * 25.4e-3;                  # radius in the middle of the shell
%! h = 2 * R * pi / 72;                   # mesh size
%! t = 0.0625 * 25.4e-3;                  # shell thickness
%! E = 28e6 * 6895;                       # Young's modulus
%! nu = 0.28;                             # Poisson ratio
%! rho = 0.000751 * 4.4482 / (25.4e-3^4); # density
%! tol_coherence = 1e-4;                  # relative tolerance
%! num_modes = 39;                        # number of modes to compute
%! fref = [5078, 6005, 6378, 6729];       # reference solution
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R=%g;\n", R - 0.5 * t);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0, h};\n");
%!     fputs(fd, "Point(2) = {0.0,0.0,-R, h};\n");
%!     fputs(fd, "Point(3) = {R,0.0,0.0, h};\n");
%!     fputs(fd, "Point(4) = {0.0,0.0,R, h};\n");
%!     fputs(fd, "Circle(1) = {2,1,3};\n");
%!     fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!     fputs(fd, "Transfinite Curve{1} = Round(R * Pi / 2. / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve{2} = Round(R * Pi / 2. / h) + 1;\n");
%!     fputs(fd, "tmp1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi}{\n");
%!     fputs(fd, "  Curve{1,2}; Layers{Round(Pi * R / h)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "tmp2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi}{\n");
%!     fputs(fd, "  Curve{tmp1[0],tmp1[4]}; Layers{Round(Pi * R / h)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Surface(\"surface\", 1) = {2, 4, 1, 3};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   ##spawn_wait(spawn("gmsh", {[filename, ".geo"]}));return;
%!   pid = spawn("gmsh", {"-format", "msh2", "-2", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"tria6h", "quad8"};
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh);
%!   [mesh.nodes, mesh.elements.iso20] = fem_pre_mesh_extrude_surf(mesh, "iso20", 1, t);
%!   [mesh.nodes, mesh.elements.penta15] = fem_pre_mesh_extrude_surf(mesh, "penta15", 1, t);
%!   [mesh, dx] = fem_pre_mesh_coherence(mesh, tol_coherence * R);
%!   mesh = fem_pre_mesh_reorder(mesh);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   load_case = struct();
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS, ...
%!                                         FEM_MAT_STIFFNESS], ...
%!                                        load_case);
%!   opt_eig.rho = sqrt(eps) * max(abs(diag(mat_ass.K))) / max(abs(diag(mat_ass.M)));
%!   opt_eig.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_eig.solver = "pardiso";
%!   opt_eig.algorithm = "shift-invert";
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, opt_eig);
%!   sol_eig.stress = fem_ass_matrix(mesh, dof_map, [FEM_SCA_STRESS_VMIS], load_case, sol_eig);
%!   assert_simple(sol_eig.f([7, 12, 19, 39]), fref, 5e-4 * max(fref));
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
