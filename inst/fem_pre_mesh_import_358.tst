## fem_pre_mesh_import.m:358
%!test
%! ### TEST 358
%! ## Code_Aster
%! ## SDLS502 - Square plate "solid" simply supported
%! ## KUDAWOO Ayaovi-Dzifa
%! ## 03/08/2011
%! ## V2.03.502
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "w");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   do_plot = false;
%!   unit_meters = 1;
%!   unit_second = 1;
%!   unit_kilograms = 1;
%!   unit_newton = unit_kilograms * unit_meters / unit_second^2;
%!   unit_pascal = unit_newton / unit_meters^2;
%!   unit_watt = unit_newton * unit_meters / unit_second;
%!   a = 10000e-3 / unit_meters;
%!   b = 10000e-3 / unit_meters;
%!   h = 1000e-3 / unit_meters;
%!   E = 210000e6 / unit_pascal;
%!   rho = 8000 / (unit_kilograms / unit_meters^3);
%!   nu = 0.3;
%!   alpha = 0 / (unit_second^-1);
%!   beta = 0 / unit_second;
%!   dx = h / 2;
%!   fref = [44.762 110.52 110.52 169.08 193.93 206.64];
%!   N = numel(fref) + 3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx=%g;\n", dx);
%!   fputs(fd, "Point(1) = {-a/2, -b/2, -h/2, dx};\n");
%!   fputs(fd, "Point(2) = {-a/2, -b/2,    0, dx};\n");
%!   fputs(fd, "Point(3) = {-a/2, -b/2,  h/2, dx};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "A1[] = Extrude{a,0,0} {\n");
%!   fputs(fd, "  Line{1}; Layers{Ceil(a / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "A2[] = Extrude{a,0,0} {\n");
%!   fputs(fd, "  Line{2}; Layers{Ceil(a / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V1[] = Extrude{0,b,0}{\n");
%!   fputs(fd, "  Surface{A1[1]}; Layers{Ceil(b / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "V2[] = Extrude{0,b,0}{\n");
%!   fputs(fd, "  Surface{A2[1]}; Layers{Ceil(b / dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Coherence;\n");
%!   fputs(fd, "Physical Volume(\"v1\",1)={V1[1],V2[1]};\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_v1_iso27 = find([mesh.groups.iso27.id] == 1);
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_v1_iso27).elements) = 1;
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.alpha = alpha;
%!   mesh.material_data.beta = beta;
%!   tol = sqrt(eps) * max([a,b,h]);
%!   node_id_support1 = find((abs(mesh.nodes(:, 1) - a/2) < tol) & (abs(mesh.nodes(:,3)) < tol));
%!   node_id_support2 = find((abs(mesh.nodes(:, 1) + a/2) < tol) & (abs(mesh.nodes(:,3)) < tol));
%!   node_id_support3 = find((abs(mesh.nodes(:, 2) - b/2) < tol) & (abs(mesh.nodes(:,3)) < tol));
%!   node_id_support4 = find((abs(mesh.nodes(:, 2) + b/2) < tol) & (abs(mesh.nodes(:,3)) < tol));
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(node_id_support1, 3) = true;
%!   load_case.locked_dof(node_id_support2, 3) = true;
%!   load_case.locked_dof(node_id_support3, 3) = true;
%!   load_case.locked_dof(node_id_support4, 3) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS, ...
%!                                         FEM_MAT_MASS], ...
%!                                        load_case);
%!   opt_sol.rho = 1;
%!   opt_sol.tolerance = sqrt(eps);
%!   opt_sol.algorithm = "shift-invert";
%!   opt_sol.solver = "pardiso";
%!   opt_sol.pre_scaling = true;
%!   opt_sol.verbose = int32(0);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, N, opt_sol);
%!   tol = 2.5e-2;
%!   assert_simple(sol_eig.f(4:end), fref, tol * max(fref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
