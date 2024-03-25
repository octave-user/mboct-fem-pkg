## fem_pre_mesh_import.m:278
%!test
%! ### TEST 278
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
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
%!     a = 10e-3;
%!     b = 10e-3;
%!     c = 10e-3;
%!     d = 2e-3;
%!     dx = 2e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "d=%g;\n", d);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {-b,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {b,0,0,dx};\n");
%!     fputs(fd, "Point(3) = {0,c,0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,1};\n");
%!     fputs(fd, "Line Loop(4) = {1,2,3};\n");
%!     fputs(fd, "Plane Surface(5) = {4};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,d} {\n");
%!     fputs(fd, "  Surface{5};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"h0\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"h1\",2) = {tmp[3]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   lambda = 50;
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh.material_data.cp = 465;
%!   thetae = [100, 200];
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   mesh.elements.convection.tria10.nodes = mesh.elements.tria10([mesh.groups.tria10.elements], :);
%!   x = mesh.nodes(:, 1);
%!   h0 = 1e11;
%!   h1 = 1e11;
%!   theta_s = thetae(1) + (thetae(2) - thetae(1)) / (2 * (a + b)) * (x + a + b);
%!   mesh.elements.convection.tria10.h = [repmat(h0, size(mesh.elements.tria10(mesh.groups.tria10(1).elements, :)));
%!                                      repmat(h1, size(mesh.elements.tria10(mesh.groups.tria10(2).elements, :)))];
%!   load_case.convection.tria10.theta = [theta_s(mesh.elements.tria10(mesh.groups.tria10(1).elements, :));
%!                                      theta_s(mesh.elements.tria10(mesh.groups.tria10(2).elements, :))];
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   e1 = [1; 0.6; -0.3];
%!   e2 = [-0.5; -0.3; 0.8];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   e1 /= norm(e1);
%!   e2 /= norm(e2);
%!   e3 /= norm(e3);
%!   R = [e1, e2, e3];
%!   mesh.nodes = [mesh.nodes(:, 1:3) * R.', mesh.nodes(:, 4:6) * R.'];
%!   [mat_ass.Kk, ...
%!    mat_ass.Qc] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_THERMAL_COND, ...
%!                                  FEM_VEC_LOAD_THERMAL], ...
%!                              load_case);
%!   opt_sol.refine_max_iter = int32(20);
%!   sol.theta = fem_sol_factor(mat_ass.Kk, opt_sol) \ mat_ass.Qc;
%!   assert_simple(sol.theta, theta_s, eps^0.5 * max(abs(thetae)));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect