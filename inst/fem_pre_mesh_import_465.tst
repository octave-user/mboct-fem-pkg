## fem_pre_mesh_import.m: 465
%!test
%! try
%! ### TEST 465
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
%!     l = 10e-3;
%!     w = 15e-3;
%!     h = 12e-3;
%!     dx = 2e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "l=%g;\n", l);
%!     fprintf(fd, "w=%g;\n", w);
%!     fprintf(fd, "h=%g;\n", h);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {0,w,0,dx};\n");
%!     fputs(fd, "Point(3) = {0,w,h,dx};\n");
%!     fputs(fd, "Point(4) = {0,0,h,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {l,0,0} {\n");
%!     fputs(fd, "  Surface{6}; Layers{Ceil(l/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"convection\",1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"source\",2) = {tmp[0]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold = 0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"penta18", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   lambda = 50;
%!   mesh.materials.penta18 = ones(rows(mesh.elements.penta18), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh.material_data.cp = 465;
%!   thetae = 100;
%!   he = lambda / l;
%!   q = 1000000;
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   mesh.elements.convection.tria6h.nodes = mesh.elements.tria6h([mesh.groups.tria6h(1).elements], :);
%!   mesh.elements.convection.tria6h.h = repmat(he, rows(mesh.elements.convection.tria6h.nodes), columns(mesh.elements.convection.tria6h.nodes));
%!   load_case.convection.tria6h.theta = repmat(thetae, rows(mesh.elements.convection.tria6h.nodes), columns(mesh.elements.convection.tria6h.nodes));
%!   load_case.heat_source.tria6h.nodes = mesh.elements.tria6h([mesh.groups.tria6h(2).elements], :);
%!   load_case.heat_source.tria6h.q = repmat(q, size(load_case.heat_source.tria6h.nodes));
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
%!   thetas = thetae + 3 * l * q / lambda;
%!   sol.theta = fem_sol_factor(mat_ass.Kk) \ mat_ass.Qc;
%!   x = mesh.nodes(:, 1:3) * R(:, 1);
%!   theta_ref = (x + l) / (3 * l) * (thetas - thetae) + thetae;
%!   assert_simple(sol.theta, theta_ref, eps^0.8 * max(abs(thetae)));
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
