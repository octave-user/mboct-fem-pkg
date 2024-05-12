## fem_pre_mesh_import.m:149
%!test
%! try
%! ## TEST 149 - 2D wave equation
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5
%! ####################################################
%! close all;
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
%!     l = 100e-3;
%!     c = 340;
%!     rho = 1.25;
%!     k = 100;
%!     omega = k * c;
%!     f = omega / (2 * pi);
%!     lambda = c / f;
%!     dx = lambda / 100;
%!     w = l;
%!     h = dx;
%!     Phi = 30 * pi / 180;
%!     A = 100 + 40j;
%!     p = @(t, x, y) A * exp(1j * (omega * t - k * x * cos(Phi) - k * y * sin(Phi)));
%!     nv = @(t, x, y, n) (A*(n(2)*sin(Phi)+n(1)*cos(Phi))*k*exp(-1i*sin(Phi)*k*y-1i*cos(Phi)*k*x+1i*omega*t))/(omega*rho);
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
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\",3) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"s1\",1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"s2\",2) = {tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"s3\",3) = {tmp[0]};\n");
%!     fputs(fd, "Physical Surface(\"s4\",4) = {tmp[3]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.5;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"iso8", "iso4"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_s1 = find([mesh.groups.iso4.id] == 1);
%!   grp_idx_s2 = find([mesh.groups.iso4.id] == 2);
%!   grp_idx_s3 = find([mesh.groups.iso4.id] == 3);
%!   grp_idx_s4 = find([mesh.groups.iso4.id] == 4);
%!   elno_s1 = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s1).elements, :);
%!   elno_s2 = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s2).elements, :);
%!   elno_s3 = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s3).elements, :);
%!   elno_s4 = mesh.elements.iso4(mesh.groups.iso4(grp_idx_s4).elements, :);
%!   x1 = mesh.nodes(:, 1)(elno_s1);
%!   y1 = mesh.nodes(:, 2)(elno_s1);
%!   x2 = mesh.nodes(:, 1)(elno_s2);
%!   y2 = mesh.nodes(:, 2)(elno_s2);
%!   x3 = mesh.nodes(:, 1)(elno_s3);
%!   y3 = mesh.nodes(:, 2)(elno_s3);
%!   x4 = mesh.nodes(:, 1)(elno_s4);
%!   y4 = mesh.nodes(:, 2)(elno_s4);
%!   n1 = [-1; 0];
%!   n2 = [0; -1];
%!   n3 = [1; 0];
%!   n4 = [0; 1];
%!   t = 0;
%!   nv1 = nv(t, x1, y1, n1);
%!   nv2 = nv(t, x2, y2, n2);
%!   nv3 = nv(t, x3, y3, n3);
%!   nv4 = nv(t, x4, y4, n4);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.elements.particle_velocity.iso4.nodes = [elno_s1; elno_s2; elno_s3; elno_s4];
%!   mesh.materials.particle_velocity.iso4 = ones(rows(mesh.elements.particle_velocity.iso4.nodes), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.iso4.vn = real([nv1; nv2; nv3; nv4]);
%!   load_case(2).particle_velocity.iso4.vn = imag([nv1; nv2; nv3; nv4]);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.Da, ...
%!    mat_ass.Ra, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = complex(-omega^2 * mat_ass.Ma + 1j * omega * mat_ass.Da + mat_ass.Ka);
%!   Reff = complex(mat_ass.Ra(:, 1) + 1j * mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   sol.t = Psi / omega;
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * omega * sol.t));
%!   pref = real(p(sol.t, mesh.nodes(:, 1), mesh.nodes(:, 2)));
%!   tol = 2e-2;
%!   sol2.p = pref;
%!   sol2.t = sol.t;
%!   assert_simple(sol.p, pref, tol * max(max(abs(pref))));
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
