## fem_pre_mesh_import.m:135
%!test
%! try
%! ### TEST 135 - 1D wave equation
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
%!   [fd, msg] = fopen([filename, ".geo"], "w");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   l = 100e-3;
%!   c = 340;
%!   rho = 1.25;
%!   k = 350;
%!   omega = k * c;
%!   f = omega / (2 * pi);
%!   lambda = c / f;
%!   dx = lambda / 30;
%!   w = dx;
%!   h = dx;
%!   pinput = 50 + 100j;
%!   poutput = pinput * exp(-1j * k * l);
%!   vx0 = k * pinput / (omega * rho);
%!   vxl = k * pinput / (omega * rho) * exp(-1j * k * l);
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l=%g;\n", l);
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx = %g;\n", dx);
%!   fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!   fputs(fd, "Point(2) = {0,w,0,dx};\n");
%!   fputs(fd, "Point(3) = {0,w,h,dx};\n");
%!   fputs(fd, "Point(4) = {0,0,h,dx};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {l,0,0} {\n");
%!   fputs(fd, "  Surface{6}; Layers{Ceil(l/dx)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface {6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",3) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!   fputs(fd, "Physical Surface(\"output\",2) = {tmp[0]};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
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
%!   opt_mesh.elem_type = {"iso20r", "quad8r"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_input = find([mesh.groups.quad8r.id] == 1);
%!   grp_idx_output = find([mesh.groups.quad8r.id] == 2);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.elements.particle_velocity.quad8r.nodes = mesh.elements.quad8r([[mesh.groups.quad8r([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.particle_velocity.quad8r = ones(rows(mesh.elements.particle_velocity.quad8r.nodes), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.quad8r.vn = [repmat(-real(vx0), numel(mesh.groups.quad8r(grp_idx_input).elements), 8);
%!                                               repmat(real(vxl), numel(mesh.groups.quad8r(grp_idx_output).elements), 8)];
%!   load_case(2).particle_velocity.quad8r.vn = [repmat(-imag(vx0), numel(mesh.groups.quad8r(grp_idx_input).elements), 8);
%!                                               repmat(imag(vxl), numel(mesh.groups.quad8r(grp_idx_output).elements), 8)];
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
%!   # opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = complex(-omega^2 * mat_ass.Ma + 1j * omega * mat_ass.Da + mat_ass.Ka);
%!   Reff = complex(mat_ass.Ra(:, 1) + 1j * mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   x = mesh.nodes(:, 1);
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.t = Psi / omega;
%!   pref = real(pinput * exp(1j * (omega * sol.t - k * x)));
%!   [x, idx] = sort(mesh.nodes(:, 1));
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(x, sol.p(idx, i), "-;p;r");
%!     plot(x, pref(idx, i), "-;pref;k");
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   tol = 1e-3;
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
