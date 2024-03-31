## fem_pre_mesh_import.m:132
%!test
%! ### TEST 132 - 1D wave equation
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5
%! ####################################################
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
%!   dx = lambda / 40;
%!   w = dx;
%!   h = dx;
%!   A = -70;
%!   B = 50;
%!   ## /* Maxima script */
%!   ## kill(all)$
%!   ## k:omega/c$
%!   ## T:2 * %pi / omega$
%!   ## pp:((A + %i * B) / 2 * exp(%i * (omega * t - k * x)))$
%!   ## pm:((A - %i * B) / 2 * exp(-%i * (omega * t - k * x)))$
%!   ## p:pp + pm$
%!   ## Phip:-integrate(pp, t)$
%!   ## Phim:-integrate(pm, t)$
%!   ## vp:diff(Phip, x) / rho$
%!   ## vm:diff(Phim, x) / rho$
%!   ## p:pp + pm$
%!   ## v:vp + vm$
%!   ## I:p * v;
%!   ## Iavg:integrate(I, t, 0, T) / T;
%!   ## p:ratsimp(rectform(p));
%!   ## I:factor(trigreduce(rectform(I)));
%!   ## v:factor(trigreduce(rectform(v)));
%!   ## Iavg:factor(trigreduce(rectform(Iavg)));
%!   I = @(x, t) (2*A*B*sin((2*omega*(x-c*t))/c)-B^2*cos((2*omega*(x-c*t))/c)+A^2*cos((2*omega*(x-c*t))/c)+B^2+A^2)/(2*c*rho);
%!   Iavg = (B^2+A^2)/(2*c*rho);
%!   v = @(x, t) (B*sin((omega*(x-c*t))/c)+A*cos((omega*(x-c*t))/c))/(c*rho);
%!   p = @(x, t) B*sin((omega*x-c*omega*t)/c)+A*cos((omega*x-c*omega*t)/c);
%!   vp = @(x, t) ((i*B+A)*e^(1i*(omega*t-(omega*x)/c)))/(2*c*rho);
%!   pp = @(x, t) ((i*B+A)*e^(1i*(omega*t-(omega*x)/c)))/2;
%!   z = rho * c;
%!   pinput = 2 * pp(0, 0);
%!   vx0 = 2 * vp(0, 0);
%!   vxl = 2 * vp(l, 0);
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
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Physical Volume(\"volume\",3) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!   fputs(fd, "Physical Surface(\"output\",2) = {tmp[0]};\n");
%!   fputs(fd, "Mesh.HighOrderIterMax = 200;\n");
%!   fputs(fd, "Mesh.HighOrderPassMax = 100;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!   fputs(fd, "Mesh.OptimizeThreshold = 0.99;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
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
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_input = find([mesh.groups.tria6h.id] == 1);
%!   grp_idx_output = find([mesh.groups.tria6h.id] == 2);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria6h = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.acoustic_boundary.tria6h = ones(rows(mesh.elements.acoustic_boundary.tria6h), 1, "int32");
%!   mesh.elements.particle_velocity.tria6h.nodes = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.particle_velocity.tria6h = ones(rows(mesh.elements.particle_velocity.tria6h.nodes), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.tria6h.vn = [repmat(-real(vx0), numel(mesh.groups.tria6h(grp_idx_input).elements), 6);
%!                                               repmat(real(vxl), numel(mesh.groups.tria6h(grp_idx_output).elements), 6)];
%!   load_case(2).particle_velocity.tria6h.vn = [repmat(-imag(vx0), numel(mesh.groups.tria6h(grp_idx_input).elements), 6);
%!                                               repmat(imag(vxl), numel(mesh.groups.tria6h(grp_idx_output).elements), 6)];
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
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   x = mesh.nodes(:, 1);
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.Phi = real(Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.PhiP = real(1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   [sol.particle_velocity, ...
%!    sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                             dof_map, ...
%!                                             [FEM_VEC_PARTICLE_VELOCITY, ...
%!                                              FEM_SCA_ACOUSTIC_INTENSITY], ...
%!                                             load_case, ...
%!                                             sol);
%!   solC.Phi = Phi(dof_map.ndof, :);
%!   solC.PhiP = 1j * omega * Phi(dof_map.ndof, :);
%!   [solC.particle_velocity, ...
%!    solC.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                              dof_map, ...
%!                                              [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                               FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                              load_case, ...
%!                                              solC);
%!   solC.f = omega / (2 * pi);
%!   sol.t = Psi / omega;
%!   pref = p(mesh.nodes(:, 1), sol.t);
%!   vnref = v(mesh.nodes(:, 1), sol.t);
%!   Iref = I(mesh.nodes(:, 1), sol.t);
%!   Pref = w * h * Iref;
%!   [~, idx] = sort(mesh.nodes(:, 1));
%!   if (do_plot)
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(mesh.nodes(idx, 1), sol.p(idx, i), "-;p;r");
%!     plot(mesh.nodes(idx, 1), pref(idx, i), "-;pref;k");
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     vx = vxC = zeros(rows(mesh.nodes), numel(sol.t));
%!     vx(mesh.elements.tet10h(:)) = sol.particle_velocity.v.tet10h(:, :, 1, i)(:);
%!     vxC(mesh.elements.tet10h(:)) = real(solC.particle_velocity.v.tet10h(:, :, 1)(:) * exp(1j * Psi(i)));
%!     plot(mesh.nodes(idx, 1), vx(idx), "-;p;r");
%!     plot(mesh.nodes(idx, 1), vxC(idx), "-;pC;g");
%!     plot(mesh.nodes(idx, 1), vnref(idx, i), "-;pref;k");
%!     xlabel("x [m]");
%!     ylabel("vx [m/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("velocity distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   figure("visible", "off");
%!   hold on;
%!   elem_id = mesh.groups.tria6h(grp_idx_output).elements;
%!   elem_no = mesh.elements.tria6h(elem_id, 1);
%!   plot(sol.t, sum(sol.acoustic_intensity.P.tria6h(elem_id, :), 1)(:), "-;P;r");
%!   plot(sol.t, mean(Pref(elem_no, :), 1).', "-;Pref;k");
%!   xlabel("t [s]");
%!   ylabel("P [W]");
%!   grid on;
%!   grid minor on;
%!   title("instantaneous output sound power versus time");
%!   endif
%!   tol = 1e-3;
%!   assert_simple(sol.p, pref, tol * max(max(abs(pref))));
%!   tol = 5e-3;
%!   for i=1:numel(mesh.groups.tria6h(grp_idx_input).elements)
%!     elem_id = mesh.groups.tria6h(grp_idx_input).elements(i);
%!     for j=1:columns(sol.particle_velocity.vn.tria6h)
%!       elem_no = mesh.elements.tria6h(elem_id, j);
%!       assert_simple(sol.particle_velocity.vn.tria6h(elem_id, j, :)(:), -vnref(elem_no, :).', tol * max(max(abs(vnref))));
%!       assert_simple(sol.acoustic_intensity.I.tria6h(elem_id, j, :)(:), -Iref(elem_no, :).', tol * max(max(abs(Iref))));
%!       assert_simple(real(solC.particle_velocity.vn.tria6h(elem_id, j) * exp(1j * omega * sol.t).'), -vnref(elem_no, :).', tol * max(max(abs(vnref))));
%!     endfor
%!   endfor
%!   elem_id = mesh.groups.tria6h(grp_idx_input).elements;
%!   elem_no = mesh.elements.tria6h(elem_id, 1);
%!   assert_simple(sum(sol.acoustic_intensity.P.tria6h(elem_id, :), 1)(:), -mean(Pref(elem_no, :),1).', tol * max(max(abs(Pref))));
%!   for i=1:numel(mesh.groups.tria6h(grp_idx_output).elements)
%!     elem_id = mesh.groups.tria6h(grp_idx_output).elements(i);
%!     for j=1:columns(sol.particle_velocity.vn.tria6h)
%!       elem_no = mesh.elements.tria6h(elem_id, j);
%!       assert_simple(sol.particle_velocity.vn.tria6h(elem_id, j, :)(:), vnref(elem_no, :).', tol * max(max(abs(vnref))));
%!       assert_simple(sol.acoustic_intensity.I.tria6h(elem_id, j, :)(:), Iref(elem_no, :).', tol * max(max(abs(Iref))));
%!       assert_simple(real(solC.particle_velocity.vn.tria6h(elem_id, j) * exp(1j * omega * sol.t).'), vnref(elem_no, :).', tol * max(max(abs(vnref))));
%!     endfor
%!   endfor
%!   elem_id = mesh.groups.tria6h(grp_idx_output).elements;
%!   elem_no = mesh.elements.tria6h(elem_id, 1);
%!   assert_simple(sum(sol.acoustic_intensity.P.tria6h(elem_id, :), 1)(:), mean(Pref(elem_no, :), 1).', tol * max(max(abs(Pref))));
%!   assert_simple(max(max(abs(solC.acoustic_intensity.I.tria6h - Iavg))) < tol * abs(Iavg));
%!   for i=1:numel(mesh.groups.tria6h)
%!     assert_simple(sum(solC.acoustic_intensity.P.tria6h(mesh.groups.tria6h(i).elements)), Iavg * w * h, tol * abs(Iavg * w * h));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
