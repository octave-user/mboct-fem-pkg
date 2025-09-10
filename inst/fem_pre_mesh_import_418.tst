## fem_pre_mesh_import.m:178
%!test
%! try
%! ### TEST 178 - 1D wave equation with medium interface
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
%!   l1 = 10e-3;
%!   l2 = 20e-3;
%!   c1 = 340;
%!   rho1 = 1.25;
%!   c2 = 1400;
%!   rho2 = 1000;
%!   k2 = 400;
%!   omega = k2 * c2;
%!   f = omega / (2 * pi);
%!   lambda2 = c2 / f;
%!   dx2 = lambda2 / 100;
%!   k1 = omega / c1;
%!   lambda1 = c1 / f;
%!   dx1 = lambda1 / 100;
%!   w = min([dx1, dx2]);
%!   h = min([dx1, dx2]);
%!   vx0 = 1;
%!   zI = rho1 * c1;
%!   zR = -rho1 * c1;
%!   zT = rho2 * c2;
%!   AI = (e^(2*i*k1*l1)*vx0*(zT-zR))/((e^(2*i*k1*l1)-1)*zT-e^(2*i*k1*l1)*zR+zI);
%!   AR = -(vx0*(zT-zI))/((e^(2*i*k1*l1)-1)*zT-e^(2*i*k1*l1)*zR+zI);
%!   AT = -(e^(i*k2*l1+i*k1*l1)*vx0*(zR-zI))/((e^(2*i*k1*l1)-1)*zT-e^(2*i*k1*l1)*zR+zI);
%!   vI = @(x, t) AI*exp(i*(omega*t-k1*x));
%!   vR = @(x, t) AR*exp(i*(omega*t+k1*x));
%!   vT = @(x, t) AT*exp(i*(omega*t-k2*x));
%!   v1x = @(x, t) vI(x, t) + vR(x, t);
%!   v2x = @(x, t) vT(x, t);
%!   p1 = @(x, t) zI * vI(x, t) + zR * vR(x, t);
%!   p2 = @(x, t) zT * vT(x, t);
%!   pref = @(x, t) p1(x, t) .* (x <= l1) + p2(x, t) .* (x > l1);
%!   vxref = @(x, t) v1x(x, t) .* (x <= l1) + v2x(x, t) .* (x > l1);
%!   vx0 = vxref(0, 0);
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l1=%g;\n", l1);
%!   fprintf(fd, "l2=%g;\n", l2);
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx1 = %g;\n", dx1);
%!   fprintf(fd, "dx2 = %g;\n", dx2);
%!   fputs(fd, "Point(1) = {0,0,0};\n");
%!   fputs(fd, "Point(2) = {0,w,0};\n");
%!   fputs(fd, "Point(3) = {0,w,h};\n");
%!   fputs(fd, "Point(4) = {0,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp1[] = Extrude {l1,0,0} {\n");
%!   fputs(fd, "  Surface{6}; Layers{Ceil(l1/dx1)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "tmp2[] = Extrude {l2 - l1,0,0} {\n");
%!   fputs(fd, "  Surface{tmp1[0]}; Layers{Ceil((l2 - l1)/dx2)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume1\",3) = {tmp1[1]};\n");
%!   fputs(fd, "Physical Volume(\"volume2\",4) = {tmp2[1]};\n");
%!   fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!   fputs(fd, "Physical Surface(\"output\",2) = {tmp2[0]};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp1[1],tmp2[1]};\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   grp_idx_volume1 = find([mesh.groups.iso20r.id] == 3);
%!   grp_idx_volume2 = find([mesh.groups.iso20r.id] == 4);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_volume1).elements) = 1;
%!   mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_volume2).elements) = 2;
%!   mesh.elements.acoustic_boundary.quad8r = mesh.elements.quad8r([[mesh.groups.quad8r([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.acoustic_boundary.quad8r = zeros(rows(mesh.elements.acoustic_boundary.quad8r), 1, "int32");
%!   mesh.materials.acoustic_boundary.quad8r(mesh.groups.quad8r(grp_idx_input).elements) = 1;
%!   mesh.materials.acoustic_boundary.quad8r(mesh.groups.quad8r(grp_idx_output).elements) = 2;
%!   mesh.elements.particle_velocity.quad8r.nodes = mesh.elements.quad8r(mesh.groups.quad8r(grp_idx_input).elements, :);
%!   mesh.materials.particle_velocity.quad8r = ones(rows(mesh.elements.particle_velocity.quad8r.nodes), 1, "int32");
%!   mesh.elements.acoustic_impedance.quad8r.nodes = mesh.elements.quad8r(mesh.groups.quad8r(grp_idx_output).elements, :);
%!   mesh.elements.acoustic_impedance.quad8r.z = repmat(zT, size(mesh.elements.acoustic_impedance.quad8r.nodes));
%!   mesh.materials.acoustic_impedance.quad8r = repmat(int32(2), rows(mesh.elements.acoustic_impedance.quad8r.nodes), 1);
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.quad8r.vn = repmat(-real(vx0), numel(mesh.groups.quad8r(grp_idx_input).elements), 8);
%!   load_case(2).particle_velocity.quad8r.vn = repmat(-imag(vx0), numel(mesh.groups.quad8r(grp_idx_input).elements), 8);
%!   mesh.material_data = struct("rho", {rho1, rho2}, "c", {c1, c2});
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.Da_re, ...
%!    mat_ass.Da_im, ...
%!    mat_ass.Ra, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pardiso";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = complex(-omega^2 * mat_ass.Ma + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + mat_ass.Ka);
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
%!   solC.Phi = complex(Phi(dof_map.ndof, :));
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
%!   [~, idx] = sort(mesh.nodes(:, 1));
%!   vxreft = real(vxref(mesh.nodes(:, 1), sol.t));
%!   preft = real(pref(mesh.nodes(:, 1), sol.t));
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(mesh.nodes(idx, 1), sol.p(idx, i), "-;p;r");
%!     plot(mesh.nodes(idx, 1), preft(idx, i), "-;pref;k");
%!     ylim([min(min(sol.p)), max(max(sol.p))]);
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     vx = vxC = zeros(rows(mesh.nodes), 1);
%!     vx(mesh.elements.iso20r(:)) = sol.particle_velocity.v.iso20r(:, :, 1, i)(:);
%!     vxC(mesh.elements.iso20r(:)) = real(solC.particle_velocity.v.iso20r(:, :, 1)(:) * exp(1j * Psi(i)));
%!     plot(mesh.nodes(idx, 1), vx(idx), "-;vn;r");
%!     plot(mesh.nodes(idx, 1), vxC(idx), "-;vC;g");
%!     plot(mesh.nodes(idx, 1), vxreft(idx, i), "-;vref;k");
%!     ylim([min(min(min(min(sol.particle_velocity.v.iso20r)))), max(max(max(max(sol.particle_velocity.v.iso20r))))]);
%!     xlabel("x [m]");
%!     ylabel("vx [m/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("velocity distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   tol = 1e-4;
%!   assert_simple(sol.p, preft, tol * max(max(abs(preft))));
%!   tol = 1e-3;
%!   vx = zeros(rows(mesh.nodes), numel(sol.t));
%!   for i=1:numel(sol.t)
%!     vx(mesh.elements.iso20r(:), i) = sol.particle_velocity.v.iso20r(:, :, 1, i)(:);
%!   endfor
%!   assert_simple(vx, vxreft, tol * max(max(abs(vxreft))));
%!   figure("visible", "off");
%!   hold on;
%!   elem_id = mesh.groups.quad8r(grp_idx_output).elements;
%!   elem_no = mesh.elements.quad8r(elem_id, 1);
%!   plot(sol.t, sum(sol.acoustic_intensity.P.quad8r(elem_id, :), 1)(:), "-;P;r");
%!   xlabel("t [s]");
%!   ylabel("P [W]");
%!   grid on;
%!   grid minor on;
%!   title("instantaneous output sound power versus time");
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
