## fem_pre_mesh_import.m:359
%!test
%! ### TEST 359 - 1D dissipative wave equation with perfectly matched layers
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
%!   lambda = 100e-3;
%!   deltaPML = 2 * lambda;
%!   l1 = lambda;
%!   l2 = l1 + deltaPML;
%!   c = 1440;
%!   alpha = 1 + 2j;
%!   sigmax = alpha * c / deltaPML;
%!   rho = 1000;
%!   f = c / lambda;
%!   eta = 100;
%!   zeta = 100;
%!   omega = 2 * pi * f;
%!   dx = lambda / 30;
%!   w = dx;
%!   h = dx;
%!   vx0 = 1 + 0.5j;
%!   tau = (4/3 * eta + zeta) / (rho * c^2);
%!   fx =  1 / (1 - 1j * sigmax / omega);
%!   k = omega / (c * fx * sqrt(1i * omega * tau + 1));
%!   z = (c * rho) / sqrt(1i * omega * tau + 1);
%!   vxref = @(x, t) vx0 * exp(1j * (omega * t - k * x));
%!   pref = @(x, t) z * vxref(x, t);
%!   e1 = [0.1; 1; -1];
%!   e2 = [-1; 0.8; 0.2];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   e1 /= norm(e1);
%!   e2 /= norm(e2);
%!   e3 /= norm(e3);
%!   R = [e1, e2, e3];
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l1=%g;\n", l1);
%!   fprintf(fd, "l2=%g;\n", l2);
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx1 = %g;\n", dx);
%!   fprintf(fd, "dx2 = %g;\n", dx);
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
%!   fputs(fd, "Physical Surface(\"output\",2) = {tmp1[0]};\n");
%!   fputs(fd, "Physical Surface(\"impedance\",3) = {tmp2[0]};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp1[1]};\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
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
%!   opt_mesh.elem_type = {"iso27", "quad9"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   [x, idx] = sort(mesh.nodes(:, 1));
%!   mesh.nodes(:, 1:3) = mesh.nodes(:, 1:3) * R.';
%!   grp_idx_input = find([mesh.groups.quad9.id] == 1);
%!   grp_idx_output = find([mesh.groups.quad9.id] == 2);
%!   grp_idx_impedance = find([mesh.groups.quad9.id] == 3);
%!   grp_idx_volume1 = find([mesh.groups.iso27.id] == 3);
%!   grp_idx_volume2 = find([mesh.groups.iso27.id] == 4);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_volume1).elements) = 1;
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_volume2).elements) = 1;
%!   mesh.elements.acoustic_boundary.quad9 = mesh.elements.quad9([[mesh.groups.quad9([grp_idx_input, grp_idx_output])].elements], :);
%!   mesh.materials.acoustic_boundary.quad9 = zeros(rows(mesh.elements.acoustic_boundary.quad9), 1, "int32");
%!   mesh.materials.acoustic_boundary.quad9(mesh.groups.quad9(grp_idx_input).elements) = 1;
%!   mesh.materials.acoustic_boundary.quad9(mesh.groups.quad9(grp_idx_output).elements) = 1;
%!   mesh.elements.particle_velocity.quad9.nodes = mesh.elements.quad9(mesh.groups.quad9(grp_idx_input).elements, :);
%!   mesh.materials.particle_velocity.quad9 = ones(rows(mesh.elements.particle_velocity.quad9.nodes), 1, "int32");
%!   mesh.elements.acoustic_impedance.quad9.nodes = mesh.elements.quad9(mesh.groups.quad9(grp_idx_impedance).elements, :);
%!   mesh.elements.acoustic_impedance.quad9.z = repmat(z * fx^2, size(mesh.elements.acoustic_impedance.quad9.nodes));
%!   mesh.materials.acoustic_impedance.quad9 = ones(rows(mesh.elements.acoustic_impedance.quad9.nodes), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.quad9.vn = repmat(-real(vx0 / fx^2), numel(mesh.groups.quad9(grp_idx_input).elements), 9);
%!   load_case(2).particle_velocity.quad9.vn = repmat(-imag(vx0 / fx^2), numel(mesh.groups.quad9(grp_idx_input).elements), 9);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.material_data.eta = eta;
%!   mesh.material_data.zeta = zeta;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   mesh.elements.perfectly_matched_layers.iso27.f = zeros(3, columns(mesh.elements.iso27), rows(mesh.elements.iso27));
%!   mesh.elements.perfectly_matched_layers.iso27.e1 = zeros(3, columns(mesh.elements.iso27), rows(mesh.elements.iso27));
%!   mesh.elements.perfectly_matched_layers.iso27.e2 = zeros(3, columns(mesh.elements.iso27), rows(mesh.elements.iso27));
%!
%!   for i=1:3
%!     mesh.elements.perfectly_matched_layers.iso27.e1(i, :, :) = e1(i);
%!     mesh.elements.perfectly_matched_layers.iso27.e2(i, :, :) = e2(i);
%!     mesh.elements.perfectly_matched_layers.iso27.f(i, :, :) = 1 ./ (1 - 1j * sigmax / omega);
%!   endfor
%!
%!   [mat_ass.Ka_re, ...
%!    mat_ass.Ka_im, ...
%!    mat_ass.Da_re, ...
%!    mat_ass.Da_im, ...
%!    mat_ass.Ma_re, ...
%!    mat_ass.Ma_im, ...
%!    mat_ass.Ra, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_STIFFNESS_ACOUSTICS_IM, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_IM, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * complex(mat_ass.Ma_re, mat_ass.Ma_im) + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + complex(mat_ass.Ka_re, mat_ass.Ka_im);
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   sol.t = Psi / omega;
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof) * exp(1j * Psi));
%!   solC.Phi = Phi(dof_map.ndof);
%!   solC.PhiP = 1j * omega * Phi(dof_map.ndof);
%!   [solC.particle_velocity, ...
%!    solC.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                              dof_map, ...
%!                                              [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                               FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                              load_case, ...
%!                                              solC);
%!   mesh.nodes(:, 1:3) = mesh.nodes(:, 1:3) * R;
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(x, sol.p(idx, i), "-;p;r");
%!     plot(x, real(pref(x, sol.t(i))), "-;pref;k");
%!     ylim([min(min(sol.p)), max(max(sol.p))]);
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   tol = 0.5e-2;
%!   assert_simple(sol.p(idx, :), real(pref(x, sol.t)), tol * max(max(abs(pref(x, Psi)))));
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     vx = zeros(rows(mesh.nodes), 1);
%!     for j=1:3
%!       vx(mesh.elements.iso27(:)) += real(solC.particle_velocity.v.iso27(:, :, j)(:) * exp(1j * Psi(i))) * e1(j);
%!     endfor
%!     assert_simple(vx(idx), real(vxref(x, sol.t(i))), tol * abs(vx0));
%!     plot(x, vx(idx), "-;vn;r");
%!     plot(x, real(vxref(x, sol.t(i))), "-;Phi;k");
%!     ylim([-abs(vx0), abs(vx0)]);
%!     xlabel("x [m]");
%!     ylabel("vx [m/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("velocity distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
