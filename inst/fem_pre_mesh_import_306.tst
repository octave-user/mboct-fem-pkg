## fem_pre_mesh_import.m:306
%!test
%! ### TEST 306 - 1D wave equation with energy dissipation
%! do_plot = false;
%! if (do_plot)
%! close all;
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
%!   c = 1400;
%!   rho = 840;
%!   k = 400;
%!   eta = 2e3;
%!   zeta = 1e3;
%!   omega = k * c;
%!   f = omega / (2 * pi);
%!   lambda = c / f;
%!   dx = lambda / 60;
%!   w = dx;
%!   h = dx;
%!   A = 1;
%!   tau = (1 / (rho * c^2)) * ((4 / 3) * eta + zeta);
%!   k = omega / (c * sqrt(1i * omega * tau + 1));
%!   z = c * rho * sqrt(1i * omega * tau + 1);
%!   pref = @(x, t) -1i * A * omega * rho * exp(1i * (omega * t - (omega * x) / (c * sqrt(1i * omega * tau + 1))));
%!   vxref = @(x, t) -(1i * A * omega * sqrt(1i * omega * tau + 1)*exp(1i * (omega * t - (omega * x) / (c * sqrt(1i * omega * tau + 1)))))/c;
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
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!   fputs(fd, "Physical Surface(\"output\",2) = {tmp[0]};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=1;\n");
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
%!   opt_mesh.elem_type = {"tet20", "tria10"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   grp_idx_input = find([mesh.groups.tria10.id] == 1);
%!   grp_idx_output = find([mesh.groups.tria10.id] == 2);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   e1 = [0.5; -0.3; 0.7];
%!   e2 = [0.1; 0.6; 0.4];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   for i=1:columns(R)
%!     R(:, i) /= norm(R(:, i));
%!   endfor
%!   mesh.nodes(:, 1:3) = mesh.nodes(:, 1:3) * R.';
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   mesh.material_data.eta = eta;
%!   mesh.material_data.zeta = zeta;
%!   mesh.elements.particle_velocity.tria10.nodes = mesh.elements.tria10(mesh.groups.tria10(grp_idx_input).elements, :);
%!   mesh.materials.particle_velocity.tria10 = ones(rows(mesh.elements.particle_velocity.tria10.nodes), 1, "int32");
%!   mesh.elements.acoustic_impedance.tria10.nodes = mesh.elements.tria10(mesh.groups.tria10(grp_idx_output).elements, :);
%!   mesh.elements.acoustic_impedance.tria10.z = repmat(z, size(mesh.elements.acoustic_impedance.tria10.nodes));
%!   mesh.materials.acoustic_impedance.tria10 = ones(rows(mesh.elements.acoustic_impedance.tria10.nodes), 1, "int32");
%!   v_input = repmat(vxref(0, 0), size(mesh.elements.particle_velocity.tria10.nodes));
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.tria10.vn = real(-v_input);
%!   load_case(2).particle_velocity.tria10.vn = imag(-v_input);
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.Da_re, ...
%!    mat_ass.Da_im, ...
%!    mat_ass.Ra, ...
%!    mat_ass.n, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS, ...
%!                                         FEM_VEC_SURFACE_NORMAL_VECTOR], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * mat_ass.Ma + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + mat_ass.Ka;
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   x = mesh.nodes(:, 1:3) * R(:, 1);
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.Phi = real(Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.PhiP = real(1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.t = Psi / omega;
%!   [sol.particle_velocity, ...
%!    sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                             dof_map, ...
%!                                             [FEM_VEC_PARTICLE_VELOCITY, ...
%!                                              FEM_SCA_ACOUSTIC_INTENSITY], ...
%!                                             load_case, ...
%!                                             sol);
%!   vx = zeros(rows(mesh.nodes), numel(sol.t));
%!   for i=1:columns(mesh.elements.tet20)
%!     vx(mesh.elements.tet20(:, i), :) = 0;
%!     for j=1:3
%!       vx(mesh.elements.tet20(:, i), :) += reshape(sol.particle_velocity.v.tet20(:, i, j, :) * R(j, 1), rows(mesh.elements.tet20), numel(sol.t));
%!     endfor
%!   endfor
%!   vx = zeros(rows(mesh.nodes), numel(sol.t));
%!   nvx = zeros(rows(mesh.nodes), 1);
%!   for i=1:columns(mesh.elements.tet20)
%!     ++nvx(mesh.elements.tet20(:, i));
%!     for j=1:3
%!       vx(mesh.elements.tet20(:, i), :) += reshape(sol.particle_velocity.v.tet20(:, i, j, :) * R(j, 1), rows(mesh.elements.tet20), numel(sol.t));
%!     endfor
%!   endfor
%!   for i=1:columns(vx)
%!     vx(:, i) ./= nvx;
%!   endfor
%!   [x, idx] = sort(x);
%!   pref_ = real(pref(x, sol.t));
%!   vxref_ = real(vxref(x, sol.t));
%!   if (do_plot)
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(x, sol.p(idx, i), "-;p;r");
%!     plot(x, pref_(:, i), "-;pref;k");
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("pressure distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   for i=1:numel(Psi)
%!     figure("visible", "off");
%!     hold on;
%!     plot(x, vx(idx, i), "-;vx;r");
%!     plot(x, vxref_(:, i), "-;vxref;k");
%!     xlabel("x [m]");
%!     ylabel("vx [m/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("velocity distribution Psi=%.1fdeg", Psi(i) * 180 / pi));
%!   endfor
%!   endif
%!   tol = 1e-3;
%!   assert_simple(vx(idx, :), vxref_, tol * max(max(abs(vxref_))));
%!   assert_simple(sol.p(idx, :), pref_, tol * max(max(abs(pref_))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
