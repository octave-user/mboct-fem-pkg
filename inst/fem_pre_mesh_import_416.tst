## fem_pre_mesh_import.m:171
%!test
%! try
%! ### TEST 171 - spherical waves
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
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     c = 340;
%!     rho = 1.25;
%!     k = 200;
%!     omega = k * c;
%!     f = omega / (2 * pi);
%!     lambda = c / f;
%!     dx = lambda / 100;
%!     R = 100e-3;
%!     r = 20e-3;
%!     A = 100 + 30j;
%!     Phi1 = -1 * pi / 180;
%!     Phi2 = 1 * pi / 180;
%!     Phi3 = 2 * pi / 180;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R=%g;\n", R);
%!     fprintf(fd, "r=%g;\n", r);
%!     fprintf(fd, "Phi1=%g;\n", Phi1);
%!     fprintf(fd, "Phi2=%g;\n", Phi2);
%!     fprintf(fd, "Phi3=%g;\n", Phi3);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Phi4 = 0.5*(Phi1+Phi2)+Pi/2;\n");
%!     fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {r * Cos(Phi1),r * Sin(Phi1),0,dx};\n");
%!     fputs(fd, "Point(3) = {R * Cos(Phi1),R * Sin(Phi1),0,dx};\n");
%!     fputs(fd, "Point(4) = {R * Cos(Phi2),R * Sin(Phi2),0,dx};\n");
%!     fputs(fd, "Point(5) = {r * Cos(Phi2),r * Sin(Phi2),0,dx};\n");
%!     fputs(fd, "Line(1) = {2,3};\n");
%!     fputs(fd, "Circle(2) = {3,1,4};\n");
%!     fputs(fd, "Line(3) = {4,5};\n");
%!     fputs(fd, "Circle(4) = {5,1,2};\n");
%!     fputs(fd, "Curve Loop(5) = {1,2,3, 4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {{Cos(Phi4),Sin(Phi4),0},{0,0,0},Phi3} {\n");
%!     fputs(fd, "  Surface{6}; Layers{Ceil(r*Phi3/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"inlet\",2) = {tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"outlet\",3) = {tmp[3]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   opt_mesh.elem_type = {"quad8r", "penta15", "iso20r"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   if (isfield(mesh.elements, "penta15"))
%!      mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   endif
%!   if (isfield(mesh.elements, "iso20r"))
%!      mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   endif
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   grp_idx_inlet = find([mesh.groups.quad8r.id] == 2);
%!   grp_idx_outlet = find([mesh.groups.quad8r.id] == 3);
%!   mesh.elements.particle_velocity.quad8r.nodes = mesh.elements.quad8r(mesh.groups.quad8r(grp_idx_inlet).elements, :);
%!   mesh.materials.particle_velocity.quad8r = ones(rows(mesh.elements.particle_velocity.quad8r.nodes), 1, "int32");
%!   xi = mesh.nodes(:, 1)(mesh.elements.particle_velocity.quad8r.nodes);
%!   yi = mesh.nodes(:, 2)(mesh.elements.particle_velocity.quad8r.nodes);
%!   zi = mesh.nodes(:, 3)(mesh.elements.particle_velocity.quad8r.nodes);
%!   ri = sqrt(xi.^2 + yi.^2 + zi.^2);
%!   p = @(r, t) A * exp(1j * (omega * t - k * r)) ./ r; ## according equation 5.11.6
%!   z = @(r) rho * c * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according equation 5.11.9
%!   vn = @(r, t) -(1 - 1j ./ (k * r)) .* p(r, t) / (rho * c);
%!   mesh.elements.acoustic_impedance.quad8r.nodes = mesh.elements.quad8r(mesh.groups.quad8r(grp_idx_outlet).elements, :);
%!   xo = mesh.nodes(:, 1)(mesh.elements.acoustic_impedance.quad8r.nodes);
%!   yo = mesh.nodes(:, 2)(mesh.elements.acoustic_impedance.quad8r.nodes);
%!   zo = mesh.nodes(:, 3)(mesh.elements.acoustic_impedance.quad8r.nodes);
%!   ro = sqrt(xo.^2 + yo.^2 + zo.^2);
%!   mesh.elements.acoustic_impedance.quad8r.z = z(ro);
%!   mesh.materials.acoustic_impedance.quad8r = ones(rows(mesh.elements.acoustic_impedance.quad8r.nodes), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.quad8r.vn = real(vn(ri, 0));
%!   load_case(2).particle_velocity.quad8r.vn = imag(vn(ri, 0));
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
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
%!   Keff = -omega^2 * mat_ass.Ma + 1j * omega * complex(mat_ass.Da_re,  mat_ass.Da_im) + mat_ass.Ka;
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   x = mesh.nodes(:, 1);
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.t = Psi / omega;
%!   xn = mesh.nodes(:, 1);
%!   yn = mesh.nodes(:, 2);
%!   zn = mesh.nodes(:, 3);
%!   rn = sqrt(xn.^2 + yn.^2 + zn.^2);
%!   sol2.p = real(p(rn, sol.t));
%!   sol2.t = sol.t;
%!   [rn, idx] = sort(rn);
%!   if (do_plot)
%!   for j=1:columns(sol.p)
%!     figure("visible", "off");
%!     hold on;
%!     plot(rn, sol.p(idx, j), "-;p(r);r");
%!     plot(rn, sol2.p(idx, j), "-;pref(r);k");
%!     xlabel("r [m]");
%!     ylabel("p [Pa]");
%!     ylim([min(min(sol2.p)), max(max(sol2.p))]);
%!     grid on;
%!     grid minor on;
%!     title(sprintf("radial pressure distribution at Psi=%.1fdeg", Psi(j)*180/pi));
%!   endfor
%!   endif
%!   tol = 1e-4;
%!   assert_simple(sol.p, sol2.p, tol * max(max(abs(sol2.p))));
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
