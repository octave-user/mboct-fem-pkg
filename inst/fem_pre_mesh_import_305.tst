## fem_pre_mesh_import.m:305
%!test
%! try
%! ### TEST 305
%! ###
%! ### Munjal, M.L. (1987). Acoustics of ducts and mufflers with application to exhaust and ventilation system design.
%! ### New York: Wiley. pp. 58-60. ISBN 0471847380.
%! ###
%! ### Leray, David and Goth, Yvon, "Acoustic Calculation With the Free Solver Code_Aster" (2012).
%! ### International Compressor Engineering Conference. Paper 2133.
%! ### https://docs.lib.purdue.edu/icec/2133
%! ###
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
%!   l1 = 10e-3;
%!   l2 = 50e-3;
%!   l3 = 10e-3;
%!   d1 = 50e-3 / 3;
%!   d2 = 50e-3;
%!   d3 = d1;
%!   c = 520;
%!   rho = 1;
%!   fc = 1.84 * c / (pi * d2);
%!   omegac = 2 * pi * fc;
%!   lambda = c / fc;
%!   dx = lambda / 6;
%!   vx0 = 1;
%!   h = d2^2 / d1^2;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l1=%g;\n", l1);
%!   fprintf(fd, "l2=%g;\n", l2);
%!   fprintf(fd, "l3=%g;\n", l3);
%!   fprintf(fd, "d1=%g;\n", d1);
%!   fprintf(fd, "d2=%g;\n", d2);
%!   fprintf(fd, "d3=%g;\n", d3);
%!   fprintf(fd, "dx = %g;\n", dx);
%!   fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!   fputs(fd, "Point(2) = {l1+l2+l3,0,0,dx};\n");
%!   fputs(fd, "Point(3) = {l1+l2+l3,0,0.5 * d3,dx};\n");
%!   fputs(fd, "Point(4) = {l1+l2,0,0.5*d3,dx};\n");
%!   fputs(fd, "Point(5) = {l1+l2,0,0.5*d2,dx};\n");
%!   fputs(fd, "Point(6) = {l1,0,0.5*d2,dx};\n");
%!   fputs(fd, "Point(7) = {l1,0,0.5*d1,dx};\n");
%!   fputs(fd, "Point(8) = {0,0,0.5*d1,dx};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Line(7) = {7,8};\n");
%!   fputs(fd, "Line(8) = {8,1};\n");
%!   fputs(fd, "Line Loop(9) = {1,2,3,4,5,6,7,8};\n");
%!   fputs(fd, "Plane Surface(10) = {9};\n");
%!   fputs(fd, "tmp[] = Extrude {{1,0,0}, {0,0,0}, Pi/2}{\n");
%!   fputs(fd, "  Surface{10};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Physical Volume(\"volume\",3) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"inlet\",1) = {tmp[8]};\n");
%!   fputs(fd, "Physical Surface(\"outlet\",2) = {tmp[2]};\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
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
%!   grp_idx_inlet = find([mesh.groups.tria10.id] == 1);
%!   grp_idx_outlet = find([mesh.groups.tria10.id] == 2);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   mesh.elements.acoustic_boundary.tria10 = mesh.elements.tria10([[mesh.groups.tria10([grp_idx_inlet])].elements], :);
%!   mesh.materials.acoustic_boundary.tria10 = ones(rows(mesh.elements.acoustic_boundary.tria10), 1, "int32");
%!   mesh.elements.particle_velocity.tria10.nodes = mesh.elements.tria10([[mesh.groups.tria10([grp_idx_inlet])].elements], :);
%!   mesh.materials.particle_velocity.tria10 = ones(rows(mesh.elements.particle_velocity.tria10.nodes), 1, "int32");
%!   mesh.elements.acoustic_impedance.tria10.nodes = mesh.elements.tria10([[mesh.groups.tria10(grp_idx_outlet)].elements], :);
%!   mesh.elements.acoustic_impedance.tria10.z = repmat(rho * c, size(mesh.elements.acoustic_impedance.tria10.nodes));
%!   mesh.materials.acoustic_impedance.tria10 = ones(rows(mesh.elements.acoustic_impedance.tria10.nodes), 1, "int32");
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   load_case(1).particle_velocity.tria10.vn = repmat(-real(vx0), numel(mesh.groups.tria10(grp_idx_inlet).elements), 10);
%!   load_case(2).particle_velocity.tria10.vn = repmat(-imag(vx0), numel(mesh.groups.tria10(grp_idx_inlet).elements), 10);
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
%!   omegai = linspace(sqrt(eps) * omegac, omegac, 50);
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   node_idx_inlet = [[mesh.groups.tria10(grp_idx_inlet)].nodes];
%!   node_idx_outlet = [[mesh.groups.tria10(grp_idx_outlet)].nodes];
%!   err = R = zeros(size(omegai));
%!   for i=1:numel(omegai)
%!     Keff = complex(-omegai(i)^2 * mat_ass.Ma + 1j * omegai(i) * mat_ass.Da + mat_ass.Ka);
%!     sol.Phi = complex(fem_sol_factor(Keff, opt_sol) \ Reff);
%!     sol.PhiP = complex(1j * omegai(i) * sol.Phi);
%!     sol.p = -sol.PhiP;
%!     err(i) = norm(Keff * sol.Phi - Reff) / norm(Keff * sol.Phi + Reff);
%!     fprintf("%d: %.1fHz %.2e\n", i, omegai(i)/(2*pi), err(i));
%!     [sol.particle_velocity, ...
%!      sol.acoustic_intensity] = fem_ass_matrix(mesh, ...
%!                                               dof_map, ...
%!                                               [FEM_VEC_PARTICLE_VELOCITY_C, ...
%!                                                FEM_SCA_ACOUSTIC_INTENSITY_C], ...
%!                                               load_case, ...
%!                                               sol);
%!     pin = sol.p(dof_map.ndof(mesh.elements.tria10(mesh.groups.tria10(grp_idx_inlet).elements, :)));
%!     vxin = -sol.particle_velocity.vn.tria10(1:numel(mesh.groups.tria10(grp_idx_inlet).elements), :);
%!     z = rho * c;
%!     R(i) = abs(mean(mean((pin - z * vxin) ./ (pin + z * vxin))));
%!   endfor
%!   T = 1 ./ (1 - R.^2);
%!   TL = 10 * log10(T);
%!   TLref = 10 * log10(1 + 1/4 * (h - 1 / h)^2 * sin(omegai / c * l2).^2);
%!   if (do_plot)
%!     figure("visible","off");
%!     hold on;
%!     plot(omegai/(2*pi), TL, "-;FEM solution;r");
%!     plot(omegai/(2*pi), TLref, "-;reference solution;k");
%!     xlabel("f [Hz]");
%!     ylabel("TL [dB]");
%!     grid on;
%!     grid minor on;
%!     title("transmission loss of duct");
%!   endif
%!   tol = 1.5;
%!   assert_simple(TL, TLref, tol);
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
