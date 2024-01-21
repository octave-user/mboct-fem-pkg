## fem_pre_mesh_import.m:143
%!test
%! ### TEST 143
%! ## The 1-D Wave Equation
%! ## 18.303 Linear Partial Differential Equations
%! ## Matthew J. Hancock
%! ## Fall 2006
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
%!     l = 1.5;
%!     dx = 3e-3;
%!     w = dx;
%!     h = dx;
%!     c = 340;
%!     rho = 1.25;
%!     pinput = 0;
%!     poutput = 0;
%!     T = l / c;
%!     dt = dx / c;
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
%!     fputs(fd, "Recombine Surface {6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"input\",1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"output\",2) = {tmp[0]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
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
%!   grp_idx_input = find([mesh.groups.iso4.id] == 1);
%!   grp_idx_output = find([mesh.groups.iso4.id] == 2);
%!   node_idx_input = mesh.groups.iso4(grp_idx_input).nodes;
%!   node_idx_output = mesh.groups.iso4(grp_idx_output).nodes;
%!   node_idx_constr = [node_idx_input, node_idx_output];
%!   p_constr = [repmat(pinput, 1, numel(node_idx_input)), ...
%!               repmat(poutput, 1, numel(node_idx_output))];
%!   mesh.elements.acoustic_constr = struct("C", mat2cell(ones(1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))), ...
%!                                   "nodes", mat2cell(node_idx_constr, 1, ones(1, numel(node_idx_constr))), ...
%!                                   "scale", mat2cell(repmat(dt, 1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))));
%!   load_case.acoustic_constr = struct("p", mat2cell(p_constr, 1, ones(1, numel(p_constr))));
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.Da, ...
%!    mat_ass.Ra] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                  FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                  FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                  FEM_VEC_LOAD_ACOUSTICS], ...
%!                                 load_case);
%!   t = 0:dt:T;
%!   Phi = dPhi_dt = dPhi_dt2 = zeros(dof_map.totdof, numel(t));
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.delta = 0.5;
%!   sol_dat = fem_sol_transient_init(mat_ass.Ma, mat_ass.Da, mat_ass.Ka, dt, opt_sol);
%!   [x, idx] = sort(mesh.nodes(:, 1));
%!   f = @(x) (cos(pi * (2 * x / l - 1)) + 1) / 2;
%!   g = @(x) c * sin(5 * pi * x / l);
%!   N = 200;
%!   alpha = zeros(1, N);
%!   beta = zeros(1, N);
%!   for n=1:N
%!     alpha(n) = 2 * quad(@(x) f(x) / l .* sin(n * pi * x / l), 0, l);
%!     beta(n) = 2 / (n * pi) * quad(@(x) g(x) * T / l .* sin(n * pi * x / l), 0, l);
%!   endfor
%!   Phi(dof_map.ndof, 1) = f(mesh.nodes(:, 1));
%!   dPhi_dt(dof_map.ndof, 1) = g(mesh.nodes(:, 1));
%!   Phiref = zeros(numel(x), numel(t));
%!   for n=1:N
%!     Phiref += (alpha(n) * cos(n * pi * t / T) + beta(n) * sin(n * pi * t / T)) .* sin(n * pi * x / l);
%!   endfor
%!   for i=2:numel(t)
%!     fprintf(stderr, "Step %d/%d\n", i, numel(t));
%!     [Phi(:, i), dPhi_dt(:, i), dPhi_dt2(:, i)] = fem_sol_transient_step(Phi(:, i - 1), dPhi_dt(:, i - 1), dPhi_dt2(:, i - 1), mat_ass.Ra, sol_dat);
%!   endfor
%!   sol.t = t;
%!   sol.p = -dPhi_dt(dof_map.ndof, :);
%!   for i=1:20:numel(t)
%!     figure("visible", "off");
%!     hold on;
%!     plot(x, Phi(idx, i), "-;Phi(x);1");
%!     plot(x, Phiref(:, i), "-;Phiref(x);0");
%!     ylim([min(min(Phi)), max(max(Phi))]);
%!     xlabel("x [m]");
%!     ylabel("Phi [m^2/s]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("t=%.2fs", t(i)));
%!   endfor
%!   tol = 5e-3;
%!   assert_simple(Phi(idx, :), Phiref, tol * max(max(abs(Phiref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
