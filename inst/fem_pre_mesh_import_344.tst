## fem_pre_mesh_import.m:344
%!test
%! ### TEST 344
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
%!     vinput = -1;
%!     poutput = 0;
%!     T = 2 * l / c;
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
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
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
%!   grp_idx_input = find([mesh.groups.quad9.id] == 1);
%!   grp_idx_output = find([mesh.groups.quad9.id] == 2);
%!   node_idx_input = mesh.groups.quad9(grp_idx_input).nodes;
%!   node_idx_output = mesh.groups.quad9(grp_idx_output).nodes;
%!   p_constr = repmat(poutput, 1, numel(node_idx_output));
%!   mesh.elements.acoustic_constr = struct("C", mat2cell(ones(1, numel(node_idx_output)), 1, ones(1, numel(node_idx_output))), ...
%!                                          "nodes", mat2cell(node_idx_output, 1, ones(1, numel(node_idx_output))), ...
%!                                          "scale", mat2cell(repmat(dt, 1, numel(node_idx_output)), 1, ones(1, numel(node_idx_output))));
%!   load_case.acoustic_constr = struct("p", mat2cell(p_constr, 1, ones(1, numel(p_constr))));
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_ACOUSTICS;
%!   mesh.elements.particle_velocity.quad9.nodes = mesh.elements.quad9(mesh.groups.quad9(grp_idx_input).elements, :);
%!   mesh.materials.particle_velocity.quad9 = ones(rows(mesh.elements.particle_velocity.quad9.nodes), 1, "int32");
%!   load_case.particle_velocity.quad9.vn = repmat(vinput, numel(mesh.groups.quad9(grp_idx_input).elements), 9);
%!   mesh.materials.iso27 = ones(rows(mesh.elements.iso27), 1, "int32");
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
%!   opt_sol.delta = 0.6;
%!   sol_dat = fem_sol_transient_init(mat_ass.Ma, mat_ass.Da, mat_ass.Ka, dt, opt_sol);
%!   [x, idx] = sort(mesh.nodes(:, 1));
%!   for i=2:numel(t)
%!     fprintf(stderr, "Step %d/%d\n", i, numel(t));
%!     [Phi(:, i), dPhi_dt(:, i), dPhi_dt2(:, i)] = fem_sol_transient_step(Phi(:, i - 1), dPhi_dt(:, i - 1), dPhi_dt2(:, i - 1), mat_ass.Ra, sol_dat);
%!   endfor
%!   sol.t = t;
%!   sol.p = -dPhi_dt(dof_map.ndof, :);
%!   for i=1:20:numel(t)
%!     figure("visible", "off");
%!     hold on;
%!     plot(x, sol.p(idx, i), "-;p(x);r");
%!     ylim([min(min(sol.p)), max(max(sol.p))]);
%!     xlabel("x [m]");
%!     ylabel("p [Pa]");
%!     grid on;
%!     grid minor on;
%!     title(sprintf("t=%.2fs", t(i)));
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
