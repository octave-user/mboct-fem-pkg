## fem_pre_mesh_import.m: 475
%!test
%! try
%! ### TEST 475
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   dx = 2.5e-3;
%!   L1 = 20e-3;
%!   L2 = 15e-3;
%!   t1 = 10e-3;
%!   t2 = 10e-3;
%!   w1 = 50e-3;
%!   h2 = 50e-3;
%!   a1 = 8e-3;
%!   a2 = 8e-3;
%!   E = 210000e6;
%!   nu = 0.3;
%!   rho = 7850;
%!   lambda = 45;
%!   cp = 478;
%!   gamma = 12.5e-6;
%!   he = 2000;
%!   thetae = 22;
%!   theta0 = 22;
%!   qs = 10000000e3;
%!   ts = 1e-3;
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
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "dx=%g;\n", dx);
%!     fprintf(fd, "L1=%g;\n", L1);
%!     fprintf(fd, "L2=%g;\n", L2);
%!     fprintf(fd, "t1=%g;\n", t1);
%!     fprintf(fd, "t2=%g;\n", t2);
%!     fprintf(fd, "w1=%g;\n", w1);
%!     fprintf(fd, "h2=%g;\n", h2);
%!     fprintf(fd, "a1=%g;\n", a1);
%!     fprintf(fd, "a2=%g;\n", a2);
%!     fputs(fd, "Point(1) = {-w1/2, t1, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(2) = {-t2/2 - a1, t1, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(3) = {-t2/2, t1, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(4) = {t2/2, t1, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(5) = {t2/2 + a1, t1, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(6) = {w1/2, t1, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(7) = {w1/2, 0, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Point(8) = {-w1/2, 0, -(L1 - L2) / 2, dx};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Curve Loop(9) = {1,2,3,4,5,6,7,8};\n");
%!     fputs(fd, "Plane Surface(10) = {9};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,L1} {\n");
%!     fputs(fd, "  Surface{10}; Layers{Ceil(L1/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{10, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"bar1\", 100) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"bar1-convection\", 101) = {10,tmp[0],tmp[2],tmp[6],tmp[7],tmp[8],tmp[9]};\n");
%!     fputs(fd, "Physical Surface(\"bar1-weldseam1\", 102) = {tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"bar1-weldseam2\", 103) = {tmp[5]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"penta18", "iso27", "quad9"};
%!   mesh_data = struct("mesh", cell(1, 4));
%!   mesh_data(1).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "dx=%g;\n", dx);
%!     fprintf(fd, "L1=%g;\n", L1);
%!     fprintf(fd, "L2=%g;\n", L2);
%!     fprintf(fd, "t1=%g;\n", t1);
%!     fprintf(fd, "t2=%g;\n", t2);
%!     fprintf(fd, "w1=%g;\n", w1);
%!     fprintf(fd, "h2=%g;\n", h2);
%!     fprintf(fd, "a1=%g;\n", a1);
%!     fprintf(fd, "a2=%g;\n", a2);
%!     fputs(fd, "Point(1) = {t2/2, t1, 0, dx};\n");
%!     fputs(fd, "Point(2) = {t2/2, t2 + a2, 0, dx};\n");
%!     fputs(fd, "Point(3) = {t2/2, t1 + h2, 0, dx};\n");
%!     fputs(fd, "Point(4) = {-t2/2, t1 + h2, 0, dx};\n");
%!     fputs(fd, "Point(5) = {-t2/2, t1 + a2, 0, dx};\n");
%!     fputs(fd, "Point(6) = {-t2/2, t1, 0, dx};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 1};\n");
%!     fputs(fd, "Curve Loop(7) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(8) = {7};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,L2} {\n");
%!     fputs(fd, "  Surface{8}; Layers{Ceil(L2/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{8, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"bar2\", 200) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"bar2-convection\", 201) = {8,tmp[0],tmp[3],tmp[4],tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"bar2-weldseam1\", 204) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"bar2-weldseam2\", 205) = {tmp[6]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"penta18", "iso27", "quad9"};
%!   mesh_data(2).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "dx=%g;\n", dx);
%!     fprintf(fd, "L1=%g;\n", L1);
%!     fprintf(fd, "L2=%g;\n", L2);
%!     fprintf(fd, "t1=%g;\n", t1);
%!     fprintf(fd, "t2=%g;\n", t2);
%!     fprintf(fd, "w1=%g;\n", w1);
%!     fprintf(fd, "h2=%g;\n", h2);
%!     fprintf(fd, "a1=%g;\n", a1);
%!     fprintf(fd, "a2=%g;\n", a2);
%!     fputs(fd, "Point(1) = {t2/2, t1, 0, dx};\n");
%!     fputs(fd, "Point(2) = {t2/2 + a1, t1, 0, dx};\n");
%!     fputs(fd, "Point(3) = {t2/2, t1 + a2, 0, dx};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 1};\n");
%!     fputs(fd, "Curve Loop(4) = {1,2,3};\n");
%!     fputs(fd, "Plane Surface(5) = {4};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,L2} {\n");
%!     fputs(fd, "  Surface{5}; Layers{Ceil(L2/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{5, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"weld-seam1\", 300) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam1-convection\", 301) = {5,tmp[0],tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam1-weldseam1\", 303) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam1-weldseam2\", 304) = {tmp[4]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"penta18", "iso27", "quad9"};
%!   mesh_data(3).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "dx=%g;\n", dx);
%!     fprintf(fd, "L1=%g;\n", L1);
%!     fprintf(fd, "L2=%g;\n", L2);
%!     fprintf(fd, "t1=%g;\n", t1);
%!     fprintf(fd, "t2=%g;\n", t2);
%!     fprintf(fd, "w1=%g;\n", w1);
%!     fprintf(fd, "h2=%g;\n", h2);
%!     fprintf(fd, "a1=%g;\n", a1);
%!     fprintf(fd, "a2=%g;\n", a2);
%!     fputs(fd, "Point(1) = {-t2/2, t1, 0, dx};\n");
%!     fputs(fd, "Point(2) = {-t2/2 - a1, t1, 0, dx};\n");
%!     fputs(fd, "Point(3) = {-t2/2, t1 + a2, 0, dx};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 1};\n");
%!     fputs(fd, "Curve Loop(4) = {1,2,3};\n");
%!     fputs(fd, "Plane Surface(5) = {4};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,L2} {\n");
%!     fputs(fd, "  Surface{5}; Layers{Ceil(L2/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{5, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"weld-seam2\", 400) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam2-convection\", 401) = {5,tmp[0],tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam2-weldseam1\", 402) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam2-weldseam2\", 405) = {tmp[4]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"penta18", "iso27", "quad9"};
%!   mesh_data(4).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh_data(1).mesh.materials.iso27 = ones(rows(mesh_data(1).mesh.elements.iso27), 1, "int32");
%!   if (isfield(mesh_data(1).mesh.elements, "penta18"))
%!     mesh_data(1).mesh.materials.penta18 = ones(rows(mesh_data(1).mesh.elements.penta18), 1, "int32");
%!   endif
%!   mesh_data(1).mesh.material_data.E = E;
%!   mesh_data(1).mesh.material_data.nu = nu;
%!   mesh_data(1).mesh.material_data.rho = rho;
%!   mesh_data(1).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(1).mesh.material_data.cp = cp;
%!   mesh_data(1).mesh.material_data.gamma = gamma;
%!   mesh_data(2).mesh.materials.iso27 = ones(rows(mesh_data(2).mesh.elements.iso27), 1, "int32");
%!   if (isfield(mesh_data(2).mesh.elements, "penta18"))
%!     mesh_data(2).mesh.materials.penta18 = ones(rows(mesh_data(2).mesh.elements.penta18), 1, "int32");
%!   endif
%!   mesh_data(2).mesh.material_data.E = E;
%!   mesh_data(2).mesh.material_data.nu = nu;
%!   mesh_data(2).mesh.material_data.rho = rho;
%!   mesh_data(2).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(2).mesh.material_data.cp = cp;
%!   mesh_data(2).mesh.material_data.gamma = gamma;
%!   mesh_data(3).mesh.materials.iso27 = ones(rows(mesh_data(3).mesh.elements.iso27), 1, "int32");
%!   if (isfield(mesh_data(3).mesh.elements, "penta18"))
%!     mesh_data(3).mesh.materials.penta18 = ones(rows(mesh_data(3).mesh.elements.penta18), 1, "int32");
%!   endif
%!   mesh_data(3).mesh.material_data.E = E;
%!   mesh_data(3).mesh.material_data.nu = nu;
%!   mesh_data(3).mesh.material_data.rho = rho;
%!   mesh_data(3).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(3).mesh.material_data.cp = cp;
%!   mesh_data(3).mesh.material_data.gamma = gamma;
%!   mesh_data(4).mesh.materials.iso27 = ones(rows(mesh_data(4).mesh.elements.iso27), 1, "int32");
%!   if (isfield(mesh_data(4).mesh.elements, "penta18"))
%!     mesh_data(4).mesh.materials.penta18 = ones(rows(mesh_data(4).mesh.elements.penta18), 1, "int32");
%!   endif
%!   mesh_data(4).mesh.material_data.E = E;
%!   mesh_data(4).mesh.material_data.nu = nu;
%!   mesh_data(4).mesh.material_data.rho = rho;
%!   mesh_data(4).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(4).mesh.material_data.cp = cp;
%!   mesh_data(4).mesh.material_data.gamma = gamma;
%!   opt_merge.group_id = "preserve";
%!   mesh = fem_post_mesh_merge(mesh_data, opt_merge);
%!   dof_stat.locked_dof = false(rows(mesh.nodes), 1);
%!   dof_stat.domain = FEM_DO_THERMAL;
%!   grp_idx_convection_tria6 = find(mod([mesh.groups.quad9.id], 100) == 1);
%!   grp_idx_convection_quad9 = find(mod([mesh.groups.quad9.id], 100) == 1);
%!   grp_idx_heat_source = find(mod([mesh.groups.quad9.id], 100) == 2 | ...
%!                              mod([mesh.groups.quad9.id], 100) == 5);
%!   grp_idx_weldseam2_master = find([mesh.groups.quad9.id] == 102);
%!   grp_idx_weldseam2_slave = find([mesh.groups.quad9.id] == 402);
%!   grp_idx_weldseam3_master = find([mesh.groups.quad9.id] == 103);
%!   grp_idx_weldseam3_slave = find([mesh.groups.quad9.id] == 303);
%!   grp_idx_weldseam4_master = find([mesh.groups.quad9.id] == 204);
%!   grp_idx_weldseam4_slave = find([mesh.groups.quad9.id] == 304);
%!   grp_idx_weldseam5_master = find([mesh.groups.quad9.id] == 205);
%!   grp_idx_weldseam5_slave = find([mesh.groups.quad9.id] == 405);
%!   mesh.elements.convection.quad9.nodes = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_convection_tria6)].elements], :);
%!   mesh.elements.convection.quad9.h = repmat(he, size(mesh.elements.convection.quad9.nodes));
%!   mesh.elements.convection.quad9.nodes = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_convection_tria6)].elements], :);
%!   mesh.elements.convection.quad9.h = repmat(he, size(mesh.elements.convection.quad9.nodes));
%!   mesh.elements.sfncon9(1).master = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_weldseam2_master)].elements], :);
%!   mesh.elements.sfncon9(1).slave = [[mesh.groups.quad9(grp_idx_weldseam2_slave)].nodes](:);
%!   mesh.elements.sfncon9(1).maxdist = 1e-6;
%!   mesh.elements.sfncon9(2).master = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_weldseam3_master)].elements], :);
%!   mesh.elements.sfncon9(2).slave = [[mesh.groups.quad9(grp_idx_weldseam3_slave)].nodes](:);
%!   mesh.elements.sfncon9(2).maxdist = 1e-6;
%!   mesh.elements.sfncon9(3).master = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_weldseam4_master)].elements], :);
%!   mesh.elements.sfncon9(3).slave = [[mesh.groups.quad9(grp_idx_weldseam4_slave)].nodes](:);
%!   mesh.elements.sfncon9(3).maxdist = 1e-6;
%!   mesh.elements.sfncon9(4).master = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_weldseam5_master)].elements], :);
%!   mesh.elements.sfncon9(4).slave = [[mesh.groups.quad9(grp_idx_weldseam5_slave)].nodes](:);
%!   mesh.elements.sfncon9(4).maxdist = 1e-6;
%!   load_case = struct("heat_source", cell(1, 2), "convection", cell(1, 2));
%!   load_case(1).heat_source.quad9.nodes = mesh.elements.quad9([[mesh.groups.quad9(grp_idx_heat_source)].elements], :);
%!   load_case(1).heat_source.quad9.q = repmat(qs, size(load_case(1).heat_source.quad9.nodes));
%!   load_case(2).heat_source = struct();
%!   for i=1:numel(load_case)
%!     load_case(i).convection.quad9.theta = repmat(thetae, size(mesh.elements.convection.quad9.nodes));
%!     load_case(i).convection.quad9.theta = repmat(thetae, size(mesh.elements.convection.quad9.nodes));
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, dof_stat);
%!   [mat_ass.C, ...
%!    mat_ass.Kk, ...
%!    mat_ass.Qc, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_HEAT_CAPACITY, ...
%!                                         FEM_MAT_THERMAL_COND, ...
%!                                         FEM_VEC_LOAD_THERMAL], ...
%!                                        load_case);
%!   dt = rho * cp * dx^2 / lambda;
%!   T = 5;
%!   alpha = 0.5;
%!   sol.t = 0:dt:T;
%!   U = zeros(columns(mat_ass.Kk), numel(sol.t));
%!   U(dof_map.idx_node, 1) = theta0;
%!   opts.number_of_threads = mbdyn_solver_num_threads_default();
%!   opts.solver = "pardiso";
%!   A = (1 / dt) * mat_ass.C + alpha * mat_ass.Kk;
%!   Afact = fem_sol_factor(A, opts);
%!   ti = [0, (2 - sqrt(eps)) * dt, 2 * dt, 3 * dt, (3 + sqrt(eps)) * dt, T];
%!   fi = [0, 0, 1, 1, 0, 0];
%!   f = interp1(ti, fi * ts / dt, sol.t, "linear");
%!   for i=2:numel(sol.t)
%!     Qci = mat_ass.Qc(:, 2) * (alpha * (1 - f(i)) + (1 - alpha) * (1 - f(i - 1))) + ...
%!           mat_ass.Qc(:, 1) * (alpha * f(i) + (1 - alpha) * f(i - 1));
%!     U(:, i) = Afact \ (mat_ass.C * (U(:, i - 1)) / dt - mat_ass.Kk * (U(:, i - 1) * (1 - alpha)) + Qci);
%!   endfor
%!   sol.theta = U(dof_map.idx_node, :);
%!   opt_post.elem_types = {"iso27", "quad9", "quad9"};
%!   opt_post.skin_only = false;
%!   if (do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t, min(sol.theta, [], 1), "-;min(theta);r");
%!     plot(sol.t, max(sol.theta, [], 1), "-;max(theta);b");
%!     xlabel("t [s]");
%!     ylabel("theta [degC]");
%!     grid on;
%!     grid minor on;
%!     title("temperature versus time");
%!   endif
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
