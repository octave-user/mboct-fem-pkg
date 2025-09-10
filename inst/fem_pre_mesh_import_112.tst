## fem_pre_mesh_import.m:112
%!test
%! try
%! ### TEST 112
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
%!     fputs(fd, "  Surface{10};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"bar1\", 100) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"bar1-convection\", 101) = {10,tmp[0],tmp[2],tmp[6],tmp[7],tmp[8],tmp[9]};\n");
%!     fputs(fd, "Physical Surface(\"bar1-weldseam1\", 102) = {tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"bar1-weldseam2\", 103) = {tmp[5]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tria6h", "tet10h"};
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
%!     fputs(fd, "  Surface{8};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"bar2\", 200) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"bar2-convection\", 201) = {8,tmp[0],tmp[3],tmp[4],tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"bar2-weldseam1\", 204) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"bar2-weldseam2\", 205) = {tmp[6]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tria6h", "tet10h"};
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
%!     fputs(fd, "  Surface{5};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"weld-seam1\", 300) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam1-convection\", 301) = {5,tmp[0],tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam1-weldseam1\", 303) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam1-weldseam2\", 304) = {tmp[4]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tria6h", "tet10h"};
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
%!     fputs(fd, "  Surface{5};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"weld-seam2\", 400) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam2-convection\", 401) = {5,tmp[0],tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam2-weldseam1\", 402) = {tmp[2]};\n");
%!     fputs(fd, "Physical Surface(\"weld-seam2-weldseam2\", 405) = {tmp[4]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_mesh.elem_type = {"tria6h", "tet10h"};
%!   mesh_data(4).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh_data(1).mesh.materials.tet10h = ones(rows(mesh_data(1).mesh.elements.tet10h), 1, "int32");
%!   mesh_data(1).mesh.material_data.E = E;
%!   mesh_data(1).mesh.material_data.nu = nu;
%!   mesh_data(1).mesh.material_data.rho = rho;
%!   mesh_data(1).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(1).mesh.material_data.cp = cp;
%!   mesh_data(1).mesh.material_data.gamma = gamma;
%!   mesh_data(2).mesh.materials.tet10h = ones(rows(mesh_data(2).mesh.elements.tet10h), 1, "int32");
%!   mesh_data(2).mesh.material_data.E = E;
%!   mesh_data(2).mesh.material_data.nu = nu;
%!   mesh_data(2).mesh.material_data.rho = rho;
%!   mesh_data(2).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(2).mesh.material_data.cp = cp;
%!   mesh_data(2).mesh.material_data.gamma = gamma;
%!   mesh_data(3).mesh.materials.tet10h = ones(rows(mesh_data(3).mesh.elements.tet10h), 1, "int32");
%!   mesh_data(3).mesh.material_data.E = E;
%!   mesh_data(3).mesh.material_data.nu = nu;
%!   mesh_data(3).mesh.material_data.rho = rho;
%!   mesh_data(3).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(3).mesh.material_data.cp = cp;
%!   mesh_data(3).mesh.material_data.gamma = gamma;
%!   mesh_data(4).mesh.materials.tet10h = ones(rows(mesh_data(4).mesh.elements.tet10h), 1, "int32");
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
%!   grp_idx_convection = find(mod([mesh.groups.tria6h.id], 100) == 1);
%!   grp_idx_heat_source = find(mod([mesh.groups.tria6h.id], 100) == 2 | ...
%!                              mod([mesh.groups.tria6h.id], 100) == 5);
%!   grp_idx_weldseam2_master = find([mesh.groups.tria6h.id] == 102);
%!   grp_idx_weldseam2_slave = find([mesh.groups.tria6h.id] == 402);
%!   grp_idx_weldseam3_master = find([mesh.groups.tria6h.id] == 103);
%!   grp_idx_weldseam3_slave = find([mesh.groups.tria6h.id] == 303);
%!   grp_idx_weldseam4_master = find([mesh.groups.tria6h.id] == 204);
%!   grp_idx_weldseam4_slave = find([mesh.groups.tria6h.id] == 304);
%!   grp_idx_weldseam5_master = find([mesh.groups.tria6h.id] == 205);
%!   grp_idx_weldseam5_slave = find([mesh.groups.tria6h.id] == 405);
%!   mesh.elements.convection.tria6h.nodes = mesh.elements.tria6h([[mesh.groups.tria6h(grp_idx_convection)].elements], :);
%!   mesh.elements.convection.tria6h.h = repmat(he, size(mesh.elements.convection.tria6h.nodes));
%!   mesh.elements.sfncon6h(1).master = mesh.elements.tria6h([[mesh.groups.tria6h(grp_idx_weldseam2_master)].elements], :);
%!   mesh.elements.sfncon6h(1).slave = [[mesh.groups.tria6h(grp_idx_weldseam2_slave)].nodes](:);
%!   mesh.elements.sfncon6h(1).maxdist = 1e-6;
%!   mesh.elements.sfncon6h(2).master = mesh.elements.tria6h([[mesh.groups.tria6h(grp_idx_weldseam3_master)].elements], :);
%!   mesh.elements.sfncon6h(2).slave = [[mesh.groups.tria6h(grp_idx_weldseam3_slave)].nodes](:);
%!   mesh.elements.sfncon6h(2).maxdist = 1e-6;
%!   mesh.elements.sfncon6h(3).master = mesh.elements.tria6h([[mesh.groups.tria6h(grp_idx_weldseam4_master)].elements], :);
%!   mesh.elements.sfncon6h(3).slave = [[mesh.groups.tria6h(grp_idx_weldseam4_slave)].nodes](:);
%!   mesh.elements.sfncon6h(3).maxdist = 1e-6;
%!   mesh.elements.sfncon6h(4).master = mesh.elements.tria6h([[mesh.groups.tria6h(grp_idx_weldseam5_master)].elements], :);
%!   mesh.elements.sfncon6h(4).slave = [[mesh.groups.tria6h(grp_idx_weldseam5_slave)].nodes](:);
%!   mesh.elements.sfncon6h(4).maxdist = 1e-6;
%!   load_case = struct("heat_source", cell(1, 2), "convection", cell(1, 2));
%!   load_case(1).heat_source.tria6h.nodes = mesh.elements.tria6h([[mesh.groups.tria6h(grp_idx_heat_source)].elements], :);
%!   load_case(1).heat_source.tria6h.q = repmat(qs, size(load_case(1).heat_source.tria6h.nodes));
%!   load_case(2).heat_source = struct();
%!   for i=1:numel(load_case)
%!     load_case(i).convection.tria6h.theta = repmat(thetae, size(mesh.elements.convection.tria6h.nodes));
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
%!   node_id_source = unique([[mesh.groups.tria6h(grp_idx_heat_source)].nodes]);
%!   theta_mean = mean(sol.theta(dof_map.ndof(node_id_source), :), 1);
%!   theta_max = max(sol.theta(dof_map.ndof(node_id_source), :), [], 1);
%!   theta_min = min(sol.theta(dof_map.ndof(node_id_source), :), [], 1);
%!   theta_ref = 900;
%!   idx_theta_ref = find(theta_mean(2:end) < theta_ref & theta_mean(1:end - 1) >= theta_ref)(1);
%!   mesh2 = mesh;
%!   load_case2.domain = FEM_DO_STRUCTURAL;
%!   constr = [FEM_CT_SLIDING, FEM_CT_FIXED, FEM_CT_FIXED, FEM_CT_SLIDING];
%!   for j=1:numel(mesh.elements.sfncon6h)
%!     mesh2.elements.sfncon6h(j).constraint = constr(j);
%!   endfor
%!   load_case2.locked_dof = false(size(mesh2.nodes));
%!   load_case2.dTheta = sol.theta(:, idx_theta_ref);
%!   dof_map2 = fem_ass_dof_map(mesh2, load_case2);
%!   [mat_ass2.K, ...
%!    mat_ass2.R, ...
%!    mat_ass2.mat_info, ...
%!    mat_ass2.mesh_info] = fem_ass_matrix(mesh2, ...
%!                                         dof_map2, ...
%!                                         [FEM_MAT_STIFFNESS, ...
%!                                          FEM_VEC_LOAD_CONSISTENT], ...
%!                                         load_case2);
%!   mat_ass2.K += speye(columns(mat_ass2.K)) * (sqrt(eps) * max(max(abs(mat_ass2.K))));
%!   sol2 = fem_sol_static(mesh2, dof_map2, mat_ass2);
%!   [sol2.stress, ...
%!    sol2.strain] = fem_ass_matrix(mesh2, ...
%!                                  dof_map2, ...
%!                                  [FEM_VEC_STRESS_CAUCH, ...
%!                                   FEM_VEC_STRAIN_TOTAL], ...
%!                                  load_case2, sol2);
%!   load_case3 = rmfield(load_case2, "dTheta");
%!   mesh3 = mesh;
%!   for i=1:numel(mesh3.elements.sfncon6h)
%!     mesh3.elements.sfncon6h(i).constraint = FEM_CT_FIXED;
%!   endfor
%!   mesh3dummy.elements.sfncon6 = mesh3.elements.sfncon6h;
%!   for i=1:numel(mesh3dummy.elements.sfncon6)
%!     mesh3dummy.elements.sfncon6(i).constraint = FEM_CT_SLIDING;
%!   endfor
%!   mesh3.elements.joints = fem_pre_mesh_constr_surf_to_node(mesh3.nodes, mesh3.elements, load_case3.domain).joints;
%!   mesh3dummy.elements.joints = fem_pre_mesh_constr_surf_to_node(mesh3.nodes, mesh3dummy.elements, load_case3.domain).joints;
%!   mesh3.elements = rmfield(mesh3.elements, "sfncon6h");
%!   load_case3.joints = struct("U", mat2cell(zeros(3, numel(mesh3.elements.joints)), ...
%!                                            3, ones(1, numel(mesh3.elements.joints))));
%!   dof_map3 = fem_ass_dof_map(mesh3, load_case3);
%!   for i=1:numel(mesh3.elements.joints)
%!     dU = zeros(numel(mesh3.elements.joints(i).nodes) * 6, 1);
%!     for j=1:numel(mesh3.elements.joints(i).nodes)
%!       for k=1:3
%!         dU((j - 1) * 6 + k) = sol2.def(mesh3.elements.joints(i).nodes(j), k);
%!       endfor
%!     endfor
%!     load_case3.joints(i).U = mesh3.elements.joints(i).C * dU;
%!   endfor
%!   [mat_ass3.K, ...
%!    mat_ass3.R, ...
%!    mat_ass3.mat_info, ...
%!    mat_ass3.mesh_info] = fem_ass_matrix(mesh3, ...
%!                                         dof_map3, ...
%!                                         [FEM_MAT_STIFFNESS, ...
%!                                          FEM_VEC_LOAD_CONSISTENT], ...
%!                                         load_case3);
%!   mat_ass3.K += speye(columns(mat_ass3.K)) * (sqrt(eps) * max(max(abs(mat_ass3.K))));
%!   sol3 = fem_sol_static(mesh3, dof_map3, mat_ass3);
%!   [sol3.stress, ...
%!    sol3.strain] = fem_ass_matrix(mesh3, ...
%!                                  dof_map3, ...
%!                                  [FEM_VEC_STRESS_CAUCH, ...
%!                                   FEM_VEC_STRAIN_TOTAL], ...
%!                                  load_case3, sol3);
%!   opt_post.elem_types = {"tet10h", "tria6h"};
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
