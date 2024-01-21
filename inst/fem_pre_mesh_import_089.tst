## fem_pre_mesh_import.m:89
%!test
%! ### TEST 89
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
%!   a = 10e-3;
%!   b = 10e-3;
%!   c = 10e-3;
%!   d = 0.25e-3;
%!   dx = 0.25e-3;
%!   lambda = 50;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "d=%g;\n", d);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {-b,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {0,0,0,dx};\n");
%!     fputs(fd, "Point(3) = {0,c,0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,1};\n");
%!     fputs(fd, "Line Loop(4) = {1,2,3};\n");
%!     fputs(fd, "Plane Surface(5) = {4};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,d} {\n");
%!     fputs(fd, "  Surface{5};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",5) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"theta\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"master\",2) = {tmp[3]};\n");
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
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh_data(1).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   mesh_data(1).mesh.materials.tet10h = ones(rows(mesh_data(1).mesh.elements.tet10h), 1, "int32");
%!   mesh_data(1).mesh.material_data.E = 210000e6;
%!   mesh_data(1).mesh.material_data.nu = 0.3;
%!   mesh_data(1).mesh.material_data.rho = 7850;
%!   mesh_data(1).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(1).mesh.material_data.cp = 465;
%!   [~] = unlink([filename, ".msh"]);
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "d=%g;\n", d);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {b,0,0,dx};\n");
%!     fputs(fd, "Point(3) = {0,c,0,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,1};\n");
%!     fputs(fd, "Line Loop(4) = {1,2,3};\n");
%!     fputs(fd, "Plane Surface(5) = {4};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0,d} {\n");
%!     fputs(fd, "  Surface{5};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",5) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"theta\",1) = {tmp[3]};\n");
%!     fputs(fd, "Physical Surface(\"slave\",3) = {tmp[4]};\n");
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
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh_data(2).mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   mesh_data(2).mesh.materials.tet10h = ones(rows(mesh_data(2).mesh.elements.tet10h), 1, "int32");
%!   mesh_data(2).mesh.material_data.E = 210000e6;
%!   mesh_data(2).mesh.material_data.nu = 0.3;
%!   mesh_data(2).mesh.material_data.rho = 7850;
%!   mesh_data(2).mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh_data(2).mesh.material_data.cp = 465;
%!   [~] = unlink([filename, ".msh"]);
%!   opt_merge.group_id = "preserve";
%!   mesh = fem_post_mesh_merge(mesh_data, opt_merge);
%!   thetae = [100, 200];
%!   x = mesh.nodes(:, 1);
%!   theta_s = thetae(1) + (thetae(2) - thetae(1)) / (2 * (a + b)) * (x + a + b);
%!   group_idx_theta = find([[mesh.groups.tria6h].id] == 1);
%!   group_idx_master = find([[mesh.groups.tria6h].id] == 2);
%!   group_idx_slave = find([[mesh.groups.tria6h].id] == 3);
%!   nodes_constr = unique([[mesh.groups.tria6h(group_idx_theta)].nodes]);
%!   elem_constr.sfncon6h.master = mesh.elements.tria6h(mesh.groups.tria6h(group_idx_master).elements, :);
%!   elem_constr.sfncon6h.slave = mesh.groups.tria6h(group_idx_slave).nodes(:);
%!   elem_constr.sfncon6h.maxdist = 1e-10 * a;
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   thermal_constr_surf = fem_pre_mesh_constr_surf_to_node(mesh.nodes, elem_constr, load_case.domain);
%!   for i=1:numel(elem_constr.sfncon6h.slave)
%!     idx = find(nodes_constr == elem_constr.sfncon6h.slave(i));
%!     nodes_constr(idx) = -1;
%!   endfor
%!   nodes_constr = nodes_constr(nodes_constr > 0);
%!   mesh.elements.thermal_constr = struct("C", mat2cell(ones(1, numel(nodes_constr)), 1, ones(1, numel(nodes_constr))), ...
%!                                   "nodes", mat2cell(nodes_constr, 1, ones(1, numel(nodes_constr))));
%!   load_case.thermal_constr = struct("theta", mat2cell(theta_s(nodes_constr).', 1, ones(1, numel(nodes_constr))));
%!   mesh.elements.thermal_constr(end + (1:numel(thermal_constr_surf))) = thermal_constr_surf;
%!   load_case.thermal_constr(end + (1:numel(thermal_constr_surf))) = struct("theta", mat2cell(zeros(1, numel(thermal_constr_surf)), 1, ones(1, numel(thermal_constr_surf))));
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   e1 = [1; 0.6; -0.3];
%!   e2 = [-0.5; -0.3; 0.8];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   e1 /= norm(e1);
%!   e2 /= norm(e2);
%!   e3 /= norm(e3);
%!   R = [e1, e2, e3];
%!   mesh.nodes = [mesh.nodes(:, 1:3) * R.', mesh.nodes(:, 4:6) * R.'];
%!   [mat_ass.Kk, ...
%!    mat_ass.Qc, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_THERMAL_COND, ...
%!                                         FEM_VEC_LOAD_THERMAL], ...
%!                                        load_case);
%!   opt_sol.refine_max_iter = int32(20);
%!   U = fem_sol_factor(mat_ass.Kk, opt_sol) \ mat_ass.Qc;
%!   sol.theta = U(dof_map.ndof);
%!   assert_simple(sol.theta, theta_s, 1e-3 * max(abs(thetae)));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
