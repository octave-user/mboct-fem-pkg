## fem_pre_mesh_import.m:48
%!test
%! ### TEST 48
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
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     r1 = 10e-3;
%!     h1 = 5e-3;
%!     h2 = 5e-3;
%!     w = 10e-3;
%!     h = 1.25e-3;
%!     r2 = r1 + h1;
%!     r3 = r2 + h2;
%!     fprintf(fd, "w = %g;\n", w);
%!     fprintf(fd, "r1 = %g;\n", r1);
%!     fprintf(fd, "r2 = %g;\n", r2);
%!     fprintf(fd, "r3 = %g;\n", r3);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0,r3,0,h};\n");
%!     fputs(fd, "Point(2) = {0,r3,w,h};\n");
%!     fputs(fd, "Point(3) = {0,r2,w,h};\n");
%!     fputs(fd, "Point(4) = {0,r2,0,h};\n");
%!     fputs(fd, "Point(5) = {0,r2,0,h};\n");
%!     fputs(fd, "Point(6) = {0,r2,w,h};\n");
%!     fputs(fd, "Point(7) = {0,r1,w,h};\n");
%!     fputs(fd, "Point(8) = {0,r1,0,h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line(5) = {5,6};\n");
%!     fputs(fd, "Line(6) = {6,7};\n");
%!     fputs(fd, "Line(7) = {7,8};\n");
%!     fputs(fd, "Line(8) = {8,5};\n");
%!     fputs(fd, "Line Loop(9) = {1,2,3,4};\n");
%!     fputs(fd, "Line Loop(10) = {5,6,7,8};\n");
%!     fputs(fd, "Plane Surface(11) = {9};\n");
%!     fputs(fd, "Plane Surface(12) = {10};\n");
%!     fputs(fd, "tmp1[] = Extrude {{0,0,1},{0,0,0},-Pi/2} { Surface{11}; };\n");
%!     fputs(fd, "tmp2[] = Extrude {{0,0,1},{0,0,0},-Pi/2} { Surface{12}; };\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-x\",3) = {11,12};\n");
%!     fputs(fd, "Physical Surface(\"clamp-y\",4) = {tmp1[0],tmp2[0]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-z\",5) = {tmp1[5],tmp2[5]};\n");
%!     fputs(fd, "Physical Surface(\"master\",6) = {tmp1[4]};\n");
%!     fputs(fd, "Physical Surface(\"slave\",7) = {tmp2[2]};\n");
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
%!   opts.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opts));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.tet10h.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.tet10h.id] == 2);
%!   grp_id_master = find([mesh.groups.tria6h.id] == 6);
%!   grp_id_slave = find([mesh.groups.tria6h.id] == 7);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.tet10h(mesh.groups.tet10h(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.tet10h(mesh.groups.tet10h(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   dT = 100;
%!   mesh.material_data(1).E = 210000e6;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).gamma = 12.5e-6;
%!   mesh.material_data(1).rho = 7850;
%!   mesh.material_data(2).E = 125000e6;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900;
%!   mesh.material_data(2).gamma = 16.7e-6;
%!   mesh.elements.sfncon6h.slave = mesh.groups.tria6h(grp_id_slave).nodes(:);
%!   mesh.elements.sfncon6h.master = mesh.elements.tria6h(mesh.groups.tria6h(grp_id_master).elements, :);
%!   mesh.elements.sfncon6h.maxdist = 1e-6;
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.tria6h].id] == 3);
%!   grp_id_clamp_y = find([[mesh.groups.tria6h].id] == 4);
%!   grp_id_clamp_z = find([[mesh.groups.tria6h].id] == 5);
%!   node_id_clamp_x = mesh.groups.tria6h(grp_id_clamp_x).nodes;
%!   node_id_clamp_y = mesh.groups.tria6h(grp_id_clamp_y).nodes;
%!   node_id_clamp_z = mesh.groups.tria6h(grp_id_clamp_z).nodes;
%!   idx_x = true(1, numel(node_id_clamp_x));
%!   idx_y = true(1, numel(node_id_clamp_y));
%!   idx_z = true(1, numel(node_id_clamp_z));
%!   for i=1:numel(mesh.elements.sfncon6h.slave)
%!     idx_x(find(node_id_clamp_x == mesh.elements.sfncon6h.slave(i))) = false;
%!     idx_y(find(node_id_clamp_y == mesh.elements.sfncon6h.slave(i))) = false;
%!     idx_z(find(node_id_clamp_z == mesh.elements.sfncon6h.slave(i))) = false;
%!   endfor
%!   node_id_clamp_x = node_id_clamp_x(idx_x);
%!   node_id_clamp_y = node_id_clamp_y(idx_y);
%!   node_id_clamp_z = node_id_clamp_z(idx_z);
%!   mesh.elements.joints = [struct("nodes", mat2cell(node_id_clamp_x,1,ones(1,numel(node_id_clamp_x))), "C", repmat({[1,0,0,0,0,0]}, 1, numel(node_id_clamp_x))), ...
%!                           struct("nodes", mat2cell(node_id_clamp_y,1,ones(1,numel(node_id_clamp_y))), "C", repmat({[0,1,0,0,0,0]}, 1, numel(node_id_clamp_y))), ...
%!                           struct("nodes", mat2cell(node_id_clamp_z,1,ones(1,numel(node_id_clamp_z))), "C", repmat({[0,0,1,0,0,0]}, 1, numel(node_id_clamp_z)))];
%!   load_case.dTheta = repmat(dT, rows(mesh.nodes), 1);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   opt_sol.refine_max_iter = 250;
%!   opt_sol.verbose = false;
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   [sol_stat.stress, ...
%!    sol_stat.strain] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH, ...
%!                                       FEM_VEC_STRAIN_TOTAL], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   load_case2.epsilon0 = sol_stat.strain.epsilon;
%!   [mat_ass2.K, ...
%!    mat_ass2.R] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS, ...
%!                                  FEM_VEC_LOAD_CONSISTENT], ...
%!                                 load_case2);
%!   [sol_stat2] = fem_sol_static(mesh, dof_map, mat_ass2, opt_sol);
%!   [sol_stat2.stress, ...
%!    sol_stat2.strain] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_VEC_STRESS_CAUCH, ...
%!                                        FEM_VEC_STRAIN_TOTAL], ...
%!                                       load_case2, ...
%!                                       sol_stat2);
%!   tol_def = 5e-5;
%!   tol_strain = 1e-3;
%!   tol_stress = 5e-3;
%!   assert_simple(sol_stat2.def, sol_stat.def, tol_def * max(max(abs(sol_stat.def))));
%!   assert_simple(sol_stat2.stress.tau.tet10h, zeros(size(sol_stat.stress.tau.tet10h)), tol_stress * max(max(max(abs(sol_stat.stress.tau.tet10h)))));
%!   assert_simple(sol_stat2.strain.epsilonm.tet10h, sol_stat.strain.epsilonm.tet10h, tol_strain * max(max(max(abs(sol_stat.strain.epsilon.tet10h)))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
