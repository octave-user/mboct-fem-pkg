## fem_pre_mesh_import.m:274
%!test
%! ### TEST 274
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
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opts.elem_type = {"tria10", "tet20"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opts));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.tet20 = zeros(rows(mesh.elements.tet20), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.tet20.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.tet20.id] == 2);
%!   grp_id_master = find([mesh.groups.tria10.id] == 6);
%!   grp_id_slave = find([mesh.groups.tria10.id] == 7);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.tet20(mesh.groups.tet20(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.tet20(mesh.groups.tet20(grp_id_mat2(i)).elements(:)) = 2;
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
%!   mesh.elements.sfncon10.slave = mesh.groups.tria10(grp_id_slave).nodes(:);
%!   mesh.elements.sfncon10.master = mesh.elements.tria10(mesh.groups.tria10(grp_id_master).elements, :);
%!   mesh.elements.sfncon10.maxdist = 1e-6;
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.tria10].id] == 3);
%!   grp_id_clamp_y = find([[mesh.groups.tria10].id] == 4);
%!   grp_id_clamp_z = find([[mesh.groups.tria10].id] == 5);
%!   node_id_clamp_x = mesh.groups.tria10(grp_id_clamp_x).nodes;
%!   node_id_clamp_y = mesh.groups.tria10(grp_id_clamp_y).nodes;
%!   node_id_clamp_z = mesh.groups.tria10(grp_id_clamp_z).nodes;
%!   idx_x = true(1, numel(node_id_clamp_x));
%!   idx_y = true(1, numel(node_id_clamp_y));
%!   idx_z = true(1, numel(node_id_clamp_z));
%!   for i=1:numel(mesh.elements.sfncon10.slave)
%!     idx_x(find(node_id_clamp_x == mesh.elements.sfncon10.slave(i))) = false;
%!     idx_y(find(node_id_clamp_y == mesh.elements.sfncon10.slave(i))) = false;
%!     idx_z(find(node_id_clamp_z == mesh.elements.sfncon10.slave(i))) = false;
%!   endfor
%!   node_id_clamp_x = node_id_clamp_x(idx_x);
%!   node_id_clamp_y = node_id_clamp_y(idx_y);
%!   node_id_clamp_z = node_id_clamp_z(idx_z);
%!   mesh.elements.joints = [struct("nodes", mat2cell(node_id_clamp_x,1,ones(1,numel(node_id_clamp_x))), "C", repmat({[1,0,0,0,0,0]}, 1, numel(node_id_clamp_x))), ...
%!                           struct("nodes", mat2cell(node_id_clamp_y,1,ones(1,numel(node_id_clamp_y))), "C", repmat({[0,1,0,0,0,0]}, 1, numel(node_id_clamp_y))), ...
%!                           struct("nodes", mat2cell(node_id_clamp_z,1,ones(1,numel(node_id_clamp_z))), "C", repmat({[0,0,1,0,0,0]}, 1, numel(node_id_clamp_z)))];
%!   load_case.dTheta = repmat(dT, rows(mesh.nodes), 1);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(100);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.verbose = int32(0);
%!   opt_sol.pre_scaling = true;
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
%!   tol_def = 1e-5;
%!   tol_strain = 1e-3;
%!   tol_stress = 5e-3;
%!   assert_simple(sol_stat2.def, sol_stat.def, tol_def * max(max(abs(sol_stat.def))));
%!   assert_simple(sol_stat2.stress.tau.tet20, zeros(size(sol_stat.stress.tau.tet20)), tol_stress * max(max(max(abs(sol_stat.stress.tau.tet20)))));
%!   assert_simple(sol_stat2.strain.epsilonm.tet20, sol_stat.strain.epsilonm.tet20, tol_strain * max(max(max(abs(sol_stat.strain.epsilon.tet20)))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
