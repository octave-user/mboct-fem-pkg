## fem_pre_mesh_import.m:50
%!test
%! try
%! ### TEST 50
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
%!     Ec = 210000e6;
%!     Ei = 125000e6;
%!     CTEc = 12.5e-6;
%!     CTEi = 16.7e-6;
%!     dT = [100, 200];
%!     H = 1e-3;
%!     r1 = 0;
%!     h1 = H;
%!     h2 = H;
%!     w = 1.25e-3;
%!     h = 1.25e-3;
%!     L = 100e-3;
%!     r2 = r1 + h1;
%!     r3 = r2 + h2;
%!     K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!     Uy_ref = -3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w = %g;\n", w);
%!     fprintf(fd, "r1 = %g;\n", r1);
%!     fprintf(fd, "r2 = %g;\n", r2);
%!     fprintf(fd, "r3 = %g;\n", r3);
%!     fprintf(fd, "h = %g;\n", h);
%!     fprintf(fd, "L = %g;\n", L);
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
%!     fputs(fd, "tmp1[] = Extrude {L,0,0}{ Surface{11}; };\n");
%!     fputs(fd, "tmp2[] = Extrude {L,0,0}{ Surface{12}; };\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {11,12};\n");
%!     fputs(fd, "Physical Surface(\"master\",6) = {tmp1[4]};\n");
%!     fputs(fd, "Physical Surface(\"slave\",7) = {tmp2[2]};\n");
%!     fputs(fd, "Physical Surface(\"end-section\",8) = {tmp1[0],tmp2[0]};\n");
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
%!   opts.elem_type = {"tria6", "tet10"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opts));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.tet10.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.tet10.id] == 2);
%!   grp_id_master = find([mesh.groups.tria6.id] == 6);
%!   grp_id_slave = find([mesh.groups.tria6.id] == 7);
%!   grp_id_end = find([mesh.groups.tria6.id] == 8);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.tet10(mesh.groups.tet10(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.tet10(mesh.groups.tet10(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).gamma = CTEc;
%!   mesh.material_data(1).rho = 7850;
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900;
%!   mesh.material_data(2).gamma = CTEi;
%!   mesh.elements.sfncon6.slave = mesh.groups.tria6(grp_id_slave).nodes(:);
%!   mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(grp_id_master).elements, :);
%!   mesh.elements.sfncon6.maxdist = 1e-6;
%!   constr.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 3);
%!   node_id_clamp = mesh.groups.tria6(grp_id_clamp).nodes;
%!   idx_clamp = true(1, numel(node_id_clamp));
%!   for i=1:numel(mesh.elements.sfncon6.slave)
%!     idx_clamp(find(node_id_clamp == mesh.elements.sfncon6.slave(i))) = false;
%!   endfor
%!   node_id_clamp = node_id_clamp(idx_clamp);
%!   mesh.elements.joints = struct("nodes", mat2cell(node_id_clamp,1,ones(1,numel(node_id_clamp))), "C", repmat({[eye(3),zeros(3,3)]}, 1, numel(node_id_clamp)));
%!   load_case = struct("dTheta", cell(1, numel(dT)));
%!   for i=1:numel(dT)
%!     load_case(i).dTheta = repmat(dT(i), rows(mesh.nodes), 1);
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, constr);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress, ...
%!    sol_stat.strain] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH, ...
%!                                       FEM_VEC_STRAIN_TOTAL], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   load_case2 = struct("epsilon0", cell(1, numel(load_case)));
%!   for k=1:numel(load_case)
%!     load_case2(k).epsilon0.tet10 = sol_stat.strain.epsilon.tet10(:, :, :, k);
%!   endfor
%!   [mat_ass2.K, ...
%!    mat_ass2.R] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS, ...
%!                                  FEM_VEC_LOAD_CONSISTENT], ...
%!                                 load_case2);
%!   [sol_stat2] = fem_sol_static(mesh, dof_map, mat_ass2);
%!   [sol_stat2.stress, ...
%!    sol_stat2.strain] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_VEC_STRESS_CAUCH, ...
%!                                        FEM_VEC_STRAIN_TOTAL], ...
%!                                       load_case2, ...
%!                                       sol_stat2);
%!   Uy = reshape(mean(sol_stat.def(mesh.groups.tria6(grp_id_end).nodes, 2, :)), 1, numel(dT));
%!   Uy2 = reshape(mean(sol_stat2.def(mesh.groups.tria6(grp_id_end).nodes, 2, :)), 1, numel(dT));
%!   tol = 2e-2;
%!   assert_simple(Uy, Uy_ref, tol * abs(Uy_ref));
%!   assert_simple(Uy2, Uy_ref, tol * abs(Uy_ref));
%!   assert_simple(sol_stat2.def, sol_stat.def, 1e-7 * max(max(max(abs(sol_stat.def)))));
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
