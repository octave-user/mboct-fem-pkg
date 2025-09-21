## fem_pre_mesh_import.m:449
%!test
%! try
%!   ## MBDyn rigid body dynamics
%!   n = [1
%!        5
%!        10
%!        20];
%!   fref = [0.34597008661255	0.345970086618602
%!           1.27458549458623	4.73428636065456
%!           2.95187224998486	9.87127398130082
%!           6.0973109171929	19.9361143740589];
%!   param.h = 100e-3;
%!   param.d = 10e-3;
%!   param.D = 280e-3;
%!   param.t = 30e-3;
%!   param.dx = 10e-3;
%!   param.g = 9.81;
%!   param.num_modes = 10;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     fd = -1;
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "h=%g;\n", param.h);
%!       fprintf(fd, "d=%g;\n", param.d);
%!       fprintf(fd, "D=%g;\n", param.D);
%!       fprintf(fd, "t=%g;\n", param.t);
%!       fprintf(fd, "dx=%g;\n", param.dx);
%!       fputs(fd, "Point(1) = {0, 0, 0};\n");
%!       fputs(fd, "Point(2) = {0.5 * d, 0, 0};\n");
%!       fputs(fd, "Point(3) = {0.5 * d, 0, 0.5 * (h - t)};\n");
%!       fputs(fd, "Point(4) = {0.5 * D, 0, 0.5 * (h - t)};\n");
%!       fputs(fd, "Point(5) = {0.5 * D, 0, 0.5 * (h + t)};\n");
%!       fputs(fd, "Point(6) = {0.5 * d, 0, 0.5 * (h + t)};\n");
%!       fputs(fd, "Point(7) = {0.5 * d, 0, h};\n");
%!       fputs(fd, "Point(8) = {0, 0, h};\n");
%!       fputs(fd, "Point(9) = {0, 0, (h+t)/2};\n");
%!       fputs(fd, "Point(10) = {0, 0, (h-t)/2};\n");
%!       fputs(fd, "Line(1) = {1, 2};\n");
%!       fputs(fd, "Line(2) = {2, 3};\n");
%!       fputs(fd, "Line(3) = {3, 10};\n");
%!       fputs(fd, "Line(4) = {10, 1};\n");
%!       fputs(fd, "Line(5) = {3, 6};\n");
%!       fputs(fd, "Line(6) = {6, 9};\n");
%!       fputs(fd, "Line(7) = {9, 10};\n");
%!       fputs(fd, "Line(8) = {3, 4};\n");
%!       fputs(fd, "Line(9) = {4,5};\n");
%!       fputs(fd, "Line(10) = {5,6};\n");
%!       fputs(fd, "Line(11) = {6,7};\n");
%!       fputs(fd, "Line(12) = {7,8};\n");
%!       fputs(fd, "Line(13) = {8,9};\n");
%!       fputs(fd, "Transfinite Curve(1) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(2) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(3) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(4) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(5) = Round(t/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(6) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(7) = Round(t/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(8) = Round((D-d)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(9) = Round(t/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(10) = Round((D-d)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(11) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(12) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(13) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Line Loop(14) = {1,2,3,4};\n");
%!       fputs(fd, "Line Loop(15) = {-3,5,6,7};\n");
%!       fputs(fd, "Line Loop(16) = {8,9,10,-5};\n");
%!       fputs(fd, "Line Loop(17) = {11,12,13,6};\n");
%!       fputs(fd, "Plane Surface(18) = {14};\n");
%!       fputs(fd, "Plane Surface(19) = {15};\n");
%!       fputs(fd, "Plane Surface(20) = {16};\n");
%!       fputs(fd, "Plane Surface(21) = {17};\n");
%!       fputs(fd, "Transfinite Surface(18)={1,2,3,10};\n");
%!       fputs(fd, "Transfinite Surface(19)={10,3,6,9};\n");
%!       fputs(fd, "Transfinite Surface(20)={3,4,5,6};\n");
%!       fputs(fd, "Transfinite Surface(21)={6,7,8,9};\n");
%!       fputs(fd, "Recombine Surface{18};\n");
%!       fputs(fd, "Recombine Surface{19};\n");
%!       fputs(fd, "Recombine Surface{20};\n");
%!       fputs(fd, "Recombine Surface{21};\n");
%!       fputs(fd, "v1[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{18};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v2[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{19};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v3[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{20};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v4[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{21};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{18,v1[0]};\n");
%!       fputs(fd, "Recombine Surface{19,v2[0]};\n");
%!       fputs(fd, "Recombine Surface{20,v3[0]};\n");
%!       fputs(fd, "Recombine Surface{21,v4[0]};\n");
%!       fputs(fd, "v5[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v1[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v6[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v2[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v7[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v3[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v8[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v4[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{v1[0],v5[0]};\n");
%!       fputs(fd, "Recombine Surface{v2[0],v6[0]};\n");
%!       fputs(fd, "Recombine Surface{v3[0],v7[0]};\n");
%!       fputs(fd, "Recombine Surface{v4[0],v8[0]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {v1[1],v2[1],v3[1],v4[1],v5[1],v6[1],v7[1],v8[1]};\n");
%!       fputs(fd, "Physical Surface(\"support\", 2) = {22, 39};\n");
%!       fputs(fd, "Mesh.ElementOrder = 2;\n");
%!       fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!       fputs(fd, "Mesh.Algorithm = 3;\n");
%!       fputs(fd, "Mesh 3;\n");
%!       fputs(fd, "Coherence Mesh;\n");
%!       fputs(fd, "Mesh.Format = 1;\n");
%!       fprintf(fd, "Save \"%s.msh\";\n", filename);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!                      #spawn_wait(spawn("gmsh", {[filename, ".geo"]}));
%!                      #return
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", "-order", "2", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     [~] = unlink([filename, ".geo"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     mesh.material_data.E = 210000e6*10;
%!     mesh.material_data.nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     mesh.materials.penta15 = ones(rows(mesh.elements.penta15),1, "int32");
%!     mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!     node_idx_tip = rows(mesh.nodes) + 1;
%!     mesh.nodes(node_idx_tip, :) = zeros(1, 6);
%!     mesh.elements.rbe2 = fem_pre_mesh_rbe2_from_surf(mesh, 2, node_idx_tip);
%!     mesh.elements.joints(1).nodes = node_idx_tip;
%!     mesh.elements.joints(1).C = eye(6)([1:3,6],:);
%!     mesh.elements.springs(1).nodes = node_idx_tip;
%!     mesh.elements.springs(1).K = eye(6);
%!     mesh.elements.springs(1).F = zeros(6,1);
%!     load_case_dof.locked_dof = false(size(mesh.nodes));
%!     load_case_stat.g = [0; 0; -param.g];
%!     opt_sol.symmetric = false;
%!     opt_sol.solver = "pardiso";
%!     opt_sol.refine_max_iter = 250;
%!     dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!     for i=1:numel(n)
%!       load_case_stat.omega = [0; 0; 2 * pi * n(i)];
%!       [mat_ass_stat.K, ...
%!        mat_ass_stat.KOMEGA, ...
%!        mat_ass_stat.R, ...
%!        mat_ass_stat.dm, ...
%!        mat_ass_stat.S, ...
%!        mat_ass_stat.J, ...
%!        mat_ass_stat.mat_info, ...
%!        mat_ass_stat.mesh_info] = fem_ass_matrix(mesh, ...
%!                                                 dof_map, ...
%!                                                 [FEM_MAT_STIFFNESS, ... ...
%!                                                                     FEM_MAT_STIFFNESS_OMEGA, ...
%!                                                  FEM_VEC_LOAD_CONSISTENT, ...
%!                                                  FEM_SCA_TOT_MASS, ...
%!                                                  FEM_VEC_INERTIA_M1, ...
%!                                                  FEM_MAT_INERTIA_J], ...
%!                                                 load_case_stat);
%!       mat_ass_stat.K += mat_ass_stat.KOMEGA;
%!       sol_stat = fem_sol_static(mesh, dof_map, mat_ass_stat, opt_sol);
%!       sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_VEC_STRESS_CAUCH], ...
%!                                        load_case_stat, ...
%!                                        sol_stat);
%!       mesh.elements.springs(1).K(:, :) = 0;
%!       load_case_dyn.g = load_case_stat.g;
%!       load_case_dyn.omega = [0; 0; 2 * pi * n(i)];
%!       load_case_dyn.tau0 = sol_stat.stress.tau;
%!       load_case_dyn.lambda = sol_stat.lambda;
%!       [mat_ass_dyn.M, ...
%!        mat_ass_dyn.K, ...
%!        mat_ass_dyn.D, ...
%!        mat_ass_dyn.KTAU0, ...
%!        mat_ass_dyn.KOMEGA, ...
%!        mat_ass_dyn.KOMEGADOT, ...
%!        mat_ass_dyn.DOMEGA, ...
%!        mat_ass_dyn.R, ...
%!        mat_ass_dyn.mat_info, ...
%!        mat_ass_dyn.mesh_info] = fem_ass_matrix(mesh, ...
%!                                                dof_map, ...
%!                                                [FEM_MAT_MASS, ...
%!                                                 FEM_MAT_STIFFNESS, ...
%!                                                 FEM_MAT_DAMPING, ...
%!                                                 FEM_MAT_STIFFNESS_TAU0, ...
%!                                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                                 FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                                load_case_dyn);
%!       mat_ass_dyn.K += mat_ass_dyn.KTAU0;
%!       mat_ass_dyn.K += mat_ass_dyn.KOMEGA;
%!       mat_ass_dyn.K += mat_ass_dyn.KOMEGADOT;
%!       mat_ass_dyn.D += mat_ass_dyn.DOMEGA;
%!       sol_eig = fem_sol_modal_damped(mesh, dof_map, mat_ass_dyn, param.num_modes, opt_sol);
%!       tol = 2e-2;
%!       printf("%5.2f\t", fref(i, :));
%!       printf("%5.2f\t", sol_eig.f((sol_eig.f > 0))(1:2));
%!       printf("\n");
%!       assert_simple(sol_eig.f(sol_eig.f > 0)(1:2), fref(i, :), tol * max(fref(i, :)));
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
