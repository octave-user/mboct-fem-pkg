## fem_pre_mesh_import.m:374
%!test
%! try
%! ### TEST 374
%! ## bi-metal bending iso20r
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
%!     N = 3;
%!     L = 100e-3;
%!     b = 0.5e-3;
%!     H = 1e-3;
%!     h = 0.5e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Round(L / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round(b / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(L / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round(b / h) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Round(H/h)}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Round(H/h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso20r", "penta15", "quad8r", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.iso20r.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.iso20r.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.iso20r(mesh.groups.iso20r(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.iso20r(mesh.groups.iso20r(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6;
%!   Ei = 125000e6;
%!   CTEc = 12.5e-6;
%!   CTEi = 16.7e-6;
%!   dT = 100;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).gamma = CTEc;
%!   mesh.material_data(1).rho = 7850;
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900;
%!   mesh.material_data(2).gamma = CTEi;
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.quad8r(find([[mesh.groups.quad8r].id] == 3)).nodes, :) = true;
%!   load_case.dTheta = repmat(dT, rows(mesh.nodes), 1);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_SCA_TOT_MASS], ...
%!                                load_case);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   grp_id_load = find([mesh.groups.quad8r.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.quad8r(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
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
