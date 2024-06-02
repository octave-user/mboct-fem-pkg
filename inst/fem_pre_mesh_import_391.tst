## fem_pre_mesh_import.m:39
%!test
%! try
%! ## TEST 39: thick walled cylinder with internal pressure iso20r
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
%!     ro = 50e-3;
%!     ri = 25e-3;
%!     L = 20e-3;
%!     dT = 100;
%!     Pi = 2e8;
%!     sigmaz = Pi * ri^2 / (ro^2 - ri^2);
%!     E = 2.1e11;
%!     nu = 0.3;
%!     gamma = 0.12e-4;
%!     rho = 7850;
%!     mesh_size = 2.5e-3;
%!     scale_def = 10e-3;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "ro = %g;\n", ro);
%!     fprintf(fd, "ri = %g;\n", ri);
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "H = %g;\n", mesh_size);
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 1;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold = 0.9;\n");
%!     fputs(fd, "Point(1) = {0,ro,0,H};\n");
%!     fputs(fd, "Point(2) = {0,ro,L,H};\n");
%!     fputs(fd, "Point(3) = {0,ri,L,H};\n");
%!     fputs(fd, "Point(4) = {0,ri,0,H};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude{{0,0,1},{0,0,0}, -Pi/2}{ Surface{6}; Layers{Ceil(ro*Pi/2/H)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-y\", 2) = {tmp[0]};\n");
%!     fputs(fd, "Physical Surface(\"clamp-x\", 3) = {6};\n");
%!     fputs(fd, "Physical Surface(\"clamp-z\", 4) = {tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"pressure-radial\", 5) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"pressure-axial\", 6) = {tmp[3]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso20r", "quad8"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_y = find([[mesh.groups.quad8].id] == 2);
%!   grp_id_clamp_x = find([[mesh.groups.quad8].id] == 3);
%!   grp_id_clamp_z = find([[mesh.groups.quad8].id] == 4);
%!   grp_id_Pi = find([mesh.groups.quad8.id] == 5);
%!   grp_id_sigmaz = find([mesh.groups.quad8.id] == 6);
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp_x).nodes, 1) = true;
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp_y).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp_z).nodes, 3) = true;
%!   elem_Pi = mesh.elements.quad8(mesh.groups.quad8(grp_id_Pi).elements, :);
%!   elem_sigmaz = mesh.elements.quad8(mesh.groups.quad8(grp_id_sigmaz).elements, :);
%!   load_case.pressure.quad8.elements = [elem_Pi; elem_sigmaz];
%!   load_case.pressure.quad8.p = [repmat(Pi, size(elem_Pi)); repmat(-sigmaz, size(elem_sigmaz))];
%!   mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.material_data.gamma = gamma;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   load_case.dTheta = repmat(dT, rows(mesh.nodes), 1);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   r = sqrt(mesh.nodes(:, 1).^2 + mesh.nodes(:, 2).^2);
%!   z = mesh.nodes(:, 3);
%!   Theta = atan2(mesh.nodes(:, 2), mesh.nodes(:, 1));
%!   sigmaR = Pi * ri^2 / (ro^2 - ri^2) * (1 - ro^2 ./ r.^2);
%!   sigmaTheta = Pi * ri^2 / (ro^2 - ri^2) * (1 + ro^2 ./ r.^2);
%!   sigmaZ = Pi * ri^2 / (ro^2 - ri^2);
%!   epsilonR = 1 / E * (sigmaR - nu * (sigmaTheta + sigmaZ));
%!   epsilonTheta = 1 / E * (sigmaTheta - nu * (sigmaR + sigmaZ));
%!   epsilonZ = 1 / E * (sigmaZ - nu * (sigmaR + sigmaTheta)) + gamma * dT;
%!   Uz = z .* epsilonZ;
%!   sigmaX = sigmaR .* cos(Theta) - sigmaTheta .* sin(Theta);
%!   sigmaY = sigmaR .* sin(Theta) + sigmaTheta .* cos(Theta);
%!   theta = Theta(mesh.elements.iso20r);
%!   sigmar = sigmatheta = sigmaz = zeros(size(mesh.elements.iso20r));
%!   idxtens = int32([1, 4, 6;
%!                    4, 2, 5;
%!                    6, 5, 3]);
%!   for i=1:rows(mesh.elements.iso20r)
%!     for j=1:columns(mesh.elements.iso20r)
%!       e1 = mesh.nodes(mesh.elements.iso20r(i, j), 1:3).';
%!       e1(3) = 0;
%!       e3 = [0; 0; 1];
%!       e2 = cross(e3, e1);
%!       e1 /= norm(e1);
%!       e2 /= norm(e2);
%!       e3 /= norm(e3);
%!       R = [e1, e2, e3];
%!       SIGMAIJ = R.' * sol_stat.stress.taum.iso20r(i, j, :)(idxtens) * R;
%!       sigmar(i, j) = SIGMAIJ(1, 1);
%!       sigmatheta(i, j) = SIGMAIJ(2, 2);
%!       sigmaz(i, j) = SIGMAIJ(3, 3);
%!     endfor
%!   endfor
%!   tol = 3e-2;
%!   assert_simple(sigmar, sigmaR(mesh.elements.iso20r), tol * max(abs(sigmaR)));
%!   assert_simple(sigmatheta, sigmaTheta(mesh.elements.iso20r), tol * max(abs(sigmaTheta)));
%!   assert_simple(sigmaz, repmat(sigmaZ, size(sigmaz)), tol * abs(sigmaZ));
%!   assert_simple(sol_stat.def(:, 3), Uz, tol * max(abs(Uz)));
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
