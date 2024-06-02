## fem_pre_mesh_import.m:108
%!test
%! try
%! ### TEST 108
%! ### Code_Aster TLV01 V4.25.001 01/02/2011
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
%!     dx = 5e-3;
%!     R = 0.1;
%!     r = 1e-6;
%!     Phi1 = -15 * pi / 180;
%!     Phi2 = 15 * pi / 180;
%!     Phi3 = 30 * pi / 180;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R=%g;\n", R);
%!     fprintf(fd, "r=%g;\n", r);
%!     fprintf(fd, "Phi1=%g;\n", Phi1);
%!     fprintf(fd, "Phi2=%g;\n", Phi2);
%!     fprintf(fd, "Phi3=%g;\n", Phi3);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Phi4 = 0.5*(Phi1+Phi2)+Pi/2;\n");
%!     fputs(fd, "Point(1) = {r * Cos(Phi1),r * Sin(Phi1),0,dx};\n");
%!     fputs(fd, "Point(2) = {R * Cos(Phi1),R * Sin(Phi1),0,dx};\n");
%!     fputs(fd, "Point(3) = {R * Cos(Phi2),R * Sin(Phi2),0,dx};\n");
%!     fputs(fd, "Point(4) = {r * Cos(Phi2),r * Sin(Phi2),0,dx};\n");
%!     fputs(fd, "Point(5) = {0, 0, 0};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Circle(2) = {2,5,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Circle(4) = {4,5,1};\n");
%!     fputs(fd, "Curve Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {{Cos(Phi4),Sin(Phi4),0},{0,0,0},Phi3} {\n");
%!     fputs(fd, "  Surface{6}; Layers{Ceil(R*Phi3/dx)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {1};\n");
%!     fputs(fd, "Physical Surface(\"convection\",2) = {8};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
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
%!   opt_mesh.elem_type = {"tria6h", "quad8", "iso20r", "penta15"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   lambda = 48.822;
%!   rho = 7200;
%!   cp = 669;
%!   mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh.material_data.cp = cp;
%!   thetae = 1000;
%!   theta0 = 20;
%!   h = 232.5;
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   grp_idx_conv = find([mesh.groups.quad8.id] == 2);
%!   mesh.elements.convection.quad8.nodes = mesh.elements.quad8([mesh.groups.quad8(grp_idx_conv).elements], :);
%!   mesh.elements.convection.quad8.h = repmat(h, size(mesh.elements.convection.quad8.nodes));
%!   load_case.convection.quad8.theta = repmat(thetae, size(mesh.elements.convection.quad8.nodes));
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   e1 = [1; 0.5; 0.3];
%!   e2 = [0.1; 1; 0.2];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   e1 /= norm(e1);
%!   e2 /= norm(e2);
%!   e3 /= norm(e3);
%!   R = [e1, e2, e3];
%!   [x, idx_x] = sort(mesh.nodes(:, 1));
%!   mesh.nodes = [mesh.nodes(:, 1:3) * R.', mesh.nodes(:, 4:6) * R.'];
%!   [mat_ass.C, ...
%!    mat_ass.Kk, ...
%!    mat_ass.Qc] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_HEAT_CAPACITY, ...
%!                                  FEM_MAT_THERMAL_COND, ...
%!                                  FEM_VEC_LOAD_THERMAL], ...
%!                                 load_case);
%!   dt = rho * cp * dx^2 / lambda;
%!   alpha = 0.5;
%!   sol.t = 0:dt:2400;
%!   sol.theta = zeros(dof_map.totdof, numel(sol.t));
%!   sol.theta(:, 1) = theta0;
%!   A = (1 / dt) * mat_ass.C + alpha * mat_ass.Kk;
%!   opts.number_of_threads = mbdyn_solver_num_threads_default();
%!   Afact = fem_sol_factor(A, opts);
%!   f = interp1([0, dt, sol.t(end)],[1, 1, 1], sol.t, "linear");
%!   for i=2:numel(sol.t)
%!     Qci = mat_ass.Qc * (f(i) * alpha + f(i - 1) * (1 - alpha));
%!     sol.theta(:, i) = Afact \ (mat_ass.C * (sol.theta(:, i - 1)) / dt - mat_ass.Kk * (sol.theta(:, i - 1) * (1 - alpha)) + Qci);
%!   endfor
%!   node_idx_center = find(sqrt(mesh.nodes(:,1).^2 + mesh.nodes(:,2).^2 + mesh.nodes(:,3).^2) == r);
%!   theta_center = mean(sol.theta(node_idx_center, :), 1);
%!   node_idx_surface = mesh.groups.quad8(grp_idx_conv).nodes;
%!   theta_surface = mean(sol.theta(node_idx_surface, :), 1);
%!   t_ref = (400:200:2400).';
%!   theta_ref_center = [334
%!                       500
%!                       618
%!                       706
%!                       774
%!                       828
%!                       872
%!                       902
%!                       923
%!                       942
%!                       956];
%!   theta_ref_surface = [461
%!                        608
%!                        696
%!                        774
%!                        828
%!                        868
%!                        902
%!                        923
%!                        942
%!                        956
%!                        962];
%!   if (do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(t_ref, theta_ref_center, "-;reference;k");
%!     plot(sol.t, theta_center, "-;solution;r");
%!     xlabel("t [s]");
%!     ylabel("theta [degC]");
%!     grid on;
%!     grid minor on;
%!     title("temperature in the center");
%!     figure("visible", "off");
%!     hold on;
%!     plot(t_ref, theta_ref_surface, "-;reference;k");
%!     plot(sol.t, theta_surface, "-;solution;r");
%!     xlabel("t [s]");
%!     ylabel("theta [degC]");
%!     grid on;
%!     grid minor on;
%!     title("temperature at the surface");
%!   endif
%!   tol = 1.5e-2;
%!   assert_simple(interp1(sol.t, theta_surface, t_ref, "linear", "extrap"), theta_ref_surface, tol * abs(thetae - theta0));
%!   assert_simple(interp1(sol.t, theta_center, t_ref, "linear", "extrap"), theta_ref_center, tol * abs(thetae - theta0));
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
