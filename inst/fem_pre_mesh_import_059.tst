## fem_pre_mesh_import.m:59
%!test
%! ### TEST 59
%! ## The 1-D Heat Equation
%! ## 18.303 Linear Partial Differential Equations
%! ## Matthew J. Hancock
%! ## Fall 2006
%! ## https://ocw.mit.edu/courses/mathematics/18-303-linear-partial-differential-equations-fall-2006/lecture-notes/heateqni.pdf
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
%!     l = 10e-3;
%!     dx = 0.5e-3;
%!     w = dx;
%!     h = dx;
%!     lambda = 50;
%!     rho = 7850;
%!     cp = 465;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "l=%g;\n", l);
%!     fprintf(fd, "w=%g;\n", w);
%!     fprintf(fd, "h=%g;\n", h);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {0,w,0,dx};\n");
%!     fputs(fd, "Point(3) = {0,w,h,dx};\n");
%!     fputs(fd, "Point(4) = {0,0,h,dx};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {l,0,0} {\n");
%!     fputs(fd, "  Surface{6};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"bnd1\",1) = {tmp[0]};\n");
%!     fputs(fd, "Physical Surface(\"bnd2\",2) = {6};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_THERMAL;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.k = diag([lambda, lambda, lambda]);
%!   mesh.material_data.cp = cp;
%!   u0 = 100;
%!   ub = 50;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.Kk, ...
%!    mat_ass.C, ...
%!    mat_ass.coll_Kk, ...
%!    mat_ass.coll_C] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_THERMAL_COND, ...
%!                                      FEM_MAT_HEAT_CAPACITY, ...
%!                                      FEM_VEC_COLL_THERMAL_COND, ...
%!                                      FEM_VEC_COLL_HEAT_CAPACITY], ...
%!                                     load_case);
%!   idx_b = [mesh.groups.tria6h.nodes];
%!   idx_i = true(rows(mesh.nodes), 1);
%!   idx_i(idx_b) = false;
%!   idx_i = find(idx_i);
%!   x = mesh.nodes(:, 1);
%!   Kk11 = mat_ass.Kk(idx_i, idx_i);
%!   Kk12 = mat_ass.Kk(idx_i, idx_b);
%!   C11 = mat_ass.C(idx_i, idx_i);
%!   theta_b = repmat(ub, numel(idx_b), 1);
%!   theta0 = repmat(u0, dof_map.totdof, 1);
%!   qref = mesh.material_data.rho * mesh.material_data.cp * l * w * h;
%!   assert_simple(sum(sum(mat_ass.C)), qref, eps^0.5 * abs(qref));
%!   dt = rho * cp * dx^2 / lambda;
%!   alpha = 0.6;
%!   T_ = l^2 / lambda * rho * cp;
%!   sol.t = 0:dt:T_;
%!   sol.theta = zeros(dof_map.totdof, numel(sol.t));
%!   sol.theta(:, 1) = theta0;
%!   A = (1 / dt) * C11 + alpha * Kk11;
%!   opts.number_of_threads = mbdyn_solver_num_threads_default();
%!   opts.solver = "pastix";
%!   Afact = fem_sol_factor(A, opts);
%!   for i=2:numel(sol.t)
%!     sol.theta(idx_i, i) = Afact \ (C11 * (sol.theta(idx_i, i - 1) / dt) - Kk11 * (sol.theta(idx_i, i - 1) * (1 - alpha)) - Kk12 * theta_b);
%!     sol.theta(idx_b, i) = theta_b;
%!   endfor
%!   u0_ = 1;
%!   n = 1:10000;
%!   Bn_ = -2 * u0_ ./ (n * pi) .* ((-1).^n - 1);
%!   x_ = x / l;
%!   u_ = zeros(numel(x), numel(sol.t));
%!   for n=1:numel(Bn_)
%!     t_ = sol.t / T_;
%!     u_ += Bn_(n) * sin(n * pi * x_) .* exp(-n^2 * pi^2 * t_);
%!   endfor
%!   sol.theta_ref = u_ * (u0 - ub) + ub;
%!   [x, idx_theta] = sort(x);
%!   if (do_plot)
%!     figure("visible", "off");
%!     hold("on");
%!     plot(sol.t, max(sol.theta_ref, [], 1), "-;max(ref);k");
%!     plot(sol.t, max(sol.theta, [], 1), "-;max(sol);r");
%!     plot(sol.t, min(sol.theta_ref, [], 1), "--;min(ref);k");
%!     plot(sol.t, min(sol.theta, [], 1), "--;min(sol);r");
%!     xlabel("t [s]");
%!     ylabel("theta [degC]");
%!     grid on;
%!     grid minor on;
%!     title("temperature versus time");
%!     for i=[1,2:5,6:100:numel(sol.t)]
%!       figure("visible", "off");
%!       hold("on");
%!       plot(x, sol.theta_ref(idx_theta, i), "-;ref;k");
%!       plot(x, sol.theta(idx_theta, i), "-;sol;r");
%!       xlabel("t [s]");
%!       ylabel("theta [degC]");
%!       title(sprintf("temperature versus x at time = %.2fs", sol.t(i)));
%!       grid on;
%!       grid minor on;
%!     endfor
%!   endif
%!   tol = 1e-2;
%!   assert_simple(sol.theta(:, 10:end), sol.theta_ref(:, 10:end), tol * abs(u0 - ub));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
