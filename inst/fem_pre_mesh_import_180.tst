## fem_pre_mesh_import.m:180
%!test
%! ### TEST 180 - 2D wave equation
%! ##
%! ## 1.138J/2.062J/18.376J, WAVE PROPAGATION
%! ## Fall, 2004 MIT
%! ## Notes by C. C. Mei
%! ## CHAPTER THREE
%! ## TWO DIMENSIONAL WAVES
%! ## Reflection and tranmission of sound at an interface
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
%!   [fd, msg] = fopen([filename, ".geo"], "w");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   l = 10e-3;
%!   w = 10e-3;
%!   c = 340;
%!   rho = 1.25;
%!   c1 = 680;
%!   rho1 = 2.5;
%!   k = 100;
%!   n = 150;
%!   omega = k * c;
%!   f = omega / (2 * pi);
%!   lambda = c / f;
%!   dx = lambda / n;
%!   k1 = omega / c1;
%!   lambda1 = c1 / f;
%!   dx1 = lambda1 / n;
%!   h = min([dx, dx1]);
%!   theta = atan2(w, l);
%!   theta1 = asin(sin(theta) * c1 / c);
%!   P0 = 1;
%!   T = 2 * rho1 * k * cos(theta) / (rho * k1 * cos(theta1) + rho1 * k * cos(theta));
%!   R = (rho1 * k * cos(theta) - rho * k1 * cos(theta1)) / (rho1 * k * cos(theta) + rho * k1 * cos(theta1));
%!   p0_theta = @(x_, z_, t) P0*(R*e.^(2*1j*k*cos(theta)*z_)+1).*e.^(-1j*k*cos(theta)*z_+1j*k*sin(theta)*x_+1j*omega*t);
%!   v0x_theta = @(x_, z_, t) -(P0*R*k*sin(2*theta).*e.^(1j*k*cos(theta)*z_+1j*k*sin(theta)*x_+1j*omega*t))/(omega*rho);
%!   v0z_theta = @(x_, z_, t) -(P0*k*(R*cos(2*theta).*e.^(2*1j*k*cos(theta)*z_)-1).*e.^(-1j*k*cos(theta)*z_+1j*k*sin(theta)*x_+1j*omega*t))/(omega*rho);
%!   z0x_theta = @(x_, z_) (c*rho*csc(2*theta)*(1j*sin((2*omega*cos(theta)*z_)/c)-cos((2*omega*cos(theta)*z_)/c)-R))/R;
%!   z0z_theta = @(x_, z_) -(c*rho*(R*e.^((2*1j*omega*cos(theta)*z_)/c)+1))./(R*cos(2*theta).*e.^((2*1j*omega*cos(theta)*z_)/c)-1);
%!   p1_theta = @(x_, z_, t) P0*T*e.^(-1j*k1*cos(theta1)*z_+1j*k1*sin(theta1)*x_+1j*omega*t);
%!   v1x_theta = @(x_, z_, t) -(P0*T*k1*sin(theta1-theta).*e.^(-1j*k1*cos(theta1)*z_+1j*k1*sin(theta1)*x_+1j*omega*t))/(omega*rho1);
%!   v1z_theta = @(x_, z_, t) (P0*T*k1*cos(theta1-theta).*e.^(-1j*k1*cos(theta1)*z_+1j*k1*sin(theta1)*x_+1j*omega*t))/(omega*rho1);
%!   z1x_theta = @(x_, z_) -c1*rho1*csc(theta1-theta);
%!   z1z_theta = @(x_, z_) c1*rho1*sec(theta1-theta);
%!   p_theta = @(x_, z_, t) p0_theta(x_, z_, t) .* (z_ <= 0) + p1_theta(x_, z_, t) .* (z_ > 0);
%!   vx_theta = @(x_, z_, t) v0x_theta(x_, z_, t) .* (z_ <= 0) + v1x_theta(x_, z_, t) .* (z_ > 0);
%!   vz_theta = @(x_, z_, t) v0z_theta(x_, z_, t) .* (z_ <= 0) + v1z_theta(x_, z_, t) .* (z_ > 0);
%!   zx_theta = @(x_, z_) z0x_theta(x_, z_) .* (z_ <= 0) + z1x_theta(x_, z_) .* (z_ > 0);
%!   zz_theta = @(x_, z_) z0z_theta(x_, z_) .* (z_ <= 0) + z1z_theta(x_, z_) .* (z_ > 0);
%!   p = @(x, z, t) p_theta(x * cos(theta) - z * sin(theta), x * sin(theta) + z * cos(theta), t);
%!   vx = @(x, z, t) vx_theta(x * cos(theta) - z * sin(theta), x * sin(theta) + z * cos(theta), t);
%!   vz = @(x, z, t) vz_theta(x * cos(theta) - z * sin(theta), x * sin(theta) + z * cos(theta), t);
%!   zx = @(x, z, t) zx_theta(x * cos(theta) - z * sin(theta), x * sin(theta) + z * cos(theta));
%!   zz = @(x, z, t) zz_theta(x * cos(theta) - z * sin(theta), x * sin(theta) + z * cos(theta));
%!
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "l=%g;\n", l);
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "h=%g;\n", h);
%!   fprintf(fd, "dx = %g;\n", dx);
%!   fprintf(fd, "dx1 = %g;\n", dx1);
%!   fputs(fd, "Point(1) = {0.5 * l, 0, 0.5 * w, h};\n");
%!   fputs(fd, "Point(2) = {-0.5 * l, 0, 0.5 * w, h};\n");
%!   fputs(fd, "Point(3) = {-0.5 * l, 0, -0.5 * w, h};\n");
%!   fputs(fd, "Point(4) = {0.5 * l, 0, -0.5 * w, h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,4};\n");
%!   fputs(fd, "Line(3) = {4,1};\n");
%!   fputs(fd, "Line(4) = {2,3};\n");
%!   fputs(fd, "Line(5) = {3,4};\n");
%!   fputs(fd, "Line(6) = {4,2};\n");
%!   fputs(fd, "Line Loop(7) = {1,2,3};\n");
%!   fputs(fd, "Line Loop(8) = {4,5,6};\n");
%!   fputs(fd, "Plane Surface(9) = {7};\n");
%!   fputs(fd, "Plane Surface(10) = {8};\n");
%!   fputs(fd, "tmp[] = Extrude {0,h,0} {\n");
%!   fputs(fd, "  Surface{9};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "tmp1[] = Extrude {0,h,0} {\n");
%!   fputs(fd, "  Surface{10};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "v1 = tmp[1];\n");
%!   fputs(fd, "v2 = tmp1[1];\n");
%!   fputs(fd, "s1 = tmp[2];\n");
%!   fputs(fd, "s2 = tmp[4];\n");
%!   fputs(fd, "s3 = tmp1[2];\n");
%!   fputs(fd, "s4 = tmp1[3];\n");
%!   fputs(fd, "Geometry.OCCBooleanPreserveNumbering = 1;\n");
%!   fputs(fd, "tmp2 = BooleanFragments {Volume{tmp[1]}; Delete; }{Volume{tmp1[1]}; Delete;};\n");
%!   fputs(fd, "Physical Volume(\"top\",1) = {v1};\n");
%!   fputs(fd, "Physical Volume(\"bottom\",2) = {v2};\n");
%!   fputs(fd, "bnd1[] = Unique(Abs(Boundary{Volume{v1};}));\n");
%!   fputs(fd, "bnd2[] = Unique(Abs(Boundary{Volume{v2};}));\n");
%!   fputs(fd, "For i In {0:#bnd1[] - 1}\n");
%!   fputs(fd, "  Physical Surface(i) = {bnd1[i]};\n");
%!   fputs(fd, "EndFor\n");
%!   fputs(fd, "For i In {0:#bnd2[] - 1}\n");
%!   fputs(fd, "  Physical Surface(i+#bnd1[]) = {bnd2[i]};\n");
%!   fputs(fd, "EndFor\n");
%!   fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!   fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!   fputs(fd, "Mesh.HighOrderOptimize=2;\n");
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
%!   opt_mesh.elem_type = {"tet10", "tria6"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   empty_cell = cell(1, 4);
%!   group_defs = struct("id", empty_cell, ...
%!                       "name", empty_cell, ...
%!                       "R", empty_cell, ...
%!                       "X0", empty_cell, ...
%!                       "type", empty_cell, ...
%!                       "geometry", empty_cell);
%!   group_defs(1).id = 1;
%!   group_defs(1).name = "top";
%!   group_defs(1).R = eye(3);
%!   group_defs(1).X0 = zeros(3, 1);
%!   group_defs(1).type = "box";
%!   group_defs(1).geometry.xmin = -inf;
%!   group_defs(1).geometry.xmax = inf;
%!   group_defs(1).geometry.ymin = -inf;
%!   group_defs(1).geometry.ymax = inf;
%!   group_defs(1).geometry.zmin = 0.5 * w;
%!   group_defs(1).geometry.zmax = 0.5 * w;
%!   group_defs(2).id = 2;
%!   group_defs(2).name = "bottom";
%!   group_defs(2).R = eye(3);
%!   group_defs(2).X0 = zeros(3, 1);
%!   group_defs(2).type = "box";
%!   group_defs(2).geometry.xmin = -inf;
%!   group_defs(2).geometry.xmax = inf;
%!   group_defs(2).geometry.ymin = -inf;
%!   group_defs(2).geometry.ymax = inf;
%!   group_defs(2).geometry.zmin = -0.5 * w;
%!   group_defs(2).geometry.zmax = -0.5 * w;
%!   group_defs(3).id = 3;
%!   group_defs(3).name = "left";
%!   group_defs(3).R = eye(3);
%!   group_defs(3).X0 = zeros(3, 1);
%!   group_defs(3).type = "box";
%!   group_defs(3).geometry.xmin = -0.5 * l;
%!   group_defs(3).geometry.xmax = -0.5 * l;
%!   group_defs(3).geometry.ymin = -inf;
%!   group_defs(3).geometry.ymax = inf;
%!   group_defs(3).geometry.zmin = -inf;
%!   group_defs(3).geometry.zmax = inf;
%!   group_defs(4).id = 4;
%!   group_defs(4).name = "right";
%!   group_defs(4).R = eye(3);
%!   group_defs(4).X0 = zeros(3, 1);
%!   group_defs(4).type = "box";
%!   group_defs(4).geometry.xmin = 0.5 * l;
%!   group_defs(4).geometry.xmax = 0.5 * l;
%!   group_defs(4).geometry.ymin = -inf;
%!   group_defs(4).geometry.ymax = inf;
%!   group_defs(4).geometry.zmin = -inf;
%!   group_defs(4).geometry.zmax = inf;
%!   tol_rel = 0;
%!   tol_abs = 1e-6;
%!   x_i = mesh.nodes(:, 1);
%!   z_i = mesh.nodes(:, 3);
%!   mesh.groups.tria6 = fem_pre_mesh_groups_create(mesh, group_defs, tol_rel, tol_abs).tria6;
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(1).elements) = 2;
%!   mesh.materials.tet10(mesh.groups.tet10(2).elements) = 1;
%!   mesh.material_data = struct("rho", {rho, rho1}, "c", {c, c1});
%!   mesh.elements.particle_velocity.tria6.nodes = mesh.elements.tria6(mesh.groups.tria6(2).elements, :);
%!   mesh.materials.particle_velocity.tria6 = repmat(int32(1), rows(mesh.elements.particle_velocity.tria6.nodes), 1);
%!   mesh.elements.acoustic_impedance.tria6.nodes = mesh.elements.tria6([[mesh.groups.tria6([1,3,4])].elements], :);
%!   x_1 = x_i(mesh.elements.tria6(mesh.groups.tria6(1).elements, :));
%!   z_1 = z_i(mesh.elements.tria6(mesh.groups.tria6(1).elements, :));
%!   x_3 = x_i(mesh.elements.tria6(mesh.groups.tria6(3).elements, :));
%!   z_3 = z_i(mesh.elements.tria6(mesh.groups.tria6(3).elements, :));
%!   x_4 = x_i(mesh.elements.tria6(mesh.groups.tria6(4).elements, :));
%!   z_4 = z_i(mesh.elements.tria6(mesh.groups.tria6(4).elements, :));
%!   mesh.elements.acoustic_impedance.tria6.z = [zz(x_1, z_1);
%!                                               -zx(x_3, z_3);
%!                                               zx(x_4, z_4)];
%!   mesh.materials.acoustic_impedance.tria6 = [repmat(int32(2), numel(mesh.groups.tria6(1).elements), 1);
%!                                              repmat(int32(1), numel(mesh.groups.tria6(3).elements), 1)
%!                                              repmat(int32(2), numel(mesh.groups.tria6(4).elements), 1)];
%!   load_case = struct("particle_velocity", cell(1, 2));
%!   vz_inlet = ones(numel(mesh.groups.tria6(1).elements), 6);
%!   load_case(1).particle_velocity.tria6.vn = -real(vz(x_i(mesh.elements.particle_velocity.tria6.nodes), z_i(mesh.elements.particle_velocity.tria6.nodes), 0));
%!   load_case(2).particle_velocity.tria6.vn = -imag(vz(x_i(mesh.elements.particle_velocity.tria6.nodes), z_i(mesh.elements.particle_velocity.tria6.nodes), 0));
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case_dof.domain = FEM_DO_ACOUSTICS;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.Da_re, ...
%!    mat_ass.Da_im, ...
%!    mat_ass.Ra, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_RE, ...
%!                                         FEM_MAT_DAMPING_ACOUSTICS_IM, ...
%!                                         FEM_VEC_LOAD_ACOUSTICS], ...
%!                                        load_case);
%!   opt_sol.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_sol.solver = "pastix";
%!   opt_sol.refine_max_iter = int32(50);
%!   Keff = -omega^2 * mat_ass.Ma + 1j * omega * complex(mat_ass.Da_re, mat_ass.Da_im) + mat_ass.Ka;
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 19);
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.Phi = real(Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.PhiP = real(1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.t = Psi / omega;
%!   [sol.particle_velocity] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_VEC_PARTICLE_VELOCITY], ...
%!                                            load_case, ...
%!                                            sol);
%!   sol2.t = Psi / omega;
%!   sol2.p = real(p(x_i, z_i, sol2.t));
%!   sol2.Phi = real(1 / (-1j * omega) * p(x_i, z_i, sol2.t));
%!   sol2.PhiP = real(-p(x_i, z_i, sol2.t));
%!   sol2.particle_velocity.v.tet10 = zeros(rows(mesh.elements.tet10), columns(mesh.elements.tet10), 3, numel(Psi));
%!   for i=1:columns(mesh.elements.tet10)
%!     sol2.particle_velocity.v.tet10(:, i, 1, :) = real(vx(x_i(mesh.elements.tet10(:, i)), z_i(mesh.elements.tet10(:, i)), sol2.t));
%!     sol2.particle_velocity.v.tet10(:, i, 3, :) = real(vz(x_i(mesh.elements.tet10(:, i)), z_i(mesh.elements.tet10(:, i)), sol2.t));
%!   endfor
%!   if (do_plot)
%!   x = linspace(-0.5 * l, 0.5 * l, 100);
%!   z = linspace(-0.5 * w, 0.5 * w, 100).';
%!   gridp = linspace(min(min(sol2.p)), max(max(sol2.p)), 25);
%!   for i=1:numel(sol.t)
%!     p_i = griddata(x_i, z_i, sol.p(:, i), x, z);
%!     p2_i = griddata(x_i, z_i, sol2.p(:, i), x, z);
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     contourf(x, z, p2_i, gridp);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("reference pressure");
%!     subplot(2, 1, 2);
%!     contourf(x, z, p_i, gridp);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("FEM solution");
%!   endfor
%!   x_i = x_i(mesh.elements.tet10(:));
%!   z_i = z_i(mesh.elements.tet10(:));
%!   gridvx = linspace(min(min(min(min(sol2.particle_velocity.v.tet10)))), max(max(max(max(sol2.particle_velocity.v.tet10)))));
%!   for i=1:numel(sol.t)
%!     vx_i = griddata(x_i, z_i, sol.particle_velocity.v.tet10(:, :, 1, i)(:), x, z);
%!     v2x_i = griddata(x_i, z_i, sol2.particle_velocity.v.tet10(:, :, 1, i)(:), x, z);
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     contourf(x, z, v2x_i, gridvx);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("reference particle velocity x");
%!     subplot(2, 1, 2);
%!     contourf(x, z, vx_i, gridvx);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("FEM solution vx");
%!   endfor
%!   gridvz = linspace(min(min(min(min(sol2.particle_velocity.v.tet10)))), max(max(max(max(sol2.particle_velocity.v.tet10)))));
%!   for i=1:numel(sol.t)
%!     vz_i = griddata(x_i, z_i, sol.particle_velocity.v.tet10(:, :, 3, i)(:), x, z);
%!     v2z_i = griddata(x_i, z_i, sol2.particle_velocity.v.tet10(:, :, 3, i)(:), x, z);
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     contourf(x, z, v2z_i, gridvz);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("reference particle velocity z");
%!     subplot(2, 1, 2);
%!     contourf(x, z, vz_i, gridvz);
%!     colormap jet;
%!     colorbar;
%!     xlabel("x [m]");
%!     ylabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("FEM solution vz");
%!   endfor
%!   endif
%!   tol = 5e-3;
%!   assert_simple(sol.p, sol2.p, tol * max(max(abs(sol2.p))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
