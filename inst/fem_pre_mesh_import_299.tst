## fem_pre_mesh_import.m:299
%!test
%! ### TEST 299 - spherical waves
%! ####################################################
%! ## Jont Allen
%! ## THE ACOUSTIC WAVE EQUATION AND SIMPLE SOLUTIONS
%! ## Chapter 5
%! ####################################################
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
%!     c = 340;
%!     rho = 1.25;
%!     k = 30;
%!     omega = k * c;
%!     f = omega / (2 * pi);
%!     lambda = c / f;
%!     dx = lambda / 100;
%!     R = 100e-3;
%!     r = 20e-3;
%!     A = 100 + 30j;
%!     Phi1 = -1 * pi / 180;
%!     Phi2 = 1 * pi / 180;
%!     Phi3 = 2 * pi / 180;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R=%g;\n", R);
%!     fprintf(fd, "r=%g;\n", r);
%!     fprintf(fd, "Phi1=%g;\n", Phi1);
%!     fprintf(fd, "Phi2=%g;\n", Phi2);
%!     fprintf(fd, "Phi3=%g;\n", Phi3);
%!     fprintf(fd, "dx = %g;\n", dx);
%!     fputs(fd, "Phi4 = 0.5*(Phi1+Phi2)+Pi/2;\n");
%!     fputs(fd, "Point(1) = {0,0,0,dx};\n");
%!     fputs(fd, "Point(2) = {r * Cos(Phi1),r * Sin(Phi1),0,dx};\n");
%!     fputs(fd, "Point(3) = {R * Cos(Phi1),R * Sin(Phi1),0,dx};\n");
%!     fputs(fd, "Point(4) = {R * Cos(Phi2),R * Sin(Phi2),0,dx};\n");
%!     fputs(fd, "Point(5) = {r * Cos(Phi2),r * Sin(Phi2),0,dx};\n");
%!     fputs(fd, "Line(1) = {2,3};\n");
%!     fputs(fd, "Circle(2) = {3,1,4};\n");
%!     fputs(fd, "Line(3) = {4,5};\n");
%!     fputs(fd, "Circle(4) = {5,1,2};\n");
%!     fputs(fd, "Curve Loop(5) = {1,2,3, 4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {{Cos(Phi4),Sin(Phi4),0},{0,0,0},Phi3} {\n");
%!     fputs(fd, "  Surface{6};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"inlet\",2) = {tmp[5]};\n");
%!     fputs(fd, "Physical Surface(\"outlet\",3) = {tmp[3]};\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=1;\n");
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
%!   opt_mesh.elem_type = {"tria10", "tet20"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.c = c;
%!   grp_idx_inlet = find([mesh.groups.tria10.id] == 2);
%!   grp_idx_outlet = find([mesh.groups.tria10.id] == 3);
%!   node_idx_inlet = [[mesh.groups.tria10(grp_idx_inlet)].nodes];
%!   node_idx_outlet = [[mesh.groups.tria10(grp_idx_outlet)].nodes];
%!   node_idx_constr = [node_idx_inlet];
%!   xi = mesh.nodes(node_idx_constr, 1).';
%!   yi = mesh.nodes(node_idx_constr, 2).';
%!   zi = mesh.nodes(node_idx_constr, 3).';
%!   ri = sqrt(xi.^2 + yi.^2 + zi.^2);
%!   p = @(r, t) A * exp(1j * (omega * t - k * r)) ./ r; ## according equation 5.11.6
%!   z = @(r) rho * c * k * r .* exp(1j * acot(k * r)) ./ sqrt(1 + (k * r).^2); ## according equation 5.11.9
%!   p_constr = p(ri, 0);
%!   mesh.elements.acoustic_impedance.tria10.nodes = mesh.elements.tria10(mesh.groups.tria10(grp_idx_outlet).elements, :);
%!   xo = mesh.nodes(:, 1)(mesh.elements.acoustic_impedance.tria10.nodes);
%!   yo = mesh.nodes(:, 2)(mesh.elements.acoustic_impedance.tria10.nodes);
%!   zo = mesh.nodes(:, 3)(mesh.elements.acoustic_impedance.tria10.nodes);
%!   ro = sqrt(xo.^2 + yo.^2 + zo.^2);
%!   mesh.elements.acoustic_impedance.tria10.z = z(ro);
%!   mesh.materials.acoustic_impedance.tria10 = ones(rows(mesh.elements.acoustic_impedance.tria10.nodes), 1, "int32");
%!   mesh.elements.acoustic_constr = struct("C", mat2cell(ones(1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))), ...
%!                                          "nodes", mat2cell(node_idx_constr, 1, ones(1, numel(node_idx_constr))), ...
%!                                          "scale", mat2cell(repmat(1/omega, 1, numel(node_idx_constr)), 1, ones(1, numel(node_idx_constr))));
%!   load_case = struct("acoustic_constr", cell(1, 2));
%!   load_case(1).acoustic_constr = struct("p", mat2cell(real(p_constr), 1, ones(1, numel(p_constr))));
%!   load_case(2).acoustic_constr = struct("p", mat2cell(imag(p_constr), 1, ones(1, numel(p_constr))));
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
%!   Keff = -omega^2 * mat_ass.Ma + 1j * omega * complex(mat_ass.Da_re,  mat_ass.Da_im) + mat_ass.Ka;
%!   Reff = complex(mat_ass.Ra(:, 1), mat_ass.Ra(:, 2));
%!   Phi = fem_sol_factor(Keff, opt_sol) \ Reff;
%!   Psi = linspace(0, 2 * pi, 37);
%!   x = mesh.nodes(:, 1);
%!   sol.p = real(-1j * omega * Phi(dof_map.ndof, :) * exp(1j * Psi));
%!   sol.t = Psi / omega;
%!   xi = mesh.nodes(:, 1);
%!   yi = mesh.nodes(:, 2);
%!   zi = mesh.nodes(:, 3);
%!   ri = sqrt(xi.^2 + yi.^2 + zi.^2);
%!   sol2.p = real(p(ri, sol.t));
%!   sol2.t = sol.t;
%!   tol = 1e-4;
%!   assert_simple(sol.p, sol2.p, tol * max(max(abs(sol2.p))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
