## fem_pre_mesh_import.m:258
%!test
%! ### TEST258
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! do_plot = false;
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
%!     a = 30e-3;
%!     b = 20e-3;
%!     c = 10e-3;
%!     h = 3.5e-3;
%!     p = 25e6;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "a=%g;\n", a);
%!     fprintf(fd, "b=%g;\n", b);
%!     fprintf(fd, "c=%g;\n", c);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!     fputs(fd, "  Surface{6};\n");
%!     fputs(fd, "};\n");
%!     if (do_rotate)
%!       fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!       fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!       fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     endif
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   if (~do_rotate)
%!     group_defs(1).id = 1;
%!     group_defs(1).name = "box1";
%!     group_defs(1).R = eye(3);
%!     group_defs(1).X0 = zeros(3, 1);
%!     group_defs(1).type = "box";
%!     group_defs(1).geometry.xmin = 0;
%!     group_defs(1).geometry.xmax = 0;
%!     group_defs(1).geometry.ymin = 0;
%!     group_defs(1).geometry.ymax = b;
%!     group_defs(1).geometry.zmin = 0;
%!     group_defs(1).geometry.zmax = c;
%!     group_defs(1).elem_type = "tria10";
%!     group_defs(2).id = 2;
%!     group_defs(2).name = "cylinder1";
%!     group_defs(2).R = [-1, 0, 0;
%!                        0, 0, 1;
%!                        0, 1, 0];
%!     group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%!     group_defs(2).type = "cylinder";
%!     group_defs(2).geometry.rmin = 0;
%!     group_defs(2).geometry.rmax = 0.5 * c;
%!     group_defs(2).geometry.zmin = -0.5 * b;
%!     group_defs(2).geometry.zmax = 0.5 * b;
%!     group_defs(2).elem_type = "tria10";
%!     group_defs(3).id = 3;
%!     group_defs(3).name = "cylinder2";
%!     group_defs(3).R = [-1, 0, 0;
%!                        0, 0, 1;
%!                        0, 1, 0];
%!     group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%!     group_defs(3).type = "cylinder";
%!     group_defs(3).geometry.rmin = 0;
%!     group_defs(3).geometry.rmax = 0.5 * c;
%!     group_defs(3).geometry.zmin = -0.5 * b;
%!     group_defs(3).geometry.zmax = 0.5 * b;
%!     group_defs(3).elem_type = "tria10";
%!     group_defs(4).id = 4;
%!     group_defs(4).name = "box2";
%!     group_defs(4).R = eye(3);
%!     group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%!     group_defs(4).type = "box";
%!     group_defs(4).geometry.xmin = 0;
%!     group_defs(4).geometry.xmax = 0;
%!     group_defs(4).geometry.ymin = -0.5 * b;
%!     group_defs(4).geometry.ymax = 0.5 * b;
%!     group_defs(4).geometry.zmin = -0.5 * c;
%!     group_defs(4).geometry.zmax = 0.5 * c;
%!     group_defs(4).elem_type = "tria10";
%!     groups = fem_pre_mesh_groups_create(mesh, group_defs, sqrt(eps));
%!     assert_simple(numel(groups.tria10), 4);
%!     assert_simple([groups.tria10.id], [group_defs.id]);
%!     assert_simple(sort(groups.tria10(1).nodes), sort(mesh.groups.tria10(1).nodes));
%!     assert_simple(sort(groups.tria10(2).nodes), sort(mesh.groups.tria10(1).nodes));
%!     assert_simple(sort(groups.tria10(3).nodes), sort(mesh.groups.tria10(2).nodes));
%!     assert_simple(sort(groups.tria10(4).nodes), sort(mesh.groups.tria10(2).nodes));
%!     assert_simple(groups.tria10(1).elements, mesh.groups.tria10(1).elements);
%!     assert_simple(groups.tria10(2).elements, mesh.groups.tria10(1).elements);
%!     assert_simple(groups.tria10(3).elements, mesh.groups.tria10(2).elements);
%!     assert_simple(groups.tria10(4).elements, mesh.groups.tria10(2).elements);
%!   endif
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria10(find([[mesh.groups.tria10].id] == 1)).nodes, :) = true;
%!   load_case.pressure.tria10.elements = mesh.elements.tria10(mesh.groups.tria10(find([mesh.groups.tria10.id] == 2)).elements, :);
%!   Xp = mesh.nodes(load_case.pressure.tria10.elements, 1:3);
%!   xp = reshape(Xp(:, 1), rows(load_case.pressure.tria10.elements), columns(load_case.pressure.tria10.elements));
%!   yp = reshape(Xp(:, 2), rows(load_case.pressure.tria10.elements), columns(load_case.pressure.tria10.elements));
%!   zp = reshape(Xp(:, 3), rows(load_case.pressure.tria10.elements), columns(load_case.pressure.tria10.elements));
%!   load_case.pressure.tria10.p = p / 2 * (yp / b + zp / c);
%!   mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.Mlumped, ...
%!    mat_ass.R, ...
%!    mat_ass.colloc_mass, ...
%!    mat_ass.colloc_stiffness, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_MAT_MASS_LUMPED, ...
%!                            FEM_VEC_LOAD_CONSISTENT, ...
%!                            FEM_VEC_COLL_MASS, ...
%!                            FEM_VEC_COLL_STIFFNESS, ...
%!                            FEM_SCA_TOT_MASS], load_case);
%!   assert_simple(mtot, a * b * c * mesh.material_data.rho, sqrt(eps) * a * b * c * mesh.material_data.rho);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   [sol_eig_lumped] = fem_sol_modal(mesh, dof_map, setfield(mat_ass, "M", mat_ass.Mlumped), number_of_modes);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   X = mesh.nodes(unique(load_case.pressure.tria10.elements), 1:3).';
%!   dof_idx = dof_map.ndof(unique(load_case.pressure.tria10.elements), 1:3);
%!   F_con = full(mat_ass.R(dof_idx)).';
%!   M_con = cross(X, F_con);
%!   Ftot_con = sum(F_con, 2);
%!   Mtot_con = sum(M_con, 2);
%!   Fx_an = -(b * c * p) / 2; ## grind(integrate(integrate(-1/2*p*(y/b+z/c),z,0,c),y,0,b));
%!   Mz_an = (7 * b^2 * c * p) / 24; ## grind(integrate(integrate(1/2*p*(y/b+z/c)*y,z,0,c),y,0,b));
%!   My_an = -(7 * b * c^2 * p) / 24; ## grind(integrate(integrate(-1/2*p*(y/b+z/c)*z,z,0,c),y,0,b));
%!   F_an = [Fx_an; 0; 0];
%!   M_an = [0; My_an; Mz_an];
%!   assert_simple(Ftot_con, F_an, eps^0.9 * norm(F_an));
%!   assert_simple(Mtot_con, M_an, eps^0.9 * norm(M_an));
%!   f = sol_eig.f(:);
%!   f_lumped = sol_eig_lumped.f(:);
%!   f_ref = [8768.74;
%!            14636.1;
%!            21145.7;
%!            39712.8;
%!            43555.5;
%!            47909;
%!            62270.4;
%!            84324.4;
%!            92665.1;
%!            94563];
%!   for i=1:length(f)
%!     fprintf(stderr, "mode %d f=%.0f f_lumped=%.0f\n", i, f(i), f_lumped(i));
%!   endfor
%!   assert_simple(all(f_lumped <= f));
%!   assert_simple(f, f_ref, tol * max(f_ref));
%!   assert_simple(f_lumped, f_ref, 5 * tol * max(f_ref));
%!   if (do_plot)
%!     figure("visible","off");
%!     hold on;
%!     fem_post_sol_plot(mesh);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title('undeformed mesh');
%!     opts_plot.elem_types = {"tria10", "tet20"};
%!     opts_plot.elem_groups.tria10 = [mesh.groups.tria10.id];
%!     opts_plot.elem_groups.tet20 = [mesh.groups.tet20.id];
%!     for i=1:min(number_of_modes, length(sol_eig.f))
%!       figure("visible", "off");
%!       hold on;
%!       fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, :, i), "rows")), i, opts_plot);
%!       view(30,30);
%!       xlabel('x [m]');
%!       ylabel('y [m]');
%!       zlabel('z [m]');
%!       grid on;
%!       grid minor on;
%!       title(sprintf("%d. eigenmode: %gHz",i,sol_eig.f(i)));
%!     endfor
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_stat, scale_eig/max(norm(sol_stat.def, "rows")), 1, opts_plot);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title("deformed mesh consistent load vector");
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_stat_lumped, scale_eig/max(norm(sol_stat_lumped.def, "rows")), 1, opts_plot);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title("deformed mesh lumped load vector");
%!     figure_list();
%!   endif
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
