## fem_pre_mesh_import.m:03
%!test
%! ## TEST3
%! number_of_modes = 3;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! plot_modes = false;
%! if (plot_modes)
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
%!     a = 30e-3;
%!     b = 20e-3;
%!     c = 10e-3;
%!     d = -5e-3;
%!     e = 35e-3;
%!     h = 4e-3;
%!     alpha = 1e-6;
%!     beta = 1e-4;
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
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"modal\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"interfaces\",2) = {tmp[2]};\n");
%!     fputs(fd, "SetFactory(\"Built-in\");\n");
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
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   cms_opt.modes.number = int32(6);
%!   cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!   cms_opt.nodes.interfaces.number = int32(rows(mesh.nodes) + 2);
%!   cms_opt.invariants = false;
%!   cms_opt.algorithm = "unsymmetric";
%!   mesh.nodes(cms_opt.nodes.modal.number, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%!   mesh.nodes([cms_opt.nodes.interfaces.number], :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];
%!   [~] = unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                    [1, 2], ...
%!                                                    [cms_opt.nodes.modal.number, ...
%!                                                     cms_opt.nodes.interfaces.number]);
%!   for i=1:numel(mesh.elements.rbe3)
%!     assert_simple(sum(mesh.elements.rbe3(i).weight), b * c, sqrt(eps) * (b * c));
%!   endfor
%!   [mesh, mat_ass_cms, dof_map_cms, sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%!   mat_ass_cms.Dred = alpha * mat_ass_cms.Mred + beta * mat_ass_cms.Kred;
%!   if (plot_modes)
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
%!     for i=1:min(number_of_modes, numel(sol_eig_cms.f))
%!       opt_plot.elem_types = {"tria6"};
%!       figure("visible", "off");
%!       hold on;
%!       fem_post_sol_plot(mesh, sol_eig_cms, scale_eig / max(norm(sol_eig_cms.def(:, 1:3, i), "rows")), i, opt_plot);
%!       view(30,30);
%!       xlabel('x [m]');
%!       ylabel('y [m]');
%!       zlabel('z [m]');
%!       grid on;
%!       grid minor on;
%!       title(sprintf("%d. eigenmode: %gHz", i, sol_eig_cms.f(i)));
%!     endfor
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
