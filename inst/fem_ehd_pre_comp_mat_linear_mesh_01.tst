%!test
%! try
%!   ## TEST 1
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 3e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(20);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       Psi = 0.8e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       mesh_size = 5e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.modal.name = "node_id_conrod_big_end";
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     cms_opt.nodes.interfaces.name = "node_id_conrod_small_end";
%!     cms_opt.element.name = "elem_id_conrod_modal";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.load_cases = "index";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.solver = "pardiso";
%!     cms_opt.verbose = int32(0);
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + 100;
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).bearing = "elem_id_piston_pin_bearing";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.number_of_nodes_x = 40;
%!     bearing_surf(1).options.number_of_nodes_z = 10;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + 100;
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).bearing = "elem_id_crank_pin_bearing";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.number_of_nodes_x = 40;
%!     bearing_surf(2).options.number_of_nodes_z = 10;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = true;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);
%!     for i=1:numel(comp_mat)
%!       fem_ehd_pre_comp_mat_export(comp_mat(i), bearing_surf(i).options, sprintf("%s_ehd_%d.dat", filename, i));
%!     endfor
%!     for k=1:4
%!       unwind_protect
%!         fd = fopen(sprintf("%s_%d.mbdyn", filename, k), "wt");
%!         if (fd == -1)
%!           error("failed to open file %s.mbdyn", filename);
%!         endif
%!         fputs(fd, " set: integer node_id_conrod_big_end = 1001;\n");
%!         fputs(fd, " set: integer node_id_conrod_small_end = 1002;\n");
%!         fputs(fd, " set: integer node_id_piston = 1004;\n");
%!         fputs(fd, " set: integer node_id_crank_pin = 1005;\n");
%!         fputs(fd, " set: integer elem_id_conrod_modal = 2001;\n");
%!         fputs(fd, " set: integer joint_id_conrod_big_end = 2002;\n");
%!         fputs(fd, " set: integer joint_id_piston = 2003;\n");
%!         fputs(fd, " set: integer joint_id_crank_pin = 2004;\n");
%!         fputs(fd, " set: integer elem_id_piston_pin_bearing = 2005;\n");
%!         fputs(fd, " set: integer elem_id_crank_pin_bearing = 2006;\n");
%!         fputs(fd, " set: integer ref_id_conrod_big_end = 3001;\n");
%!         fputs(fd, " set: integer ref_id_conrod_small_end = 3002;\n");
%!         fputs(fd, " set: integer ref_id_conrod_modal = 3003;\n");
%!         fprintf(fd, " set: real piston_pin_bearing_width = %g;\n", w);
%!         fprintf(fd, " set: real piston_pin_bearing_journal_diameter = %g;\n", d * (1 - Psi));
%!         fprintf(fd, " set: real piston_pin_bearing_shell_diameter = %g;\n", d);
%!         fprintf(fd, " set: real crank_pin_bearing_width = %g;\n", w);
%!         fprintf(fd, " set: real crank_pin_bearing_journal_diameter = %g;\n", d * (1 - Psi));
%!         fprintf(fd, " set: real crank_pin_bearing_shell_diameter = %g;\n", d);
%!         fputs(fd, " set: real rho = 800.;\n");
%!         fputs(fd, " set: real eta = 5e-3;\n");
%!         fprintf(fd, "set: real conrod_length = %g;\n", l);
%!         fprintf(fd, "set: integer piston_pin_bearing_shell_number_of_nodes_x = %d;\n", numel(comp_mat(1).bearing_surf.grid_x) + 1);
%!         fprintf(fd, "set: integer piston_pin_bearing_shell_number_of_nodes_z = %d;\n", numel(comp_mat(1).bearing_surf.grid_z));
%!         fprintf(fd, "set: integer crank_pin_bearing_shell_number_of_nodes_x = %d;\n", numel(comp_mat(2).bearing_surf.grid_x) + 1);
%!         fprintf(fd, "set: integer crank_pin_bearing_shell_number_of_nodes_z = %d;\n", numel(comp_mat(2).bearing_surf.grid_z));
%!         fputs(fd, " begin: data;\n");
%!         fputs(fd, "     problem: initial value; # the default\n");
%!         fputs(fd, " end: data;\n");
%!         fputs(fd, " begin: initial value;\n");
%!         fputs(fd, "     initial time: 0;\n");
%!         fputs(fd, "     final time: 1;\n");
%!         fputs(fd, "     time step: 0.01;\n");
%!         fputs(fd, "     max iterations: 100;\n");
%!         fputs(fd, "   linear solver: umfpack, grad, scale, row max column max, always, max iterations, 10;\n");
%!         fputs(fd, "   method: implicit euler;\n");
%!         fputs(fd, "     threads: disable;\n");
%!         fputs(fd, "   nonlinear solver: linesearch, default solver options, heavy nonlinear, divergence check, no;\n");
%!         fputs(fd, "   derivatives max iterations: 10;\n");
%!         fputs(fd, "   tolerance: 1e-6, 1e-6;\n");
%!         fputs(fd, "   output: iterations;\n");
%!         fputs(fd, " end: initial value;\n");
%!         fputs(fd, " begin: control data;\n");
%!         fputs(fd, "        use automatic differentiation;\n");
%!         fputs(fd, "        max iterations: 1;\n");
%!         fputs(fd, "        structural nodes: 4;\n");
%!         fputs(fd, "        joints: 4;\n");
%!         fputs(fd, "        loadable elements: 2;\n");
%!         fputs(fd, "        default output: all;\n");
%!         fputs(fd, "        default output: reference frames;\n");
%!         fputs(fd, "        output meter: closest next, 1., forever, const, 0.1;\n");
%!         fputs(fd, " end: control data;\n");
%!         fputs(fd, " begin: nodes;\n");
%!         fputs(fd, " reference: ref_id_conrod_big_end,\n");
%!         fputs(fd, "     position,\n");
%!         fputs(fd, "     reference, global, null,\n");
%!         fputs(fd, "     orientation,\n");
%!         fputs(fd, "     reference, global, eye,\n");
%!         fputs(fd, "     velocity,\n");
%!         fputs(fd, "     reference, global, null,\n");
%!         fputs(fd, "     angular velocity,\n");
%!         fputs(fd, "     reference, global, null;\n");
%!         fputs(fd, " reference: ref_id_conrod_small_end,\n");
%!         fputs(fd, "     position,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end,\n");
%!         fputs(fd, "                        conrod_length,\n");
%!         fputs(fd, "                        0.,\n");
%!         fputs(fd, "                        0.,\n");
%!         fputs(fd, "     orientation,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "     velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "     angular velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null;\n");
%!         fputs(fd, " structural: node_id_conrod_big_end, modal,\n");
%!         fputs(fd, "     position,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "     orientation,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "     velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "     angular velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null;\n");
%!         fputs(fd, " structural: node_id_conrod_small_end, static,\n");
%!         fputs(fd, "     position,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "     orientation,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "     velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "     angular velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, null;\n");
%!         fputs(fd, " structural: node_id_piston, static,\n");
%!         fputs(fd, "     position,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "     orientation,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "     velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "     angular velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_small_end, null;\n");
%!         fputs(fd, " structural: node_id_crank_pin, static,\n");
%!         fputs(fd, "     position,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "     orientation,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "     velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "     angular velocity,\n");
%!         fputs(fd, "     reference, ref_id_conrod_big_end, null;\n");
%!         fputs(fd, " end: nodes;\n");
%!         fputs(fd, " begin: elements;\n");
%!         fprintf(fd, " include: \"%s_cms.elm\";\n", filename);
%!         fputs(fd, "  joint: joint_id_conrod_big_end, total pin joint,\n");
%!         fputs(fd, "       node_id_conrod_big_end,\n");
%!         fputs(fd, "       position, reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "       position orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       rotation orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       position, reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "       position orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       rotation orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       position constraint,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        null,\n");
%!         fputs(fd, "       orientation constraint,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         null;\n");
%!         fputs(fd, "  joint: joint_id_piston, total pin joint,\n");
%!         fputs(fd, "       node_id_piston,\n");
%!         fputs(fd, "       position, reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "       position orientation, reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "       rotation orientation, reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "       position, reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "       position orientation, reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "       rotation orientation, reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "       position constraint,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        null,\n");
%!         fputs(fd, "       orientation constraint,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         null;\n");
%!         fputs(fd, "  joint: joint_id_crank_pin, total pin joint,\n");
%!         fputs(fd, "       node_id_crank_pin,\n");
%!         fputs(fd, "       position, reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "       position orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       rotation orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       position, reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "       position orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       rotation orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "       position constraint,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        active,\n");
%!         fputs(fd, "        null,\n");
%!         fputs(fd, "       orientation constraint,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         active,\n");
%!         fputs(fd, "         null;\n");
%!         fputs(fd, "     user defined: elem_id_piston_pin_bearing, hydrodynamic plain bearing2,\n");
%!         fputs(fd, "             hydraulic fluid, incompressible,\n");
%!         fputs(fd, "                     density, rho,\n");
%!         fputs(fd, "                     viscosity, eta,\n");
%!         fputs(fd, "                     pressure, 0,\n");
%!         fputs(fd, "             mesh, linear finite difference,\n");
%!         fputs(fd, "             geometry, cylindrical,\n");
%!         fputs(fd, "                     mesh position, at bearing,\n");
%!         fputs(fd, "                     bearing width, piston_pin_bearing_width,\n");
%!         fputs(fd, "                     shaft diameter, piston_pin_bearing_journal_diameter,\n");
%!         fputs(fd, "                     bearing diameter, piston_pin_bearing_shell_diameter,\n");
%!         fputs(fd, "             shaft node, node_id_piston, \n");
%!         fputs(fd, "             offset, reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "             orientation, reference, ref_id_conrod_small_end, eye,\n");
%!         fputs(fd, "             bearing node, node_id_conrod_big_end,\n");
%!         fputs(fd, "             offset, reference, ref_id_conrod_small_end, null,\n");
%!         fputs(fd, "             orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "                     number of nodes z, piston_pin_bearing_shell_number_of_nodes_z,\n");
%!         fputs(fd, "                     number of nodes Phi, piston_pin_bearing_shell_number_of_nodes_x,\n");
%!         fputs(fd, "                     boundary conditions, \n");
%!         fputs(fd, "                                     pressure, const, 0.,\n");
%!         fputs(fd, "                                     pressure, const, 0.,\n");
%!         switch (k)
%!           case 1
%!             fputs(fd, "                  lubrication grooves, 1,\n");
%!             fputs(fd, "                     at bearing,\n");
%!             fprintf(fd, "                     pressure, mult, time, const, %g,\n", p1 * bearing_surf(1).options.reference_pressure);
%!             fputs(fd, "                       position, 0., 0.,\n");
%!             fputs(fd, "                       complete surface,\n");
%!           case 2
%!             fprintf(fd, "                       lubrication grooves, %d,\n", (numel(comp_mat(1).bearing_surf.grid_x) - 1) * numel(comp_mat(1).bearing_surf.grid_z));
%!             for i=1:numel(comp_mat(1).bearing_surf.grid_x) - 1
%!               for j=1:numel(comp_mat(1).bearing_surf.grid_z)
%!                 fputs(fd, "                             at bearing,\n");
%!                 fprintf(fd, "                                   pressure, mult, time, const, %g,\n", p1 * (sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).grid_z(j) / w + 1) * bearing_surf(1).options.reference_pressure);
%!                 fprintf(fd, "                                   position, %g, %g,\n", comp_mat(1).bearing_surf.grid_x(i), comp_mat(1).bearing_surf.grid_z(j));
%!                 fprintf(fd, "                                   rectangle, width, %g, height, %g,\n", comp_mat(1).bearing_surf.grid_x(2) - comp_mat(1).bearing_surf.grid_x(1), comp_mat(1).bearing_surf.grid_z(2) - comp_mat(1).bearing_surf.grid_z(1));
%!               endfor
%!             endfor
%!         otherwise
%!           fputs(fd, "                            lubrication grooves, 1,\n");
%!           fputs(fd, "                              at bearing,\n");
%!           fputs(fd, "                              pressure, const, 0.,\n");
%!           fputs(fd, "                              position, 0., 0.,\n");
%!           fputs(fd, "                              complete surface,\n");
%!         endswitch
%!         fputs(fd, "             compliance model,\n");
%!         fprintf(fd, "                   matrix, 1, from file, \"%s_ehd_1.dat\",\n", filename);
%!         fprintf(fd, "                   E1, %g,\n", mesh.material_data(1).E);
%!         fprintf(fd, "                   nu1, %g,\n", mesh.material_data(1).nu);
%!         fputs(fd, "                         modal element, elem_id_conrod_modal,\n");
%!         fputs(fd, "                 output total deformation, yes,\n");
%!         fputs(fd, "                 output pressure, yes,\n");
%!         fputs(fd, "                 output, yes;\n");
%!         fputs(fd, "     user defined: elem_id_crank_pin_bearing, hydrodynamic plain bearing2,\n");
%!         fputs(fd, "             hydraulic fluid, incompressible,\n");
%!         fputs(fd, "                     density, rho,\n");
%!         fputs(fd, "                     viscosity, eta,\n");
%!         fputs(fd, "                     pressure, 0,\n");
%!         fputs(fd, "             mesh, linear finite difference,\n");
%!         fputs(fd, "             geometry, cylindrical,\n");
%!         fputs(fd, "                     mesh position, at bearing,\n");
%!         fputs(fd, "                     bearing width, crank_pin_bearing_width,\n");
%!         fputs(fd, "                     shaft diameter, crank_pin_bearing_journal_diameter,\n");
%!         fputs(fd, "                     bearing diameter, crank_pin_bearing_shell_diameter,\n");
%!         fputs(fd, "             shaft node, node_id_crank_pin, \n");
%!         fputs(fd, "             offset, reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "             orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "             bearing node, node_id_conrod_big_end,\n");
%!         fputs(fd, "             offset, reference, ref_id_conrod_big_end, null,\n");
%!         fputs(fd, "             orientation, reference, ref_id_conrod_big_end, eye,\n");
%!         fputs(fd, "                     number of nodes z, crank_pin_bearing_shell_number_of_nodes_z,\n");
%!         fputs(fd, "                     number of nodes Phi, crank_pin_bearing_shell_number_of_nodes_x,\n");
%!         fputs(fd, "                     boundary conditions, \n");
%!         fputs(fd, "                                     pressure, const, 0.,\n");
%!         fputs(fd, "                                     pressure, const, 0.,\n");
%!         switch (k)
%!           case 3
%!             fputs(fd, "                  lubrication grooves, 1,\n");
%!             fputs(fd, "                     at bearing,\n");
%!             fprintf(fd, "                     pressure, mult, time, const, %g,\n", p2 * bearing_surf(2).options.reference_pressure);
%!             fputs(fd, "                       position, 0., 0.,\n");
%!             fputs(fd, "                       complete surface,\n");
%!           case 4
%!             fprintf(fd, "                       lubrication grooves, %d,\n", (numel(comp_mat(2).bearing_surf.grid_x) - 1) * numel(comp_mat(2).bearing_surf.grid_z));
%!             for i=1:numel(comp_mat(2).bearing_surf.grid_x) - 1
%!               for j=1:numel(comp_mat(2).bearing_surf.grid_z)
%!                 fputs(fd, "                             at bearing,\n");
%!                 fprintf(fd, "                                   pressure, mult, time, const, %g,\n", p2 * (cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r).^2 * bearing_surf(2).grid_z(j) / w + 1) * bearing_surf(2).options.reference_pressure);
%!                 fprintf(fd, "                                   position, %g, %g,\n", comp_mat(2).bearing_surf.grid_x(i), comp_mat(2).bearing_surf.grid_z(j));
%!                 fprintf(fd, "                                   rectangle, width, %g, height, %g,\n", comp_mat(2).bearing_surf.grid_x(2) - comp_mat(2).bearing_surf.grid_x(1), comp_mat(2).bearing_surf.grid_z(2) - comp_mat(2).bearing_surf.grid_z(1));
%!               endfor
%!             endfor
%!         otherwise
%!           fputs(fd, "                            lubrication grooves, 1,\n");
%!           fputs(fd, "                              at bearing,\n");
%!           fputs(fd, "                              pressure, const, 0.,\n");
%!           fputs(fd, "                              position, 0., 0.,\n");
%!           fputs(fd, "                              complete surface,\n");
%!         endswitch
%!         fputs(fd, "             compliance model,\n");
%!         fprintf(fd, "                   matrix, 1, from file, \"%s_ehd_2.dat\",\n", filename);
%!         fprintf(fd, "                   E1, %g,\n", mesh.material_data(1).E);
%!         fprintf(fd, "                   nu1, %g,\n", mesh.material_data(1).nu);
%!         fputs(fd, "                         modal element, elem_id_conrod_modal,\n");
%!         fputs(fd, "                 output total deformation, yes,\n");
%!         fputs(fd, "                 output pressure, yes,\n");
%!         fputs(fd, "                 output, yes;\n");
%!         fputs(fd, " end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!         fd = -1;
%!       end_unwind_protect
%!       opt_mbd.output_file = sprintf("%s_%d_mbdyn", filename, k);
%!       shell(sprintf("nl \"%s\"", sprintf("%s_%d.mbdyn", filename, k)));
%!       info = mbdyn_solver_run(sprintf("%s_%d.mbdyn", filename, k), opt_mbd);
%!       log_dat(k) = mbdyn_post_load_log(opt_mbd.output_file);
%!       bearings(k, :) = mbdyn_post_ehd_load_output(opt_mbd.output_file, log_dat(k));
%!       sol_mbd(k) = fem_post_cms_sol_import(opt_mbd.output_file, struct("cms_opt", cms_opt));
%!       sol_mbd_exp.def(:, :, k) = fem_post_cms_expand_body(mesh, dof_map, mat_ass, sol_mbd(k).bodies(1).q);
%!     endfor
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = 1;
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * (sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).grid_z / w + 1) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     p2red1 = repmat(p2 * bearing_surf(1).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * (cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r).^2 * bearing_surf(2).grid_z / w + 1) * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end - nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     qred = mat_ass.Kred \ Fred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements.joints = mesh_post.elements.joints(1:end - 2);
%!     load_case_post = fem_pre_load_case_create_empty(4);
%!     load_case_post_dof.locked_dof = false(size(mesh_post.nodes));
%!     x1 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     z1 = mesh.nodes(:, 3)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(3);
%!     Phi1 = atan2(y1, x1);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     z2 = mesh.nodes(:, 3)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(3);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(1).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(1).pressure.tria6.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%!     load_case_post(2).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(2).pressure.tria6.p = p1 * (sin(Phi1) .* z1 / w + 1) * bearing_surf(1).options.reference_pressure;
%!     load_case_post(3).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(3).pressure.tria6.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6);
%!     load_case_post(4).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(4).pressure.tria6.p = p2 * (cos(Phi1).^2 .* z2 / w + 1) * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints(end + 1).C = eye(6);
%!     mesh_post.elements.joints(end).nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post_dof);
%!     [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     mesh_data(3).mesh = mesh;
%!     mesh_data(3).mesh.nodes(:, 2) += 50e-3;
%!     mesh_data(3).dof_map = dof_map;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(3) + (1:rows(sol_mbd_exp.def)), :, :) = sol_mbd_exp.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(:, :, i) - sol_red.def(:, :, i)))) / max(max(abs(sol_post.def(:, :, i))));
%!     endfor
%!     printf("difference between modal and FEM:\n");
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     err_red2 = zeros(1, size(sol_post.def, 3));
%!     for i=1:size(sol_post.def, 3)
%!       err_red2(i) = max(max(abs(sol_red.def(:, :, i) - sol_mbd_exp.def(:, :, i)))) / max(max(abs(sol_red.def(:, :, i))));
%!     endfor
%!     printf("difference modal and MBDyn:\n");
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red2(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     assert_simple(all(err_red2 < tol_red));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
