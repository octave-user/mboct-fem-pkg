## fem_pre_mesh_import.m:373
%!test
%! try
%! ## TEST 373
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.E = 55e6 / SI_unit_pascal;
%! param.delta = 0.0;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.p1 = 0.006e6 / SI_unit_pascal;
%! param.bound_cond = "edge";
%! param.analysis = "plain strain";
%! options.verbose = false;
%! elem_types = {"tet20", "tet10", "tet10h", "iso8", "penta15", "iso20", "iso20r", "iso27", "penta18"};
%! nu_val = [0.3];
%! for idx_elem_type=1:numel(elem_types)
%! param.elem_type = elem_types{idx_elem_type};
%! switch (param.elem_type)
%! case {"iso8", "penta15", "iso20", "iso20r", "iso27", "penta18"}
%!   param.transfinite = true;
%! otherwise
%!   param.transfinite = false;
%! endswitch
%! for idx_nu = 1:numel(nu_val)
%! param.nu = nu_val(idx_nu);
%! switch (param.elem_type)
%! case {"iso8"}
%! param.h = 10e-3 / 32 / SI_unit_meter;
%! case {"tet20"}
%! param.h = 10e-3 / 8 / SI_unit_meter;
%! otherwise
%! param.h = 10e-3 / 16 / SI_unit_meter;
%! endswitch
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fprintf(fd, "SetFactory(\"Built-in\");\n");
%!     fprintf(fd, "h = %g;\n", param.h);
%!     fprintf(fd, "Point(1) = {0,0,20};\n");
%!     fprintf(fd, "Point(2) = {5,0,20};\n");
%!     fprintf(fd, "Point(3) = {10,0,20};\n");
%!     fprintf(fd, "Point(4) = {10,0,10};\n");
%!     fprintf(fd, "Point(5) = {15 - 5 * Cos(Pi/4), 0, 10 - 5 * Sin(Pi/4)};\n");
%!     fprintf(fd, "Point(6) = {15, 0, 5};\n");
%!     fprintf(fd, "Point(7) = {65, 0, 5};\n");
%!     fprintf(fd, "Point(8) = {65, 0, 0};\n");
%!     fprintf(fd, "Point(9) = {65, 0, -5};\n");
%!     fprintf(fd, "Point(10) = {15, 0, -5};\n");
%!     fprintf(fd, "Point(11) = {15 - 5 * Cos(Pi/4), 0, -(10 - 5 * Sin(Pi/4))};\n");
%!     fprintf(fd, "Point(12) = {10, 0, -10};\n");
%!     fprintf(fd, "Point(13) = {10,0,-20};\n");
%!     fprintf(fd, "Point(14) = {5, 0, -20};\n");
%!     fprintf(fd, "Point(15) = {0, 0, -20};\n");
%!     fprintf(fd, "Point(16) = {0, 0, -10};\n");
%!     fprintf(fd, "Point(17) = {0, 0, 0};\n");
%!     fprintf(fd, "Point(18) = {0, 0, 10};\n");
%!     fprintf(fd, "Point(19) = {5, 0, 10};\n");
%!     fprintf(fd, "Point(20) = {5, 0, 0};\n");
%!     fprintf(fd, "Point(21) = {15, 0, 0};\n");
%!     fprintf(fd, "Point(22) = {5, 0, -10};\n");
%!     fprintf(fd, "Point(23) = {15,0,10};\n");
%!     fprintf(fd, "Point(25) = {15,0,-10};\n");
%!     fprintf(fd, "Line(1) = {1,2};\n");
%!     fprintf(fd, "Line(2) = {2,3};\n");
%!     fprintf(fd, "Line(3) = {3,4};\n");
%!     fprintf(fd, "Circle(4) = {4,23,5};\n");
%!     fprintf(fd, "Circle(5) = {5,23,6};\n");
%!     fprintf(fd, "Line(6) = {6,7};\n");
%!     fprintf(fd, "Line(7) = {7,8};\n");
%!     fprintf(fd, "Line(8) = {8,9};\n");
%!     fprintf(fd, "Line(9) = {9,10};\n");
%!     fprintf(fd, "Circle(10)={10,25,11};\n");
%!     fprintf(fd, "Circle(11)={11,25,12};\n");
%!     fprintf(fd, "Line(12) = {12,13};\n");
%!     fprintf(fd, "Line(13) = {13,14};\n");
%!     fprintf(fd, "Line(14) = {14,15};\n");
%!     fprintf(fd, "Line(15) = {15,16};\n");
%!     fprintf(fd, "Line(16) = {16,17};\n");
%!     fprintf(fd, "Line(17) = {17,18};\n");
%!     fprintf(fd, "Line(18) = {18,1};\n");
%!     fprintf(fd, "Line(19) = {2,19};\n");
%!     fprintf(fd, "Line(20) = {19,20};\n");
%!     fprintf(fd, "Line(21) = {20,22};\n");
%!     fprintf(fd, "Line(22) = {22,14};\n");
%!     fprintf(fd, "Line(23) = {18,19};\n");
%!     fprintf(fd, "Line(24) = {19, 4};\n");
%!     fprintf(fd, "Line(25) = {17,20};\n");
%!     fprintf(fd, "Line(26) = {16,22};\n");
%!     fprintf(fd, "Line(27) = {22,12};\n");
%!     fprintf(fd, "Line(28) = {20,21};\n");
%!     fprintf(fd, "Line(29) = {21,8};\n");
%!     fprintf(fd, "Line(30) = {6,21};\n");
%!     fprintf(fd, "Line(31) = {21,10};\n");
%!     fprintf(fd, "Line(32) = {20,5};\n");
%!     fprintf(fd, "Line(33) = {20,11};\n");
%!     if (param.transfinite)
%!     fprintf(fd, "Transfinite Curve(1) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(2) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(3) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(4) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(5) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(6) = Max(1, Round(50 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(7) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(8) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(9) = Max(1, Round(50 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(10) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(11) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(12) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(13) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(14) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(15) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(16) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(17) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(18) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(19) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(20) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(21) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(22) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(23) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(24) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(25) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(26) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(27) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(28) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(29) = Max(1, Round(50 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(30) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(31) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(32) = Max(1, Round(5 / h)) + 1;\n");
%!     fprintf(fd, "Transfinite Curve(33) = Max(1, Round(5 / h)) + 1;\n");
%!     endif
%!     fprintf(fd, "Line Loop(1) = {1,19,-23,18};\n");
%!     fprintf(fd, "Line Loop(2) = {2,3,-24,-19};\n");
%!     fprintf(fd, "Line Loop(3) = {23,20,-25,17};\n");
%!     fprintf(fd, "Line Loop(4) = {24,4,-32,-20};\n");
%!     fprintf(fd, "Line Loop(5) = {32,5,30,-28};\n");
%!     fprintf(fd, "Line Loop(6) = {6,7,-29,-30};\n");
%!     fprintf(fd, "Line Loop(7) = {25,21,-26,16};\n");
%!     fprintf(fd, "Line Loop(8) = {33,11,-27,-21};\n");
%!     fprintf(fd, "Line Loop(9) = {28,31,10,-33};\n");
%!     fprintf(fd, "Line Loop(10) = {29,8,9,-31};\n");
%!     fprintf(fd, "Line Loop(11) = {26, 22, 14, 15};\n");
%!     fprintf(fd, "Line Loop(12) = {27, 12, 13, -22};\n");
%!     fprintf(fd, "Plane Surface(1) = {1};\n");
%!     fprintf(fd, "Plane Surface(2) = {2};\n");
%!     fprintf(fd, "Plane Surface(3) = {3};\n");
%!     fprintf(fd, "Plane Surface(4) = {4};\n");
%!     fprintf(fd, "Plane Surface(5) = {5};\n");
%!     fprintf(fd, "Plane Surface(6) = {6};\n");
%!     fprintf(fd, "Plane Surface(7) = {7};\n");
%!     fprintf(fd, "Plane Surface(8) = {8};\n");
%!     fprintf(fd, "Plane Surface(9) = {9};\n");
%!     fprintf(fd, "Plane Surface(10) = {10};\n");
%!     fprintf(fd, "Plane Surface(11) = {11};\n");
%!     fprintf(fd, "Plane Surface(12) = {12};\n");
%!     if (param.transfinite)
%!     fprintf(fd, "Transfinite Surface(1) = {PointsOf{Surface{1};}};\n");
%!     fprintf(fd, "Transfinite Surface(2) = {PointsOf{Surface{2};}};\n");
%!     fprintf(fd, "Transfinite Surface(3) = {PointsOf{Surface{3};}};\n");
%!     fprintf(fd, "Transfinite Surface(4) = {PointsOf{Surface{4};}};\n");
%!     fprintf(fd, "Transfinite Surface(5) = {PointsOf{Surface{5};}};\n");
%!     fprintf(fd, "Transfinite Surface(6) = {PointsOf{Surface{6};}};\n");
%!     fprintf(fd, "Transfinite Surface(7) = {PointsOf{Surface{7};}};\n");
%!     fprintf(fd, "Transfinite Surface(8) = {PointsOf{Surface{8};}};\n");
%!     fprintf(fd, "Transfinite Surface(9) = {PointsOf{Surface{9};}};\n");
%!     fprintf(fd, "Transfinite Surface(10) = {PointsOf{Surface{10};}};\n");
%!     fprintf(fd, "Transfinite Surface(11) = {PointsOf{Surface{11};}};\n");
%!     fprintf(fd, "Transfinite Surface(12) = {PointsOf{Surface{12};}};\n");
%!     endif
%!     switch (param.elem_type)
%!     case {"iso8", "iso20", "iso20r", "iso27", "penta15", "penta18"}
%!       extrude_opt = "Layers{1}; Recombine;";
%!     otherwise
%!       extrude_opt = "";
%!     endswitch
%!     fprintf(fd, "v1[] = Extrude {0, h, 0}{ Surface{1}; %s };\n", extrude_opt);
%!     fprintf(fd, "v2[] = Extrude {0, h, 0}{ Surface{2}; %s };\n", extrude_opt);
%!     fprintf(fd, "v3[] = Extrude {0, h, 0}{ Surface{3}; %s };\n", extrude_opt);
%!     fprintf(fd, "v4[] = Extrude {0, h, 0}{ Surface{4}; %s };\n", extrude_opt);
%!     fprintf(fd, "v5[] = Extrude {0, h, 0}{ Surface{5}; %s };\n", extrude_opt);
%!     fprintf(fd, "v6[] = Extrude {0, h, 0}{ Surface{6}; %s };\n", extrude_opt);
%!     fprintf(fd, "v7[] = Extrude {0, h, 0}{ Surface{7}; %s };\n", extrude_opt);
%!     fprintf(fd, "v8[] = Extrude {0, h, 0}{ Surface{8}; %s };\n", extrude_opt);
%!     fprintf(fd, "v9[] = Extrude {0, h, 0}{ Surface{9}; %s };\n", extrude_opt);
%!     fprintf(fd, "v10[] = Extrude {0, h, 0}{ Surface{10}; %s };\n", extrude_opt);
%!     fprintf(fd, "v11[] = Extrude {0, h, 0}{ Surface{11}; %s };\n", extrude_opt);
%!     fprintf(fd, "v12[] = Extrude {0, h, 0}{ Surface{12}; %s };\n", extrude_opt);
%!     switch (param.elem_type)
%!     case {"iso8", "iso20", "iso20r", "iso27"}
%!       fprintf(fd, "Recombine Surface{1, v1[0]};\n");
%!       fprintf(fd, "Recombine Surface{2, v2[0]};\n");
%!       fprintf(fd, "Recombine Surface{3, v3[0]};\n");
%!       fprintf(fd, "Recombine Surface{4, v4[0]};\n");
%!       fprintf(fd, "Recombine Surface{5, v5[0]};\n");
%!       fprintf(fd, "Recombine Surface{6, v6[0]};\n");
%!       fprintf(fd, "Recombine Surface{7, v7[0]};\n");
%!       fprintf(fd, "Recombine Surface{8, v8[0]};\n");
%!       fprintf(fd, "Recombine Surface{9, v9[0]};\n");
%!       fprintf(fd, "Recombine Surface{10, v10[0]};\n");
%!       fprintf(fd, "Recombine Surface{11, v11[0]};\n");
%!       fprintf(fd, "Recombine Surface{12, v12[0]};\n");
%!     endswitch
%!     fprintf(fd, "Physical Volume(\"volume\", 1) = {1, 2, 3, 4, 7, 8, 9, 5, 10, 6, 11, 12};\n");
%!     fprintf(fd, "Physical Surface(\"clamp\", 2) = {54, 98, 186, 274};\n");
%!     fprintf(fd, "Physical Surface(\"pressure\", 3) = {68, 112, 134, 152};\n");
%!     fprintf(fd, "Physical Surface(\"displacement\", 4) = {156, 244};\n");
%!     fprintf(fd, "Physical Surface(\"stress\", 5) = {112, 134};\n");
%!     fprintf(fd, "Physical Surface(\"top\", 6) = {1,2,3,4,5,6,7,8,9,10,11,12};\n");
%!     fprintf(fd, "Physical Surface(\"bottom\", 7) = {v1[0],v2[0],v3[0],v4[0],v5[0],v6[0],v7[0],v8[0],v9[0],v10[0],v11[0],v12[0]};\n");
%!     fprintf(fd, "Physical Curve(\"clampedge\", 8) = {40, 269};\n");
%!     if (~param.transfinite)
%!     fputs(fd, "MeshSize{PointsOf{Volume{1,2,3,4,5,6,7,8,9,10,11,12};}} = h;\n");
%!     endif
%!     switch (param.elem_type)
%!     case {"penta15", "iso20", "iso20r"}
%!       fprintf(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     otherwise
%!       fprintf(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     endswitch
%!     if (~param.transfinite)
%!     switch (param.elem_type)
%!     case {"penta15", "tet10", "tet10h", "tet20"}
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     endswitch
%!     endif
%!     switch (param.elem_type)
%!     case {"iso8"}
%!       fprintf(fd, "Mesh.ElementOrder = 1;\n");
%!     case {"tet20"}
%!       fprintf(fd, "Mesh.ElementOrder = 3;\n");
%!     otherwise
%!       fprintf(fd, "Mesh.ElementOrder = 2;\n");
%!     endswitch
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   #spawn_wait(spawn("gmsh", {[filename, ".geo"]}));
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   switch (param.elem_type)
%!   case {"iso8"}
%!     param.elem_type_surf = {"iso4"};
%!   case "iso20"
%!     param.elem_type_surf = {"quad8"};
%!   case "iso20r"
%!     param.elem_type_surf = {"quad8r"};
%!   case {"iso27"}
%!     param.elem_type_surf = {"quad9"};
%!   case {"penta15"}
%!     param.elem_type_surf = {"quad8", "tria6h"};
%!   case "tet10"
%!     param.elem_type_surf = {"tria6"};
%!   case "tet10h"
%!     param.elem_type_surf = {"tria6h"};
%!   case "tet20"
%!     param.elem_type_surf = {"tria10"};
%!   case "penta18"
%!     param.elem_type_surf = {"quad9", "tria6h"};
%!   endswitch
%!   switch (param.elem_type)
%!   case {"iso8"}
%!     param.elem_type_line = "line2";
%!   case {"tet20"}
%!     param.elem_type_line = "line4";
%!   otherwise
%!     param.elem_type_line = "line3";
%!   endswitch
%!   opt_msh.elem_type = {param.elem_type, param.elem_type_surf{:}, param.elem_type_line, "point1"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.pressure = struct();
%!   mesh.materials = struct();
%!   for i=1:numel(param.elem_type_surf)
%!   switch (param.bound_cond)
%!   case "surface"
%!     grp_id_clamp = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 2);
%!     if (~isempty(grp_id_clamp))
%!       load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_surf{i})(grp_id_clamp).nodes, 1:3) = true;
%!     endif
%!   case "edge"
%!     grp_id_clamp = find([[getfield(mesh.groups, param.elem_type_line)].id] == 8);
%!     if (~isempty(grp_id_clamp))
%!       load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_line)(grp_id_clamp).nodes, 1:3) = true;
%!     endif
%!   endswitch
%!   switch (param.analysis)
%!   case "plain strain"
%!   if (isfield(mesh.groups, param.elem_type_surf{i}))
%!   grp_id_top = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 6);
%!   grp_id_bottom = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 7);
%!   if (~isempty(grp_id_top))
%!     load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_surf{i})(grp_id_top).nodes, 2) = true;
%!   endif
%!   if (~isempty(grp_id_bottom))
%!     load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_surf{i})(grp_id_bottom).nodes, 2) = true;
%!   endif
%!   endif
%!   endswitch
%!   grp_id_p1 = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 3);
%!   if (~isempty(grp_id_p1))
%!     elem_id_p1 = getfield(mesh.groups, param.elem_type_surf{i})(grp_id_p1).elements;
%!     elno_p1 = getfield(mesh.elements, param.elem_type_surf{i})(elem_id_p1, :);
%!     press_load.elements = elno_p1;
%!     press_load.p = [repmat(param.p1, rows(elno_p1), columns(elno_p1))];
%!     load_case.pressure = setfield(load_case.pressure, param.elem_type_surf{i}, press_load);
%!   endif
%!   endfor
%!   mesh.materials = setfield(mesh.materials, param.elem_type, ones(rows(getfield(mesh.elements, param.elem_type)), 1, "int32"));
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).rho = param.rho;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
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
%!   delta = nan(1, numel(param.elem_type_surf));
%!   sigma1_max = nan(1, numel(param.elem_type_surf));
%!   for k=1:numel(param.elem_type_surf)
%!   grp_id_displacement = find([[getfield(mesh.groups, param.elem_type_surf{k})].id] == 4);
%!   if (~isempty(grp_id_displacement))
%!   if (isfield(mesh.groups, param.elem_type_surf{k}))
%!     elem_id_displacement = getfield(mesh.groups, param.elem_type_surf{k})(grp_id_displacement).elements;
%!     elno_id_displacement = getfield(mesh.elements, param.elem_type_surf{k})(elem_id_displacement, :);
%!     delta(k) = mean(sol_stat.def(elno_id_displacement, 3, end));
%!   endif
%!   endif
%!   grp_id_stress = find([[getfield(mesh.groups, param.elem_type_surf{k})].id] == 5);
%!   if(~isempty(grp_id_stress))
%!   elem_id_stress = getfield(mesh.groups, param.elem_type_surf{k})(grp_id_stress).elements;
%!   elno_id_stress = getfield(mesh.elements, param.elem_type_surf{k})(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(getfield(mesh.elements, param.elem_type) == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(getfield(sol_stat.stress.taum, param.elem_type)(ridx(j), cidx(j), :, end), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max(k) = max(sigma1_max(k), max(eig(TAU)));
%!   endfor
%!   endif
%!   endfor
%!   delta = mean(delta(isfinite(delta)));
%!   sigma1_max = mean(sigma1_max(isfinite(sigma1_max)));
%!   fprintf(stderr, "mesh size=%.1f\n", param.h);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   ## K.J.Bathe page 329 chaper 4.4
%!   sigma1_max_ref = 0.6056e6 / SI_unit_pascal;
%!   delta_ref = -1.669e-3 / SI_unit_meter;
%!   tol_sigma = 2.5e-2;
%!   tol_delta = 1.25e-2;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   assert_simple(sigma1_max, sigma1_max_ref, tol_sigma * abs(sigma1_max_ref));
%!   assert_simple(delta, delta_ref, tol_delta * abs(delta_ref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
