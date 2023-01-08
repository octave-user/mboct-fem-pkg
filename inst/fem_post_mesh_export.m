## Copyright (C) 2019(-2023) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} fem_post_mesh_export(@var{filename}, @var{mesh})
## @deftypefnx {} fem_post_mesh_export(@var{filename}, @var{mesh}, @var{options})
## @deftypefnx {} fem_post_mesh_export(@var{filename}, @var{mesh}, @var{options}, @var{load_case})
## @deftypefnx {} fem_post_mesh_export(@var{filename}, @var{mesh}, @var{options}, @var{load_case}, @var{dof_map})
## Export a finite element mesh to Gmsh- or APDL-format.
##
## @var{filename} @dots{} Finite element mesh output filename.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{options}.format @dots{} One of "gmsh", "apdl"
##
## @var{options}.elem_types @dots{} Subset of element types to be exported (used only for format "gmsh")
##
## @var{load_case} @dots{} Scalar struct for loads and constraints (used only for format "apdl")
##
## @var{dof_map} @dots{} Degree of freedom mapping
##
## @seealso{fem_pre_mesh_import}
## @end deftypefn

function fem_post_mesh_export(filename, mesh, options, load_case, dof_map)
  if (nargin < 2 || nargin > 5 || nargout > 0)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "format"))
    options.format = "gmsh";
  endif

  if (nargin < 4)
    load_case = struct();
  endif

  if (~isstruct(load_case) || numel(load_case) ~= 1)
    error("invalid argument for load_case");
  endif

  if (nargin < 5)
    switch (options.format)
      case "apdl"
        dof_map = fem_ass_dof_map(mesh, load_case);
      otherwise
        dof_map = struct();
    endswitch
  endif

  fd = -1;

  unwind_protect
    [fd] = fopen(filename, "w");

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    switch (options.format)
      case "gmsh"
        fem_export_gmsh(fd, filename, mesh, options, load_case, dof_map);
      case "apdl"
        fem_export_apdl(fd, filename, mesh, options, load_case, dof_map);
      otherwise
        error("unknown format \"%s\"", options.format);
    endswitch
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

function fem_export_apdl(fd, filename, mesh, options, load_case, dof_map)
  fprintf(fd, "/CLEAR\n");
  fprintf(fd, "/BATCH\n");
  fprintf(fd, "/TITLE,\"%s\"\n", filename);
  fprintf(fd, "/PREP7\n");

  fprintf(fd, "\n!! MATERIAL DATA\n\n");

  for i=1:numel(mesh.material_data)
    if (isfield(mesh.material_data(i), "E") && ~isempty(mesh.material_data(i).E))
      fprintf(fd, "MP,EX,%d,%.16e\n", i, mesh.material_data(i).E);
    endif

    if (isfield(mesh.material_data(i), "nu") && ~isempty(mesh.material_data(i).nu))
      fprintf(fd, "MP,NUXY,%d,%.16e\n", i, mesh.material_data(i).nu);
    endif

    if (isfield(mesh.material_data(i), "rho") && ~isempty(mesh.material_data(i).rho))
      fprintf(fd, "MP,DENS,%d,%.16e\n", i, mesh.material_data(i).rho);
    endif
  endfor

  fprintf(fd, "\n!! NODES\n\n");
  fprintf(fd, "N,%d,%.16e,%.16e,%.16e\n", [1:rows(mesh.nodes); mesh.nodes(:, 1:3).']);

  ielem_type = int32(0);
  ielem_num = int32(0);

  if (isfield(mesh.elements, "tet10"))
    fprintf(fd, "\n!! TET10\n\n");

    fprintf(fd, "ET,%d,SOLID187\n", ++ielem_type);
    fprintf(fd, "TYPE,%d\n", ielem_type);
    fprintf(fd, "KEYO,%d,6,0\n", ielem_type);
    fprintf(fd, "MAT,%d\nEN,%d,%d,%d,%d,%d,%d,%d,%d,%d\nEMORE,%d,%d\n", [mesh.materials.tet10.'; ielem_num + (1:rows(mesh.elements.tet10)); mesh.elements.tet10.']);

    ielem_num += rows(mesh.elements.tet10);
  endif

  if (isfield(mesh.elements, "rbe3"))
    fprintf(fd, "\n!! RBE3\n\n");

    for i=1:numel(mesh.elements.rbe3)
      cid = ++ielem_type;
      tid = ++ielem_type;

      fprintf(fd, "et,%d,174\n", cid);
      fprintf(fd, "et,%d,170\n", tid);
      fprintf(fd, "keyo,%d,2,1\n", tid);
      fprintf(fd, "keyo,%d,4,111111\n", tid);
      fprintf(fd, "keyo,%d,5,5\n", tid);
      fprintf(fd, "keyo,%d,12,5\n", cid);
      fprintf(fd, "keyo,%d,4,1\n", cid);
      fprintf(fd, "keyo,%d,2,2\n", cid); ## MPC constraint
      fprintf(fd, "eblock,%d,,,%d\n", 10, numel(mesh.elements.rbe3(i).elements.tria6));
      fprintf(fd, "(15i9)\n");
      fprintf(fd, "%9d%9d%9d%9d%9d%9d%9d%9d%9d%9d%9d%9d%9d\n", ...
              [ielem_num + (1:numel(mesh.elements.rbe3(i).elements.tria6));
               repmat(cid, 3, numel(mesh.elements.rbe3(i).elements.tria6));
               zeros(1, numel(mesh.elements.rbe3(i).elements.tria6));
               mesh.elements.tria6(mesh.elements.rbe3(i).elements.tria6, [1, 2, 3, 3, 4, 5, 3, 6]).']);
      fprintf(fd, "%d\n", -1);
      fprintf(fd, "tshape\n");

      ielem_num += numel(mesh.elements.rbe3(i).elements.tria6);

      fprintf(fd, "type,%d\n", tid);
      fprintf(fd, "mat,%d\n", cid);
      fprintf(fd, "real,%d\n", cid);
      fprintf(fd, "tshape,pilo\n");
      fprintf(fd, "en,%d,%d\n", ++ielem_num, mesh.elements.rbe3(i).nodes(1));
      fprintf(fd, "tshape\n\n");
    endfor
  endif

  if (isfield(mesh.elements, "joints"))
    fprintf(fd, "\n!! JOINTS\n\n");

    for i=1:numel(mesh.elements.joints)
      joint_type = "simple";

      for j=1:rows(mesh.elements.joints(i).C)
        [dof, fact] = find(mesh.elements.joints(i).C(j, :));

        if (numel(fact) > 1)
          joint_type = "generic";
          break;
        endif
      endfor

      switch (joint_type)
        case "simple"
          for j=1:rows(mesh.elements.joints(i).C)
            if (isfield(load_case, "joints"))
              D = load_case.joints(i).U(j);
            else
              D = 0;
            endif

            [equ, dof, fact] = find(mesh.elements.joints(i).C(j, :));
            inode = floor((dof - 1) / 6) + 1;
            idof = mod((dof - 1), 6) + 1;

            fprintf(fd, "D,%d,%s,%.16e\n", mesh.elements.joints(i).nodes(inode), {"UX","UY","UZ","ROTX","ROTY","ROTZ"}{idof}, D / fact);
          endfor
        case "generic"
          for j=1:rows(mesh.elements.joints(i).C)
            if (isfield(load_case, "joints"))
              D = load_case.joints(i).U(j);
            else
              D = 0;
            endif

            fprintf(fd, "CE,NEXT,%.16e", D);

            inode = int32(0);

            for k=1:numel(mesh.elements.joints(i).nodes)
              for l=1:6
                C = mesh.elements.joints(i).C(j, (k - 1) * 6 + l);
                if (C ~= 0 && dof_map.ndof(mesh.elements.joints(i).nodes(k), l) > 0)
                  if (++inode > 3)
                    fprintf(fd, "\nCE,HIGH,%.16e", D);
                    inode = int32(0);
                  endif
                  fprintf(fd, ",%d,%s,%.16e", mesh.elements.joints(i).nodes(k), {"UX","UY","UZ","ROTX","ROTY","ROTZ"}{l}, C);
                endif
              endfor
            endfor
            fprintf(fd, "\n");
          endfor
      endswitch
    endfor
  endif

  if (isfield(load_case, "pressure"))
    if (isfield(load_case.pressure, "tria6"))
      fprintf(fd, "\n!! PRESSURE LOADS\n\n");

      fprintf(fd, "ET,%d,SURF154\n", ++ielem_type);
      fprintf(fd, "TYPE,%d\n", ielem_type);

      fprintf(fd, "EN,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", [ielem_num + (1:rows(load_case.pressure.tria6.elements));
                                                      load_case.pressure.tria6.elements(:, [1, 2, 3, 3, 4, 5, 3, 6]).']);

      fprintf(fd, "\n");

      fprintf(fd, "SFE,%d,1,PRES,0,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n", [double(ielem_num) + (1:rows(load_case.pressure.tria6.p));
                                                                load_case.pressure.tria6.p(:, [1, 2, 3, 3, 4, 5, 3, 6]).']);

      ielem_num += rows(load_case.pressure.tria6.elements);
    endif
  endif

  if (isfield(load_case, "loads"))
    fprintf(fd, "\n!! NODAL LOADS\n\n");

    for i=1:rows(load_case.loads)
      for j=1:columns(load_case.loads)
        fprintf(fd, "F,%d,%s,%.16e\n", load_case.loaded_nodes(i), {"FX","FY","FZ","MX","MY","MZ"}{j}, load_case.loads(i, j));
      endfor
    endfor

    fprintf(fd, "\n");
  endif

  if (isfield(load_case, "locked_dof"))
    [inode, idof, flag] = find(load_case.locked_dof);

    if (numel(inode))
      fprintf(fd, "\n!! NODAL CONSTRAINTS\n\n");

      for i=1:numel(inode)
        if (dof_map.ndof(inode(i), idof(i)) ~= -1)
          fprintf(fd, "D,%d,%s,0\n", inode(i), {"UX", "UY", "UZ", "ROTX", "ROTY", "ROTZ"}{idof(i)});
        endif
      endfor
    endif

    fprintf(fd, "\n");
  endif

  fprintf(fd, "ALLSEL,ALL\n");
  fprintf(fd, "/SOLU\n");
  fprintf(fd, "ANTYPE,STATIC\n");
  fprintf(fd, "NLGEO,OFF\n");
  fprintf(fd, "SOLVE\n");
  fprintf(fd, "SAVE\n");
  fprintf(fd, "/POST1\n");
  fprintf(fd, "PRNSOL,U,COMP\n");
  fprintf(fd, "PRNSOL,ROT,COMP\n");
  fprintf(fd, "PRNLD\n");
  fprintf(fd, "/EXIT\n");
endfunction

function fem_export_gmsh(fd, filename, mesh, options, load_case, dof_map)
  persistent gmsh_elem_types = {"iso8", "iso20", "iso20r", "iso4", "quad8", "tet4", "tet10", ...
                                "tria6", "tria3", "beam2", "beam3", "penta15", "tet10h", "tet20", "tria6h", "tria10"};

  if (~isfield(options, "elem_types"))
    options.elem_types = gmsh_elem_types;
  endif

  valid_elem_type = false(1, numel(options.elem_types));

  for i=1:numel(options.elem_types)
    switch (options.elem_types{i})
      case gmsh_elem_types
        valid_elem_type(i) = true;
    endswitch
  endfor

  options.elem_types = {options.elem_types{find(valid_elem_type)}};

  fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

  inumelem = int32(0);
  elem_types = fieldnames(mesh.elements);

  for i=1:numel(elem_types)
    switch (elem_types{i})
      case options.elem_types
        elem_i = getfield(mesh.elements, elem_types{i});
        if (isstruct(elem_i))
          inumelem += numel(elem_i);
        elseif (ismatrix(elem_i))
          inumelem += rows(elem_i);
        endif
    endswitch
  endfor

  inumgroups = int32(0);

  if (isfield(mesh, "groups"))
    group_types = fieldnames(mesh.groups);

    for i=1:numel(group_types)
      switch (group_types{i})
        case options.elem_types
          inumgroups += numel(getfield(mesh.groups, group_types{i}));
      endswitch
    endfor
  else
    group_types = {};
  endif

  fputs(fd, "$PhysicalNames\n");
  fprintf(fd, "%d\n", inumgroups);

  for i=1:numel(group_types)
    groups = getfield(mesh.groups, group_types{i});

    switch (group_types{i})
      case {"iso8", "iso20", "iso20r", "tet10", "penta15", "tet10h", "tet20"}
        group_dim = 3;
      case {"iso4", "quad8", "tria6", "tria3", "tria6h", "tria10"}
        group_dim = 2;
      case {"beam2", "beam3"}
        group_dim = 1;
      otherwise
        continue;
    endswitch

    for j=1:numel(groups)
      if (any(isspace(groups(j).name)))
        error("spaces in group names (%s) are not allowed for gmsh format", groups(j).name);
      endif
      fprintf(fd, "%d %d \"%s\"\n", group_dim, groups(j).id, groups(j).name);
    endfor
  endfor

  fputs(fd, "$EndPhysicalNames\n");

  fprintf(fd, "$Nodes\n%d\n", rows(mesh.nodes));
  fprintf(fd, "%d %.16e %.16e %.16e\n", [1:rows(mesh.nodes); mesh.nodes(:, 1:3).']);
  fputs(fd, "$EndNodes\n");

  fprintf(fd, "$Elements\n%d\n", inumelem);

  inumelem = int32(0);

  for i=1:numel(elem_types)
    switch (elem_types{i})
      case "iso8"
        elem_type_id = 5;
        elem_node_order([5:8, 1:4]) = 1:8;
      case "iso20"
        elem_type_id = 17;
        elem_node_order([5:8, 1:4, 17, 19, 20, 18, 9, 12, 14, 10, 11, 13, 15, 16]) = 1:20;
      case "iso20r"
        elem_type_id = 17;
        elem_node_order([1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16]) = 1:20;
      case "penta15"
        elem_type_id = 18;
        elem_node_order([1, 2, 3, 4, 5, 6, 7, 10, 8, 13, 15, 14, 9, 11, 12]) = 1:15;
      case "iso4"
        elem_type_id = 3;
        elem_node_order = 1:4;
      case "quad8"
        elem_type_id = 16;
        elem_node_order = 1:8;
      case {"tria6", "tria6h"}
        elem_type_id = 9;
        elem_node_order = 1:6;
      case "tria3"
        elem_type_id = 2;
        elem_node_order = 1:3;
      case "tria10"
        elem_type_id = 21;
        elem_node_order = 1:10;
      case {"tet10", "tet10h"}
        elem_type_id = 11;
        elem_node_order([1:8, 10, 9]) = 1:10;
      case "tet20"
        elem_type_id = 29;
        elem_node_order([1,5,6,2,7,8,3,9,10,17,16,20,14,19,12,18,15,13,11,4]) = 1:20;
      case "tet4"
        elem_type_id = 4;
        elem_node_order = 1:4;
      case "beam2"
        elem_type_id = 1;
        elem_node_order = 1:2;
      case "beam3"
        elem_type_id = 8;
        elem_node_order = 1:3;
      otherwise
        continue
    endswitch

    switch (elem_types{i})
      case {"beam2", "beam3"}
        beam_elem = getfield(mesh.elements, elem_types{i});
        if (isstruct(beam_elem))
          elem_nodes = reshape([beam_elem.nodes], numel(elem_node_order), numel(beam_elem)).';
        else
          elem_nodes = beam_elem;
        endif
      otherwise
        elem_nodes = getfield(mesh.elements, elem_types{i});
    endswitch

    elem_groups = zeros(rows(elem_nodes), 1, "int32");
    elem_tags = [repmat([elem_type_id; 2], 1, rows(elem_nodes)); zeros(2, rows(elem_nodes))];

    if (isfield(mesh, "groups") && isfield(mesh.groups, elem_types{i}))
      groups = getfield(mesh.groups, elem_types{i});
      for j=1:numel(groups)
        elem_tags(3:4, groups(j).elements) = groups(j).id;
      endfor
    endif

    numcols = rows(elem_tags) + columns(elem_nodes) + 1;
    format = "%d";

    for j=1:numcols - 1
      format = [format, " %d"];
    endfor

    format = [format, "\n"];
    fprintf(fd, format, [inumelem + (1:rows(elem_nodes));
                         elem_tags;
                         elem_nodes(:, elem_node_order).']);
    inumelem += rows(elem_nodes);
    clear elem_node_order elem_type_id;
  endfor

  fputs(fd, "$EndElements\n");
endfunction

%!demo
%! ## DEMO 1
%! ## K.J.Bathe 2002, page 328 4.20a
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!  filename(filename == "\\") = "/";
%! endif
%! fd = -1;
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "w");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! mesh_size = 20;
%! p1 = 0.006;
%! E = 55;
%! nu = 0.3;
%! rho = 1000e-12;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%! fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%! fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%! fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%! fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%! fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%! fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%! fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%! fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%! fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%! fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%! fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Circle(3) = {3,11,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Circle(7) = {7,12,8};\n");
%! fputs(fd, "Line(8) = {8,9};\n");
%! fputs(fd, "Line(9) = {9,10};\n");
%! fputs(fd, "Line(10) = {10,1};\n");
%! fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%! fputs(fd, "Plane Surface(14) = {11};\n");
%! fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; };\n", mesh_size);
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%! fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%! fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%! fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");
%! [~] = unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! [~] = unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! [~] = unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);

%! load_case.pressure.tria6.elements = elno_p1;
%! load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.rho = rho;
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! opts.format = "apdl";
%! fem_post_mesh_export([filename, ".dat"], mesh, opts, load_case, dof_map);
%! spawn_wait(spawn("cat", {[filename, ".dat"]}));
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
