## Copyright (C) 2019(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{load_case}] = fem_pre_mesh_import(@var{filename})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_import(@dots{}, @var{format})
## Load a mesh and a load case from a file. Supported values for @var{format} are:
##
## "gmsh" from C. Geuzaine and J.-F. Remacle (*.msh)
##
## "apdl" from ANSYS inc (*.dat, *.inp)
##
## "eossp" from Erich Payer (*.xfe)
##
## If no argument @var{format} is provided, it will be determined from the file extension.
##
## @end deftypefn

function [mesh, load_case] = fem_pre_mesh_import(filename, format)
  if (nargin < 1 || nargin > 2 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    [fdir, fname, fext] = fileparts(filename);

    switch (fext)
      case ".msh"
	format = "gmsh";
      case {".dat", ".inp"}
	format = "apdl";
      case ".xfe"
	format = "eossp";
      otherwise
	format = "*unknown*";
    endswitch
  endif

  switch format
    case "gmsh"
      mesh = fem_load_mesh_gmsh(filename, format);
    case "apdl"
      mesh = fem_load_mesh_apdl(filename, format);
    case "eossp"
      [mesh, load_case] = fem_load_mesh_eossp(filename, format);
    otherwise
      error("mesh format \"%s\" not supported", format);
  endswitch
endfunction

function mesh = fem_load_mesh_gmsh(filename, format)
  mesh.nodes = zeros(0, 6);
  mesh.elements = struct();
  mesh.groups = struct();

  [fd, msg] = fopen(filename, "rt");

  if (fd == -1)
    error("failed to open file \"%s\": %s", filename, msg);
  endif

  nodes = zeros(0, 4);
  elements = zeros(0, 15, "int32");
  p_name = {};
  p_dim = [];
  p_id = [];

  unwind_protect
    while (true)
      line = fgetl(fd);

      if (~ischar(line))
	break
      endif

      switch (line)
	case "$PhysicalNames"
	  [p_names_count, count, msg] = fscanf(fd, "%d\n", "C");

	  p_dim = p_id = zeros(p_names_count, 1, "int32");

          for i=1:p_names_count
            [p_dim(i), p_id(i), p_name{i}, count, msg] = fscanf(fd, "%d %d %s\n", "C");
            p_name{i} = p_name{i}(2:end - 1);
          endfor
        case "$Nodes"
          [num_nodes, count, msg] = fscanf(fd, "%d\n", "C");

          if (count ~= 1)
            error("failed to parse number of nodes in file \"%s\"", filename);
          endif

          [nodes, count, msg] = fscanf(fd, "%g", [4, num_nodes]);

          if (count ~= num_nodes * 4)
            error("failed to parse nodes in file \"%s\"", filename);
          endif

          mesh.nodes = [nodes(2:4, nodes(1, :)).', zeros(columns(nodes), 3)];
        case "$Elements"
          [num_elements, count, msg] = fscanf(fd, "%d\n", "C");

          if (count ~= 1)
            error("failed to parse number of elements in file \"%s\"", filename);
          endif

          elem_type = zeros(num_elements, 1, "int32");
          elem_tags = zeros(num_elements, 2, "int32");
          elements = zeros(num_elements, 10, "int32");

          for i=1:num_elements
            line = fgetl(fd);

            if (~ischar(line))
              error("unexpected end of file in \"%s\"", filename);
            endif

            [val, count, msg] = sscanf(line, "%d");

            elem_no = val(1);
            elem_type(elem_no) = val(2);
            elem_tags(elem_no, :) = val(4:val(3) + 3);
            node_idx = val(3) + 4:rows(val);
            elements(elem_no, 1:length(node_idx)) = val(node_idx);
          endfor

          groups = unique(elem_tags(:, 1));

          persistent eltype = struct("name", {"tet10", "tria6", "tet4", "prism6", "tria3", "iso8", "iso4"}, ...
                                     "promote", {-1, -1, 6, 6, 7, -1, -1}, ...
                                     "id", {11, 9, 4, 6, 2, 5, 3}, ...
                                     "nnodes", {10, 6, 8, 8, 4, 8, 4}, ...
                                     "dim", {3, 2, 3, 3, 2, 3, 2}, ...
                                     "norder", {[1:8, 10, 9], ...
                                                1:6, ...
                                                [4, 4, 4, 4, 1:3, 3], ...
                                                [6,4,4,5,3,1,1,2], ...
                                                [1:3, 3], ...
                                                [5:8, 1:4], ...
                                                1:4});

          for k=1:numel(eltype)
            idx_eltype = find(elem_type == eltype(k).id);

            if (~numel(idx_eltype))
              continue;
            endif

            idx_eltype_elem(idx_eltype) = 1:numel(idx_eltype);
            elnodes = elements(idx_eltype, eltype(k).norder);

            if (any(any(elnodes > rows(mesh.nodes) | elnodes < 1)))
              error("invalid node index in file \"%s\" for element type %s", filename, eltype(k).name);
            endif

            if (isfield(mesh.elements, eltype(k).name))
              elnodes = [getfield(mesh.elements, eltype(k).name);
                         elnodes];
            endif

            mesh.elements = setfield(mesh.elements, ...
                                     eltype(k).name, ...
                                     elnodes);

            for i=1:numel(groups)
              idx_elgrp = find(elem_type == eltype(k).id & elem_tags(:, 1) == groups(i));

              if (numel(idx_elgrp))
                if (isfield(mesh.groups, eltype(k).name))
                  mshgrp = getfield(mesh.groups, eltype(k).name);
                else
                  mshgrp = struct()([]);
                endif

                mshgrp(end + 1).id = groups(i);
                idx_grp_name = find((p_id == groups(i)) & (p_dim == eltype(k).dim));

                if (numel(idx_grp_name))
                  mshgrp(end).name = p_name{idx_grp_name};
                else
                  mshgrp(end).name = sprintf("unnamed-group[%d]", mshgrp(end).id);
                endif

                mshgrp(end).elements = idx_eltype_elem(idx_elgrp);
                mshgrp(end).nodes = unique(reshape(elnodes(mshgrp(end).elements, :), ...
                                                   1, ...
                                                   eltype(k).nnodes * numel(mshgrp(end).elements)));

                mesh.groups = setfield(mesh.groups, eltype(k).name, mshgrp);
              endif
            endfor
          endfor

          for k=1:numel(eltype)
            idp = eltype(k).promote;

            if (idp < 0)
              continue;
            endif

            if (~isfield(mesh.elements, eltype(k).name))
              continue;
            endif

            elnodes = getfield(mesh.elements, eltype(k).name);
            eloffset = int32(0);

            if (isfield(mesh.elements, eltype(idp).name))
              elnodesp = getfield(mesh.elements, eltype(idp).name);
              eloffset = rows(elnodesp);
              elnodes = [elnodesp; elnodes];
            endif

            if (isfield(mesh.groups, eltype(k).name))
              mshgrp = getfield(mesh.groups, eltype(k).name);
              for i=1:numel(mshgrp)
                mshgrp(i).elements += eloffset;
              endfor

              if (isfield(mesh.groups, eltype(idp).name))
                mshgrp = [getfield(mesh.groups, eltype(idp).name), ...
                          mshgrp];
              endif

              mesh.groups = setfield(rmfield(mesh.groups, eltype(k).name), eltype(idp).name, mshgrp);
            endif

            mesh.elements = setfield(rmfield(mesh.elements, eltype(k).name), eltype(idp).name, elnodes);
          endfor
      endswitch
    endwhile
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect

  if (~isfield(mesh, "nodes"))
    error("no nodes in file \"%s\"", filename);
  endif

  if (~isfield(mesh, "elements"))
    error("no elements in file \"%s\"", filename);
  endif
endfunction

function mesh = fem_load_mesh_apdl(filename, format)
  mesh = struct();

  [fd, msg] = fopen(filename, "rt");

  if (fd == -1)
    error("failed to open file \"%s\": %s", filename, msg);
  endif

  nodes = zeros(0, 4);
  elements = zeros(0, 15, "int32");
  elemtype = -1;

  unwind_protect
    while (true)
      line = fgetl(fd);

      if (~ischar(line))
	break
      endif

      idx = strfind(line, "nblock,");

      if (numel(idx) && idx == 1)
	[numcols, count, msg] = sscanf(line, "nblock,%d", "C");

	if (count ~= 1)
	  error("parse error in file \"%s\": %s (%s)", filename, line, msg);
	endif

	line = fgetl(fd);

	[nodes, count, msg] = fscanf(fd, "%g", [numcols + 1, inf]);

	if (count == 0)
	  error("parse error in file \"%s\": %s", filename, msg);
	endif

	idx_node = find(nodes(1, :) > 0);

	if (~numel(idx_node))
	  error("parse error in file \"%s\": invalid nblock", filename);
	endif

	mesh.nodes = [nodes(2:4, nodes(1, idx_node)).', zeros(numel(idx_node), 3)];

	continue;
      endif

      idx = strfind(line, "et,");

      if (numel(idx) && idx == 1)
	[matid, elemtype, count, msg] = sscanf(line, "et,%d,%d", "C");

	if (count ~= 2)
	  error("parse error in file \"%s\": %s (%s)", filename, line, msg);
	endif

	continue;
      endif

      idx = strfind(line, "eblock");

      if (numel(idx) && idx == 1)
	[numcols, count, msg] = sscanf(line, "eblock,%d", "C");

	if (count ~= 1)
	  error("parse error in file \"%s\": %s (%s)", filename, line, msg);
	endif

	switch (elemtype)
	  case SOLID187
	    numcols += 2;
	  case SURF154
	    numcols += 3;
	endswitch

	line = fgetl(fd);

	[elements, count, msg] = fscanf(fd, "%d", [numcols, inf]);

	if (count == 0)
	  error("parse error in file \"%s\": %s", filename, msg);
	endif

	idx_elem = find(elements(1, :) > 0);

	if (~numel(idx_elem))
	  error("parse error in file \"%s\": invalid eblock", filename);
	endif

	switch (elemtype)
	  case SOLID187
	    if (~isfield(mesh, "elements") || ~isfield(mesh.elements, "tet10"))
	      mesh.elements.tet10 = zeros(0, 10, "int32");
	      mesh.materials.tet10 = zeros(0, 1, "int32");
	      mesh.groups.tet10 = struct("id", [], "name", [], "nodes", [], "elements", [])([]);
	    endif

	    mesh.groups.tet10(end + 1).id = matid;
	    mesh.groups.tet10(end).name = sprintf("solid%d[%d]", elemtype, matid);
	    mesh.groups.tet10(end).elements = int32(rows(mesh.elements.tet10) + (1:numel(idx_elem)));
	    mesh.groups.tet10(end).nodes = unique(reshape(elements(12:21, idx_elem), 1, 10 * numel(idx_elem)));

	    mesh.elements.tet10 = [mesh.elements.tet10;
				   int32(elements(12:21, idx_elem)).'];

	    mesh.materials.tet10 = [mesh.materials.tet10;
				    int32(elements(1, idx_elem)).'];
	  case SURF154
	    ## ANSYS Mechanical Workbench is creating degenerated quad8 elements instead of tria6 elements
	    idx_tria6 = find(elements(8, idx_elem) == elements(9, idx_elem) & ...
			     elements(12, idx_elem) == elements(8, idx_elem));

	    if (numel(idx_tria6))
	      if (~isfield(mesh, "elements") || ~isfield(mesh.elements, "tria6"))
		mesh.elements.tria6 = zeros(0, 6, "int32");
		mesh.groups.tria6 = struct("id", [], "name", [], "nodes", [], "elements", [])([]);
	      endif

	      mesh.groups.tria6(end + 1).id = matid;
	      mesh.groups.tria6(end).name = sprintf("surf%d[%d]", elemtype, matid);
	      mesh.groups.tria6(end).elements = int32(rows(mesh.elements.tria6) + (1:numel(idx_elem)));
	      mesh.groups.tria6(end).nodes = unique(reshape(elements([6,7,8,10,11,13], idx_elem(idx_tria6)), 1, 6 * numel(idx_elem)));

	      mesh.elements.tria6 = [mesh.elements.tria6;
				     int32(elements([6,7,8,10,11,13], idx_elem(idx_tria6))).'];
	    endif
	  otherwise
	    warning("unknown element type %d\n", elemtype);
	endswitch

	continue;
      endif

      idx = strfind(line, "CMBLOCK,");

      if (numel(idx) && idx == 1)
	cm_rem = line;
	cm_idx = 0;
	cm_name = "";
	cm_type = "";
	cm_size = 0;

	while numel(cm_rem)
	  [cm_tok, cm_rem] = strtok(cm_rem, ",");

	  switch ++cm_idx
	    case 2
	      cm_name = cm_tok;
	    case 3
	      cm_type = cm_tok;
	    case 4
	      [cm_size, count, msg] = sscanf(cm_tok, "%d", "C");

	      if (count ~= 1)
		error("parse error in file \"%s\": %s (%s)", filename, line, msg);
	      endif
	  endswitch
	endwhile

	line = fgetl(fd);

	if (cm_size > 0)
	  [nodes, count, msg] = fscanf(fd, "%d", cm_size);

	  switch cm_type
	    case "NODE"
	      if (~isfield(mesh, "groups") || ~isfield(mesh.groups, "nodes"))
		mesh.groups.nodes = struct("id",[],"name",[],"nodes",[])([]);
	      endif

	      mesh.groups.nodes(end + 1).id = numel(mesh.groups.nodes) + 1;
	      mesh.groups.nodes(end).name = cm_name;
	      mesh.groups.nodes(end).nodes = int32(nodes);
	  endswitch
	endif

	continue;
      endif
    endwhile
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect

  if (~isfield(mesh, "nodes"))
    error("no nodes in file \"%s\"", filename);
  endif

  if (~isfield(mesh, "elements"))
    error("no elements in file \"%s\"", filename);
  endif
endfunction

function id = SOLID187()
  id = 187;
endfunction

function id = SURF154()
  id = 154;
endfunction

function [mesh, load_case] = fem_load_mesh_eossp(filename, format)
  persistent delim = " \t0123456789";

  [fd, msg] = fopen(filename, "rt");

  if (fd == -1)
    error("failed to open file \"%s\"", filename);
  endif

  nodes = zeros(4, 0);
  node_id = zeros(1, 0, "int32");
  material = zeros(4, 0);
  material_id = zeros(1, 0, "int32");

  mesh.nodes = zeros(0, 6);
  mesh.material_data = struct("E", [], "nu", [], "rho", [], "C", [])([]);

  load_case.locked_dof = false(0, 6);
  load_case.loads = zeros(0, 6);
  load_case.loaded_nodes = zeros(0, 1, "int32");

  unwind_protect
    while (true)
      line = fgetl(fd);

      if (~ischar(line))
	break;
      endif

      [tok1, rem] = strtok(line, delim);
      [tok2, rem] = strtok(rem, delim);

      switch ([tok1, " ", tok2])
        case "KOOR KNR"
          nodes(:, end + 1) = sscanf(line, "KOOR KNR %d X %g Y %g  Z %g\n");
        case "MATE MNR"
          if (numel(mesh.nodes) == 0)
            mesh.node_id = int32(nodes(1, :)).';
            node_id(nodes(1, :)) = int32(1:columns(nodes));
            mesh.nodes = [nodes(2:4, :).', zeros(columns(nodes), 3)];
          endif

          material(:, end + 1) = sscanf(line, "MATE MNR %d E %g NUE %g RHO %g\n");
        case "STRU BNR"
          if (numel(mesh.material_data) == 0)
            material_id(material(1, :)) = int32(1:columns(material));

            for i=1:columns(material)
              mesh.material_data(i).E = material(2, i);
              mesh.material_data(i).nu = material(3, i);
              mesh.material_data(i).rho = material(4, i);
              mesh.material_data(i).C = fem_pre_mat_isotropic(mesh.material_data(i).E, mesh.material_data(i).nu);
            endfor
          endif

          while (ischar(line))
            [elem1, count] = sscanf(line, "STRU BNR %d ENR %d TNR %d MNR %d KNR1 %d KNR2 %d KNR3 %d KNR4 %d");

            if (count)
              line = fgetl(fd);

              if (~ischar(line))
                break;
              endif

              switch (elem1(3))
                case 11
                  [elem2, count] = sscanf(line, " KNR5 %d KNR6 %d");
                case 12
                  [elem2, count] = sscanf(line, " KNR5 %d KNR6 %d KNR7 %d KNR8 %d");
                otherwise
                  elem2 = [];
              endswitch
            else
              break;
            endif

            elem_nodes = [elem1(5:8).', elem2.'];

            switch (elem1(3))
              case 8
                mesh.materials.iso4(end + 1, 1) = material_id(elem1(4));
                mesh.elements.iso4(end + 1, :) = node_id(elem_nodes(1:4));
              case 11
                mesh.materials.iso8(end + 1, 1) = material_id(elem1(4));
                mesh.elements.iso8(end + 1, :) = node_id(elem_nodes([4:6, 4, 1:3, 1]));
              case 12
                mesh.materials.iso8(end + 1, 1) = material_id(elem1(4));
                mesh.elements.iso8(end + 1, :) = node_id(elem_nodes([5:8, 1:4]));
            endswitch

            line = fgetl(fd);
          endwhile
        case "RBED KNR"
          while (ischar(line))
            [constr_node_id, vx, vy, vz, dx, dy, dz, count] = sscanf(line, "RBED KNR %d VX %d VY %d VZ %d DX %d DY %d DZ %d", "C");

            if (count ~= 7)
              break;
            endif

            load_case(1).locked_dof(node_id(constr_node_id), 1:6) = logical([vx, vy, vz, dx, dy, dz]);

            line = fgetl(fd);
          endwhile
        case "LAST LFNR"
          while (ischar(line))
            [load_case_id, loaded_node_id, Fx, Fy, Fz, count] = sscanf(line, "LAST LFNR %d KNR %d FX %g FY %g FZ %g", "C");

            if (count ~= 5)
              break;
            endif

            if (numel(load_case) < load_case_id)
              load_case(end + 1:load_case_id) = repmat(load_case(1), 1, load_case_id - numel(load_case));
            endif

	    load_case(load_case_id).loads(end + 1, 1:3) = [Fx, Fy, Fz];
	    load_case(load_case_id).loaded_nodes(end + 1, 1) = node_id(loaded_node_id);

            line = fgetl(fd);
          endwhile
      endswitch
    endwhile
  unwind_protect_cleanup
    fclose(fd);
  end_unwind_protect
endfunction

%!test
%! ### TEST1
%! close all;
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 30e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 3.5e-3;
%! p = 25e6;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! if (do_rotate)
%! fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! endif
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! if (~do_rotate)
%! group_defs(1).id = 1;
%! group_defs(1).name = "box1";
%! group_defs(1).R = eye(3);
%! group_defs(1).X0 = zeros(3, 1);
%! group_defs(1).type = "box";
%! group_defs(1).geometry.xmin = 0;
%! group_defs(1).geometry.xmax = 0;
%! group_defs(1).geometry.ymin = 0;
%! group_defs(1).geometry.ymax = b;
%! group_defs(1).geometry.zmin = 0;
%! group_defs(1).geometry.zmax = c;
%! group_defs(1).elem_type = "tria6";
%! group_defs(2).id = 2;
%! group_defs(2).name = "cylinder1";
%! group_defs(2).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%! group_defs(2).type = "cylinder";
%! group_defs(2).geometry.rmin = 0;
%! group_defs(2).geometry.rmax = 0.5 * c;
%! group_defs(2).geometry.zmin = -0.5 * b;
%! group_defs(2).geometry.zmax = 0.5 * b;
%! group_defs(2).elem_type = "tria6";
%! group_defs(3).id = 3;
%! group_defs(3).name = "cylinder2";
%! group_defs(3).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(3).type = "cylinder";
%! group_defs(3).geometry.rmin = 0;
%! group_defs(3).geometry.rmax = 0.5 * c;
%! group_defs(3).geometry.zmin = -0.5 * b;
%! group_defs(3).geometry.zmax = 0.5 * b;
%! group_defs(3).elem_type = "tria6";
%! group_defs(4).id = 4;
%! group_defs(4).name = "box2";
%! group_defs(4).R = eye(3);
%! group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(4).type = "box";
%! group_defs(4).geometry.xmin = 0;
%! group_defs(4).geometry.xmax = 0;
%! group_defs(4).geometry.ymin = -0.5 * b;
%! group_defs(4).geometry.ymax = 0.5 * b;
%! group_defs(4).geometry.zmin = -0.5 * c;
%! group_defs(4).geometry.zmax = 0.5 * c;
%! group_defs(4).elem_type = "tria6";
%! groups = fem_pre_mesh_groups_create(mesh, group_defs, sqrt(eps));
%! assert(numel(groups.tria6), 4);
%! assert([groups.tria6.id], [group_defs.id]);
%! assert(groups.tria6(1).nodes, mesh.groups.tria6(1).nodes);
%! assert(groups.tria6(2).nodes, mesh.groups.tria6(1).nodes);
%! assert(groups.tria6(3).nodes, mesh.groups.tria6(2).nodes);
%! assert(groups.tria6(4).nodes, mesh.groups.tria6(2).nodes);
%! assert(groups.tria6(1).elements, mesh.groups.tria6(1).elements);
%! assert(groups.tria6(2).elements, mesh.groups.tria6(1).elements);
%! assert(groups.tria6(3).elements, mesh.groups.tria6(2).elements);
%! assert(groups.tria6(4).elements, mesh.groups.tria6(2).elements);
%! endif
%! fprintf(stderr, "assembling matrices ...\n");
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%! load_case.pressure.tria6.elements = mesh.elements.tria6(mesh.groups.tria6(find([mesh.groups.tria6.id] == 2)).elements, :);
%! Xp = mesh.nodes(load_case.pressure.tria6.elements, 1:3);
%! xp = reshape(Xp(:, 1), rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%! yp = reshape(Xp(:, 2), rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%! zp = reshape(Xp(:, 3), rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%! load_case.pressure.tria6.p = p / 2 * (yp / b + zp / c); #repmat(p, rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.Mlumped, ...
%!  mat_ass.R, ...
%!  mat_ass.Rlumped, ...
%!  mtot] = fem_ass_matrix(mesh, ...
%!                         dof_map, ...
%!                         [FEM_MAT_STIFFNESS, ...
%!                          FEM_MAT_MASS, ...
%!                          FEM_MAT_MASS_LUMPED, ...
%!                          FEM_VEC_LOAD_CONSISTENT, ...
%!                          FEM_VEC_LOAD_LUMPED, ...
%!                          FEM_SCA_TOT_MASS], load_case);
%! assert(mtot, a * b * c * mesh.material_data.rho, sqrt(eps) * a * b * c * mesh.material_data.rho);
%! fprintf(stderr, "eigenanalysis ...\n");
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%! [sol_eig_lumped] = fem_sol_modal(mesh, dof_map, setfield(mat_ass, "M", mat_ass.Mlumped), number_of_modes);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_stat);
%! sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%! X = mesh.nodes(unique(load_case.pressure.tria6.elements), 1:3).';
%! dof_idx = dof_map.ndof(unique(load_case.pressure.tria6.elements), 1:3);
%! F_con = full(mat_ass.R(dof_idx)).';
%! F_lumped = full(mat_ass.Rlumped(dof_idx)).';
%! M_con = cross(X, F_con);
%! M_lumped = cross(X, F_lumped);
%! Ftot_con = sum(F_con, 2);
%! Mtot_con = sum(M_con, 2);
%! Ftot_lumped = sum(F_lumped, 2);
%! Mtot_lumped = sum(M_lumped, 2);
%! Fx_an = -(b * c * p) / 2; ## grind(integrate(integrate(-1/2*p*(y/b+z/c),z,0,c),y,0,b));
%! Mz_an = (7 * b^2 * c * p) / 24; ## grind(integrate(integrate(1/2*p*(y/b+z/c)*y,z,0,c),y,0,b));
%! My_an = -(7 * b * c^2 * p) / 24; ## grind(integrate(integrate(-1/2*p*(y/b+z/c)*z,z,0,c),y,0,b));
%! F_an = [Fx_an; 0; 0];
%! M_an = [0; My_an; Mz_an];
%! assert(Ftot_con, F_an, eps^0.9 * norm(F_an));
%! assert(Ftot_lumped, F_an, eps^0.9 * norm(F_an));
%! assert(Mtot_con, M_an, eps^0.9 * norm(M_an));
%! assert(Mtot_lumped, M_an, 5e-3 * norm(M_an));
%! f = sol_eig.f(:);
%! f_lumped = sol_eig_lumped.f(:);
%! f_ref = [8768.74;
%!          14636.1;
%!          21145.7;
%!          39712.8;
%!          43555.5;
%!          47909;
%!          62270.4;
%!          84324.4;
%!          92665.1;
%!          94563];
%! for i=1:length(f)
%!  fprintf(stderr, "mode %d f=%.0f f_lumped=%.0f\n", i, f(i), f_lumped(i));
%! endfor
%! assert(all(f_lumped <= f));
%! assert(f, f_ref, tol * max(f_ref));

%! fprintf(stderr, "plotting ...\n");
%! figure("visible","off");
%! hold on;
%! fem_post_sol_plot(mesh);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('undeformed mesh');
%! opts_plot.elem_types = {"tria6", "tet10"};
%! opts_plot.elem_groups.tria6 = [mesh.groups.tria6.id];
%! opts_plot.elem_groups.tet10 = [mesh.groups.tet10.id];
%! for i=1:min(number_of_modes, length(sol_eig.f))
%!      figure("visible", "off");
%!      hold on;
%!      fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, :, i), "rows")), i, opts_plot);
%!      view(30,30);
%!      xlabel('x [m]');
%!      ylabel('y [m]');
%!      zlabel('z [m]');
%!      grid on;
%!      grid minor on;
%!      title(sprintf("%d. eigenmode: %gHz",i,sol_eig.f(i)));
%! endfor

%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_stat, scale_eig/max(norm(sol_stat.def, "rows")), 1, opts_plot);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title("deformed mesh consistent load vector");

%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_stat_lumped, scale_eig/max(norm(sol_stat_lumped.def, "rows")), 1, opts_plot);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title("deformed mesh lumped load vector");

%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ## TEST2
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! p1 = 25.79e6;
%! p2 = 7.83e6;
%! p3 = 1.3758e6;
%! scale_def = 5e-3;
%! mesh_size = 3e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! grp_id_p3 = find([[mesh.groups.tria6].id] == 4);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elem_id_p2 = mesh.groups.tria6(grp_id_p2).elements;
%! elem_id_p3 = mesh.groups.tria6(grp_id_p3).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%! elno_p2 = mesh.elements.tria6(elem_id_p2, :);
%! elno_p3 = mesh.elements.tria6(elem_id_p3, :);

%! load_case.pressure.tria6.elements = [elno_p1; elno_p2; elno_p3];
%! load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1));
%!                               repmat(p2, rows(elno_p2), columns(elno_p2));
%!                               repmat(p3, rows(elno_p3), columns(elno_p3))];

%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS, ...
%!                                     FEM_VEC_LOAD_CONSISTENT, ...
%!                                     FEM_VEC_LOAD_LUMPED], ...
%!                                    load_case);
%!
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("static deflection - consistent pressure load");
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat_lumped, scale_def / max(norm(sol_stat_lumped.def(:, 1:3), "rows")));
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("static deflection - lumped pressure load");

%! X1 = mesh.nodes(unique(elno_p1), 1:3).';
%! X2 = mesh.nodes(unique(elno_p2), 1:3).';
%! X3 = mesh.nodes(unique(elno_p3), 1:3).';
%! dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%! dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%! dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%! F1_con = full(mat_ass.R(dof1)).';
%! F2_con = full(mat_ass.R(dof2)).';
%! F3_con = full(mat_ass.R(dof3)).';
%! M1_con = cross(X1, F1_con);
%! M2_con = cross(X2, F2_con);
%! M3_con = cross(X3, F3_con);
%! F1_lumped = full(mat_ass.Rlumped(dof1)).';
%! F2_lumped = full(mat_ass.Rlumped(dof2)).';
%! F3_lumped = full(mat_ass.Rlumped(dof3)).';
%! M1_lumped = cross(X1, F1_lumped);
%! M2_lumped = cross(X2, F2_lumped);
%! M3_lumped = cross(X3, F3_lumped);
%! Ftot1_con = sum(F1_con, 2);
%! Ftot2_con = sum(F2_con, 2);
%! Ftot3_con = sum(F3_con, 2);
%! Mtot1_con = sum(M1_con, 2);
%! Mtot2_con = sum(M2_con, 2);
%! Mtot3_con = sum(M3_con, 2);
%!
%! Ftot1_lumped = sum(F1_lumped, 2);
%! Ftot2_lumped = sum(F2_lumped, 2);
%! Ftot3_lumped = sum(F3_lumped, 2);
%! Mtot1_lumped = sum(M1_lumped, 2);
%! Mtot2_lumped = sum(M2_lumped, 2);
%! Mtot3_lumped = sum(M3_lumped, 2);

%! ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%! F1_an = [ri * b * p1;
%!          ri * b * p1;
%!          0];

%! M1_an = [-ri * b * p1 * (c + b/2);
%!           ri * b * p1 * (c + b/2);
%!           0];

%! F2_an = [-ro * b * p2;
%!          -ro * b * p2;
%!           0];

%! M2_an = [ ro * b * p2 * (c + b/2);
%!          -ro * b * p2 * (c + b/2);
%!           0];

%! F3_an = [0;
%!          0;
%!          -p3 * (ro^2 - ri^2) * pi / 4];

%! M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!           ys * p3 * (ro^2 - ri^2) * pi / 4;
%!          0];

%! assert(Ftot1_con, F1_an, eps^0.9 * norm(F1_an));
%! assert(Ftot2_con, F2_an, eps^0.9 * norm(F2_an));
%! assert(Ftot1_lumped, F1_an, eps^0.9 * norm(F1_an));
%! assert(Ftot2_lumped, F2_an, eps^0.9 * norm(F2_an));

%! assert(Mtot1_con, M1_an, eps^0.9 * norm(M1_an));
%! assert(Mtot2_con, M2_an, eps^0.9 * norm(M2_an));
%! assert(Mtot1_lumped, M1_an, eps^0.2 * norm(M1_an));
%! assert(Mtot2_lumped, M2_an, eps^0.2 * norm(M2_an));

%! assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%! assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%! assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%! assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ## TEST3
%! close all;
%! number_of_modes = 3;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! plot_modes = true;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif

%! a = 30e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = -5e-3;
%! e = 35e-3;
%! h = 4e-3;
%! alpha = 1e-6;
%! beta = 1e-4;
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"modal\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"interfaces\",2) = {tmp[2]};\n");
%! fputs(fd, "//+\n");
%! fputs(fd, "SetFactory(\"Built-in\");\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-ho_min", "0.95", "-ho_max", "1.05", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%! cms_opt.nodes.interfaces.number = int32(rows(mesh.nodes) + 2);
%! cms_opt.invariants = false;
%! cms_opt.algorithm = "unsymmetric";
%! mesh.nodes(cms_opt.nodes.modal.number, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%! mesh.nodes([cms_opt.nodes.interfaces.number], :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! fprintf(stderr, "building rbe3 elements ...\n");
%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                [1, 2], ...
%!                                                [cms_opt.nodes.modal.number, ...
%!                                                 cms_opt.nodes.interfaces.number]);
%! for i=1:numel(mesh.elements.rbe3)
%!   assert(sum(mesh.elements.rbe3(i).weight), b * c, sqrt(eps) * (b * c));
%! endfor
%! fprintf(stderr, "building cms element ...\n");
%! [mesh, mat_ass_cms, dof_map_cms, sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%! mat_ass_cms.Dred = alpha * mat_ass_cms.Mred + beta * mat_ass_cms.Kred;
%!
%! mbdyn_pre_write_fem_data([filename, ".fem"], mat_ass_cms.Mred, mat_ass_cms.Dred, mat_ass_cms.Kred, mat_ass_cms.Phi, mesh.nodes(:, 1:3).', zeros(columns(mat_ass_cms.Mred), 1), zeros(columns(mat_ass_cms.Mred), 1), mat_ass_cms.diagM);
%! fprintf(stderr, "file \"%s\" created ...\n", [filename, ".fem"]);
%! if (plot_modes)
%! fprintf(stderr, "plotting ...\n");
%! figure("visible","off");
%! hold on;
%! fem_post_sol_plot(mesh);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('undeformed mesh');

%! for i=1:min(number_of_modes, numel(sol_eig_cms.f))
%!      opt_plot.elem_types = {"tria6"};
%!      figure("visible", "off");
%!      hold on;
%!      fem_post_sol_plot(mesh, sol_eig_cms, scale_eig / max(norm(sol_eig_cms.def(:, 1:3, i), "rows")), i, opt_plot);
%!      view(30,30);
%!      xlabel('x [m]');
%!      ylabel('y [m]');
%!      zlabel('z [m]');
%!      grid on;
%!      grid minor on;
%!      title(sprintf("%d. eigenmode: %gHz", i, sol_eig_cms.f(i)));
%! endfor
%! figure_list();
%! endif
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ## TEST4
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! p1 = 25.79e6;
%! p2 = 7.83e6;
%! p3 = 1.3758e6;
%! scale_def = 5e-3;
%! mesh_size = 1e-3;
%! enable_linear_dist = false;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! grp_id_p3 = find([[mesh.groups.tria6].id] == 4);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elem_id_p2 = mesh.groups.tria6(grp_id_p2).elements;
%! elem_id_p3 = mesh.groups.tria6(grp_id_p3).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%! elno_p2 = mesh.elements.tria6(elem_id_p2, :);
%! elno_p3 = mesh.elements.tria6(elem_id_p3, :);
%! x1 = mesh.nodes(:, 1)(elno_p1);
%! y1 = mesh.nodes(:, 2)(elno_p1);
%! z1 = mesh.nodes(:, 3)(elno_p1);
%! Phi1 = atan2(y1, x1);
%!
%! x2 = mesh.nodes(:, 1)(elno_p2);
%! y2 = mesh.nodes(:, 2)(elno_p2);
%! z2 = mesh.nodes(:, 3)(elno_p2);
%! Phi2 = atan2(y2, x2);

%! load_case.pressure.tria6.elements = [elno_p1; elno_p2; elno_p3];
%! load_case.pressure.tria6.p = [p1 * Phi1 / (pi / 2) .* z1 / h;
%!                               p2 * Phi2 / (pi / 2) .* z2 / h;
%!                               repmat(p3, rows(elno_p3), columns(elno_p3))];
%! if (enable_linear_dist)
%! p_mid = load_case.pressure.tria6.p(:, 1:3);
%! load_case.pressure.tria6.p(:, 4:6) = [0.5 * (p_mid(:, 1) + p_mid(:, 2)), 0.5 * (p_mid(:, 2) + p_mid(:, 3)), 0.5 * (p_mid(:, 1) + p_mid(:, 3))];
%! endif
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS, ...
%!                                     FEM_VEC_LOAD_CONSISTENT, ...
%!                                     FEM_VEC_LOAD_LUMPED], ...
%!                                    load_case);
%!
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("static deflection - consistent pressure load");
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat_lumped, scale_def / max(norm(sol_stat_lumped.def(:, 1:3), "rows")));
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("static deflection - lumped pressure load");

%! X1 = mesh.nodes(unique(elno_p1), 1:3).';
%! X2 = mesh.nodes(unique(elno_p2), 1:3).';
%! X3 = mesh.nodes(unique(elno_p3), 1:3).';
%! dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%! dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%! dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%! F1_con = full(mat_ass.R(dof1)).';
%! F2_con = full(mat_ass.R(dof2)).';
%! F3_con = full(mat_ass.R(dof3)).';
%! M1_con = cross(X1, F1_con);
%! M2_con = cross(X2, F2_con);
%! M3_con = cross(X3, F3_con);
%! F1_lumped = full(mat_ass.Rlumped(dof1)).';
%! F2_lumped = full(mat_ass.Rlumped(dof2)).';
%! F3_lumped = full(mat_ass.Rlumped(dof3)).';
%! M1_lumped = cross(X1, F1_lumped);
%! M2_lumped = cross(X2, F2_lumped);
%! M3_lumped = cross(X3, F3_lumped);
%! Ftot1_con = sum(F1_con, 2);
%! Ftot2_con = sum(F2_con, 2);
%! Ftot3_con = sum(F3_con, 2);
%! Mtot1_con = sum(M1_con, 2);
%! Mtot2_con = sum(M2_con, 2);
%! Mtot3_con = sum(M3_con, 2);
%!
%! Ftot1_lumped = sum(F1_lumped, 2);
%! Ftot2_lumped = sum(F2_lumped, 2);
%! Ftot3_lumped = sum(F3_lumped, 2);
%! Mtot1_lumped = sum(M1_lumped, 2);
%! Mtot2_lumped = sum(M2_lumped, 2);
%! Mtot3_lumped = sum(M3_lumped, 2);

%! F1_an = [(2*(pi/2-1)*((h^2-2*c*h+c^2)/2-c^2/2)*p1*ri)/(pi*h);
%!          (2*((h^2-2*c*h+c^2)/2-c^2/2)*p1*ri)/(pi*h);
%!          0];

%! F2_an = [-(2*(pi/2-1)*((h^2-2*c*h+c^2)/2-c^2/2)*p2*ro)/(pi*h);
%!          -(2*((h^2-2*c*h+c^2)/2-c^2/2)*p2*ro)/(pi*h);
%!          0];

%! M1_an = [-(2*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p1*ri)/(pi*h);
%!          (2*(pi/2-1)*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p1*ri)/(pi*h);
%!          0];

%! M2_an = [(2*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p2*ro)/(pi*h);
%!          -(2*(pi/2-1)*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p2*ro)/(pi*h);
%!          0];

%! ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%! F3_an = [0;
%!          0;
%!          -p3 * (ro^2 - ri^2) * pi / 4];

%! M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!           ys * p3 * (ro^2 - ri^2) * pi / 4;
%!          0];

%! assert(Ftot1_con, F1_an, 1e-4 * norm(F1_an));
%! assert(Ftot2_con, F2_an, 1e-4 * norm(F2_an));
%! assert(Ftot1_lumped, F1_an, 2e-3 * norm(F1_an));
%! assert(Ftot2_lumped, F2_an, 2e-3 * norm(F2_an));

%! assert(Mtot1_con, M1_an, 1e-4 * norm(M1_an));
%! assert(Mtot2_con, M2_an, 1e-4 * norm(M2_an));
%! assert(Mtot1_lumped, M1_an, 5e-3 * norm(M1_an));
%! assert(Mtot2_lumped, M2_an, 5e-3 * norm(M2_an));

%! assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%! assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%! assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%! assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));

%! fprintf(stderr, "err(F1_con)=%e\n", norm(Ftot1_con - F1_an) / norm(F1_an));
%! fprintf(stderr, "err(F2_con)=%e\n", norm(Ftot2_con - F2_an) / norm(F2_an));
%! fprintf(stderr, "err(F3_con)=%e\n", norm(Ftot3_con - F3_an) / norm(F3_an));
%! fprintf(stderr, "err(M1_con)=%e\n", norm(Mtot1_con - M1_an) / norm(M1_an));
%! fprintf(stderr, "err(M2_con)=%e\n", norm(Mtot2_con - M2_an) / norm(M2_an));
%! fprintf(stderr, "err(M3_con)=%e\n", norm(Mtot3_con - M3_an) / norm(M3_an));
%! fprintf(stderr, "err(F1_lumped)=%e\n", norm(Ftot1_lumped - F1_an) / norm(F1_an));
%! fprintf(stderr, "err(F2_lumped)=%e\n", norm(Ftot2_lumped - F2_an) / norm(F2_an));
%! fprintf(stderr, "err(F3_lumped)=%e\n", norm(Ftot3_lumped - F3_an) / norm(F3_an));
%! fprintf(stderr, "err(M1_lumped)=%e\n", norm(Mtot1_lumped - M1_an) / norm(M1_an));
%! fprintf(stderr, "err(M2_lumped)=%e\n", norm(Mtot2_lumped - M2_an) / norm(M2_an));
%! fprintf(stderr, "err(M3_lumped)=%e\n", norm(Mtot3_lumped - M3_an) / norm(M3_an));

%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ## TEST5
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! p1 = 25.79e6;
%! p2 = 7.83e6;
%! p3 = 1.3758e6;
%! scale_def = 100;
%! mesh_size = 2e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! grp_id_p3 = find([[mesh.groups.tria6].id] == 4);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elem_id_p2 = mesh.groups.tria6(grp_id_p2).elements;
%! elem_id_p3 = mesh.groups.tria6(grp_id_p3).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%! elno_p2 = mesh.elements.tria6(elem_id_p2, :);
%! elno_p3 = mesh.elements.tria6(elem_id_p3, :);
%! x1 = mesh.nodes(:, 1)(elno_p1);
%! y1 = mesh.nodes(:, 2)(elno_p1);
%! z1 = mesh.nodes(:, 3)(elno_p1);
%! Phi1 = atan2(y1, x1);
%!
%! x2 = mesh.nodes(:, 1)(elno_p2);
%! y2 = mesh.nodes(:, 2)(elno_p2);
%! z2 = mesh.nodes(:, 3)(elno_p2);
%! Phi2 = atan2(y2, x2);

%! load_case.pressure.tria6.elements = [elno_p1; elno_p2; elno_p3];
%! load_case.pressure.tria6.p = [p1 * cos(Phi1) .* sin(pi * z1 / h);
%!                               p2 * cos(Phi2) .* sin(pi * z2 / h);
%!                               repmat(p3, rows(elno_p3), columns(elno_p3))];

%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS, ...
%!                                     FEM_VEC_LOAD_CONSISTENT, ...
%!                                     FEM_VEC_LOAD_LUMPED], ...
%!                                    load_case);
%!
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat, scale_def);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("static deflection - consistent pressure load");
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat_lumped, scale_def);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("static deflection - lumped pressure load");

%! X1 = mesh.nodes(unique(elno_p1), 1:3).';
%! X2 = mesh.nodes(unique(elno_p2), 1:3).';
%! X3 = mesh.nodes(unique(elno_p3), 1:3).';
%! dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%! dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%! dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%! F1_con = full(mat_ass.R(dof1)).';
%! F2_con = full(mat_ass.R(dof2)).';
%! F3_con = full(mat_ass.R(dof3)).';
%! M1_con = cross(X1, F1_con);
%! M2_con = cross(X2, F2_con);
%! M3_con = cross(X3, F3_con);
%! F1_lumped = full(mat_ass.Rlumped(dof1)).';
%! F2_lumped = full(mat_ass.Rlumped(dof2)).';
%! F3_lumped = full(mat_ass.Rlumped(dof3)).';
%! M1_lumped = cross(X1, F1_lumped);
%! M2_lumped = cross(X2, F2_lumped);
%! M3_lumped = cross(X3, F3_lumped);
%! Ftot1_con = sum(F1_con, 2);
%! Ftot2_con = sum(F2_con, 2);
%! Ftot3_con = sum(F3_con, 2);
%! Mtot1_con = sum(M1_con, 2);
%! Mtot2_con = sum(M2_con, 2);
%! Mtot3_con = sum(M3_con, 2);
%!
%! Ftot1_lumped = sum(F1_lumped, 2);
%! Ftot2_lumped = sum(F2_lumped, 2);
%! Ftot3_lumped = sum(F3_lumped, 2);
%! Mtot1_lumped = sum(M1_lumped, 2);
%! Mtot2_lumped = sum(M2_lumped, 2);
%! Mtot3_lumped = sum(M3_lumped, 2);

%! F1_an = [(pi*((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p1*ri)/4;
%!          (((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p1*ri)/2;
%!          0];
%! F2_an = [-(pi*((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p2*ro)/4;
%!          -(((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p2*ro)/2;
%!          0];
%! M1_an = [-(((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p1*ri)/2;
%!          (pi*((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p1*ri)/4;
%!          0];
%! M2_an = [(((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p2*ro)/2;
%!          -(pi*((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p2*ro)/4;
%!          0];
%! ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%! F3_an = [0;
%!          0;
%!          -p3 * (ro^2 - ri^2) * pi / 4];

%! M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!           ys * p3 * (ro^2 - ri^2) * pi / 4;
%!          0];

%! fprintf(stderr, "err(F1_con)=%e\n", norm(Ftot1_con - F1_an) / norm(F1_an));
%! fprintf(stderr, "err(F2_con)=%e\n", norm(Ftot2_con - F2_an) / norm(F2_an));
%! fprintf(stderr, "err(F3_con)=%e\n", norm(Ftot3_con - F3_an) / norm(F3_an));
%! fprintf(stderr, "err(M1_con)=%e\n", norm(Mtot1_con - M1_an) / norm(M1_an));
%! fprintf(stderr, "err(M2_con)=%e\n", norm(Mtot2_con - M2_an) / norm(M2_an));
%! fprintf(stderr, "err(M3_con)=%e\n", norm(Mtot3_con - M3_an) / norm(M3_an));
%! fprintf(stderr, "err(F1_lumped)=%e\n", norm(Ftot1_lumped - F1_an) / norm(F1_an));
%! fprintf(stderr, "err(F2_lumped)=%e\n", norm(Ftot2_lumped - F2_an) / norm(F2_an));
%! fprintf(stderr, "err(F3_lumped)=%e\n", norm(Ftot3_lumped - F3_an) / norm(F3_an));
%! fprintf(stderr, "err(M1_lumped)=%e\n", norm(Mtot1_lumped - M1_an) / norm(M1_an));
%! fprintf(stderr, "err(M2_lumped)=%e\n", norm(Mtot2_lumped - M2_an) / norm(M2_an));
%! fprintf(stderr, "err(M3_lumped)=%e\n", norm(Mtot3_lumped - M3_an) / norm(M3_an));

%! assert(Ftot1_con, F1_an, eps^0.3 * norm(F1_an));
%! assert(Ftot2_con, F2_an, eps^0.3 * norm(F2_an));
%! assert(Ftot1_lumped, F1_an, 2e-2 * norm(F1_an));
%! assert(Ftot2_lumped, F2_an, 2e-2 * norm(F2_an));

%! assert(Mtot1_con, M1_an, eps^0.3 * norm(M1_an));
%! assert(Mtot2_con, M2_an, eps^0.3 * norm(M2_an));
%! assert(Mtot1_lumped, M1_an, 1e-2 * norm(M1_an));
%! assert(Mtot2_lumped, M2_an, 1e-2 * norm(M2_an));

%! assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%! assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%! assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%! assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));

%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST6
%! close all;
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 15e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 3.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! if (do_rotate)
%! fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! endif
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh1.nodes));
%! unlink([filename, ".msh"]);
%! if (~do_rotate)
%! group_defs(1).id = 1;
%! group_defs(1).name = "box1";
%! group_defs(1).R = eye(3);
%! group_defs(1).X0 = zeros(3, 1);
%! group_defs(1).type = "box";
%! group_defs(1).geometry.xmin = 0;
%! group_defs(1).geometry.xmax = 0;
%! group_defs(1).geometry.ymin = 0;
%! group_defs(1).geometry.ymax = b;
%! group_defs(1).geometry.zmin = 0;
%! group_defs(1).geometry.zmax = c;

%! group_defs(2).id = 2;
%! group_defs(2).name = "cylinder1";
%! group_defs(2).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%! group_defs(2).type = "cylinder";
%! group_defs(2).geometry.rmin = 0;
%! group_defs(2).geometry.rmax = 0.5 * c;
%! group_defs(2).geometry.zmin = -0.5 * b;
%! group_defs(2).geometry.zmax = 0.5 * b;

%! group_defs(3).id = 3;
%! group_defs(3).name = "cylinder2";
%! group_defs(3).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(3).type = "cylinder";
%! group_defs(3).geometry.rmin = 0;
%! group_defs(3).geometry.rmax = 0.5 * c;
%! group_defs(3).geometry.zmin = -0.5 * b;
%! group_defs(3).geometry.zmax = 0.5 * b;

%! group_defs(4).id = 4;
%! group_defs(4).name = "box2";
%! group_defs(4).R = eye(3);
%! group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(4).type = "box";
%! group_defs(4).geometry.xmin = 0;
%! group_defs(4).geometry.xmax = 0;
%! group_defs(4).geometry.ymin = -0.5 * b;
%! group_defs(4).geometry.ymax = 0.5 * b;
%! group_defs(4).geometry.zmin = -0.5 * c;
%! group_defs(4).geometry.zmax = 0.5 * c;
%! endif
%! fprintf(stderr, "assembling matrices ...\n");
%! mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%! mesh1.material_data.E = 210000e6;
%! mesh1.material_data.nu = 0.3;
%! mesh1.material_data.rho = 7850;
%! mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%! mesh2 = mesh1;
%! mesh2.nodes(:, 1) += a;
%! for i=1:numel(mesh2.groups.tria6)
%!   mesh2.groups.tria6(i).id += 100;
%! endfor
%! data(1).mesh = mesh1;
%! data(2).mesh = mesh2;
%! [mesh] = fem_post_mesh_merge(data);
%! mesh.elements.sfncon6.slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 2)).nodes(:);
%! mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 101)).elements, :);
%! mesh.elements.sfncon6.maxdist = sqrt(eps) * max(abs([a,b,c]));
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mtot] = fem_ass_matrix(mesh, ...
%!                         dof_map, ...
%!                         [FEM_MAT_STIFFNESS, ...
%!                          FEM_MAT_MASS, ...
%!                          FEM_SCA_TOT_MASS], load_case);
%! assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%! fprintf(stderr, "eigenanalysis ...\n");
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%! [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_SCA_STRESS_VMIS], ...
%!                                   load_case, ...
%!                                   sol_eig);

%! for i=1:numel(sol_eig.f)
%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%! endfor

%! f_ref = [8768.74;
%!          14636.1;
%!          21145.7;
%!          39712.8;
%!          43555.5;
%!          47909;
%!          62270.4;
%!          84324.4;
%!          92665.1;
%!          94563];
%! for i=1:length(sol_eig)
%!  fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig(i).f);
%! endfor
%! assert(sol_eig.f(:), f_ref, tol * max(f_ref));
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST7
%! close all;
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 10e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 3.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! if (do_rotate)
%! fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! endif
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh1.nodes));
%! unlink([filename, ".msh"]);
%! if (~do_rotate)
%! group_defs(1).id = 1;
%! group_defs(1).name = "box1";
%! group_defs(1).R = eye(3);
%! group_defs(1).X0 = zeros(3, 1);
%! group_defs(1).type = "box";
%! group_defs(1).geometry.xmin = 0;
%! group_defs(1).geometry.xmax = 0;
%! group_defs(1).geometry.ymin = 0;
%! group_defs(1).geometry.ymax = b;
%! group_defs(1).geometry.zmin = 0;
%! group_defs(1).geometry.zmax = c;

%! group_defs(2).id = 2;
%! group_defs(2).name = "cylinder1";
%! group_defs(2).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%! group_defs(2).type = "cylinder";
%! group_defs(2).geometry.rmin = 0;
%! group_defs(2).geometry.rmax = 0.5 * c;
%! group_defs(2).geometry.zmin = -0.5 * b;
%! group_defs(2).geometry.zmax = 0.5 * b;

%! group_defs(3).id = 3;
%! group_defs(3).name = "cylinder2";
%! group_defs(3).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(3).type = "cylinder";
%! group_defs(3).geometry.rmin = 0;
%! group_defs(3).geometry.rmax = 0.5 * c;
%! group_defs(3).geometry.zmin = -0.5 * b;
%! group_defs(3).geometry.zmax = 0.5 * b;

%! group_defs(4).id = 4;
%! group_defs(4).name = "box2";
%! group_defs(4).R = eye(3);
%! group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(4).type = "box";
%! group_defs(4).geometry.xmin = 0;
%! group_defs(4).geometry.xmax = 0;
%! group_defs(4).geometry.ymin = -0.5 * b;
%! group_defs(4).geometry.ymax = 0.5 * b;
%! group_defs(4).geometry.zmin = -0.5 * c;
%! group_defs(4).geometry.zmax = 0.5 * c;
%! endif
%! fprintf(stderr, "assembling matrices ...\n");
%! mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%! mesh1.material_data.E = 210000e6;
%! mesh1.material_data.nu = 0.3;
%! mesh1.material_data.rho = 7850;
%! mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%! data.mesh = mesh1;
%! data = repmat(data, 1, 3);
%! for i=2:3
%!  data(i).mesh.nodes(:, 1) += a * (i - 1);
%!  for j=1:numel(data(i).mesh.groups.tria6)
%!    data(i).mesh.groups.tria6(j).id += 100 * (i - 1);
%!  endfor
%! endfor
%! [mesh] = fem_post_mesh_merge(data);
%!
%! for i=1:2
%!   grp_idx_slave = find([[mesh.groups.tria6].id] == (i - 1) * 100 + 2);
%!   grp_idx_master = mesh.groups.tria6(find([[mesh.groups.tria6].id] == i * 100 + 1)).elements;
%!   mesh.elements.sfncon6(i).slave = mesh.groups.tria6(grp_idx_slave).nodes(:);
%!   mesh.elements.sfncon6(i).master = mesh.elements.tria6(grp_idx_master, :);
%!   mesh.elements.sfncon6(i).maxdist = sqrt(eps) * max(abs([a,b,c]));
%! endfor
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mtot] = fem_ass_matrix(mesh, ...
%!                         dof_map, ...
%!                         [FEM_MAT_STIFFNESS, ...
%!                          FEM_MAT_MASS, ...
%!                          FEM_SCA_TOT_MASS], ...
%!                         load_case);
%! assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%! fprintf(stderr, "eigenanalysis ...\n");
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%! sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_VEC_STRESS_CAUCH], ...
%!                                 load_case, ...
%!                                 sol_eig);
%! for i=1:numel(sol_eig.f)
%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3,i), "rows")), i);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%! endfor
%! f_ref = [8768.74;
%!          14636.1;
%!          21145.7;
%!          39712.8;
%!          43555.5;
%!          47909;
%!          62270.4;
%!          84324.4;
%!          92665.1;
%!          94563];
%! for i=1:length(sol_eig.f)
%!  fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig.f(i));
%! endfor
%! assert(sol_eig.f(:), f_ref, tol * max(f_ref));
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST8
%! close all;
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 10e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 3.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! if (do_rotate)
%! fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! endif
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh1.nodes));
%! unlink([filename, ".msh"]);
%! if (~do_rotate)
%! group_defs(1).id = 1;
%! group_defs(1).name = "box1";
%! group_defs(1).R = eye(3);
%! group_defs(1).X0 = zeros(3, 1);
%! group_defs(1).type = "box";
%! group_defs(1).geometry.xmin = 0;
%! group_defs(1).geometry.xmax = 0;
%! group_defs(1).geometry.ymin = 0;
%! group_defs(1).geometry.ymax = b;
%! group_defs(1).geometry.zmin = 0;
%! group_defs(1).geometry.zmax = c;

%! group_defs(2).id = 2;
%! group_defs(2).name = "cylinder1";
%! group_defs(2).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%! group_defs(2).type = "cylinder";
%! group_defs(2).geometry.rmin = 0;
%! group_defs(2).geometry.rmax = 0.5 * c;
%! group_defs(2).geometry.zmin = -0.5 * b;
%! group_defs(2).geometry.zmax = 0.5 * b;

%! group_defs(3).id = 3;
%! group_defs(3).name = "cylinder2";
%! group_defs(3).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(3).type = "cylinder";
%! group_defs(3).geometry.rmin = 0;
%! group_defs(3).geometry.rmax = 0.5 * c;
%! group_defs(3).geometry.zmin = -0.5 * b;
%! group_defs(3).geometry.zmax = 0.5 * b;

%! group_defs(4).id = 4;
%! group_defs(4).name = "box2";
%! group_defs(4).R = eye(3);
%! group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(4).type = "box";
%! group_defs(4).geometry.xmin = 0;
%! group_defs(4).geometry.xmax = 0;
%! group_defs(4).geometry.ymin = -0.5 * b;
%! group_defs(4).geometry.ymax = 0.5 * b;
%! group_defs(4).geometry.zmin = -0.5 * c;
%! group_defs(4).geometry.zmax = 0.5 * c;
%! endif
%! fprintf(stderr, "assembling matrices ...\n");
%! mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%! mesh1.material_data.E = 210000e6;
%! mesh1.material_data.nu = 0.3;
%! mesh1.material_data.rho = 7850;
%! mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%! data.mesh = mesh1;
%! data = repmat(data, 1, 3);
%! for i=2:3
%!  data(i).mesh.nodes(:, 1) += a * (i - 1);
%!  for j=1:numel(data(i).mesh.groups.tria6)
%!    data(i).mesh.groups.tria6(j).id += 100 * (i - 1);
%!  endfor
%! endfor
%! [mesh] = fem_post_mesh_merge(data);
%!
%! for i=1:2
%!   grp_idx_slave = find([[mesh.groups.tria6].id] == i * 100 + 1);
%!   grp_idx_master = mesh.groups.tria6(find([[mesh.groups.tria6].id] == (i - 1) * 100 + 2)).elements;
%!   mesh.elements.sfncon6(i).slave = mesh.groups.tria6(grp_idx_slave).nodes(:);
%!   mesh.elements.sfncon6(i).master = mesh.elements.tria6(grp_idx_master, :);
%!   mesh.elements.sfncon6(i).maxdist = sqrt(eps) * max(abs([a,b,c]));
%! endfor
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mtot] = fem_ass_matrix(mesh, ...
%!                         dof_map, ...
%!                         [FEM_MAT_STIFFNESS, ...
%!                          FEM_MAT_MASS, ...
%!                          FEM_SCA_TOT_MASS], ...
%!                         load_case);
%! assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%! fprintf(stderr, "eigenanalysis ...\n");
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);

%! for i=1:numel(sol_eig.f)
%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%! endfor
%! f_ref = [8768.74;
%!          14636.1;
%!          21145.7;
%!          39712.8;
%!          43555.5;
%!          47909;
%!          62270.4;
%!          84324.4;
%!          92665.1;
%!          94563];
%! for i=1:length(sol_eig.f)
%!  fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig.f(i));
%! endfor
%! assert(sol_eig.f(:), f_ref, tol * max(f_ref));
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST9
%! close all;
%! number_of_modes = 10;
%! scale_eig = 1e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 15e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 3.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! if (do_rotate)
%! fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%! endif
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh1.nodes));
%! unlink([filename, ".msh"]);
%! if (~do_rotate)
%! group_defs(1).id = 1;
%! group_defs(1).name = "box1";
%! group_defs(1).R = eye(3);
%! group_defs(1).X0 = zeros(3, 1);
%! group_defs(1).type = "box";
%! group_defs(1).geometry.xmin = 0;
%! group_defs(1).geometry.xmax = 0;
%! group_defs(1).geometry.ymin = 0;
%! group_defs(1).geometry.ymax = b;
%! group_defs(1).geometry.zmin = 0;
%! group_defs(1).geometry.zmax = c;

%! group_defs(2).id = 2;
%! group_defs(2).name = "cylinder1";
%! group_defs(2).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(2).X0 = [0; 0.5 * b; 0.5 * c];
%! group_defs(2).type = "cylinder";
%! group_defs(2).geometry.rmin = 0;
%! group_defs(2).geometry.rmax = 0.5 * c;
%! group_defs(2).geometry.zmin = -0.5 * b;
%! group_defs(2).geometry.zmax = 0.5 * b;

%! group_defs(3).id = 3;
%! group_defs(3).name = "cylinder2";
%! group_defs(3).R = [-1, 0, 0;
%!                     0, 0, 1;
%!                     0, 1, 0];
%! group_defs(3).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(3).type = "cylinder";
%! group_defs(3).geometry.rmin = 0;
%! group_defs(3).geometry.rmax = 0.5 * c;
%! group_defs(3).geometry.zmin = -0.5 * b;
%! group_defs(3).geometry.zmax = 0.5 * b;

%! group_defs(4).id = 4;
%! group_defs(4).name = "box2";
%! group_defs(4).R = eye(3);
%! group_defs(4).X0 = [a; 0.5 * b; 0.5 * c];
%! group_defs(4).type = "box";
%! group_defs(4).geometry.xmin = 0;
%! group_defs(4).geometry.xmax = 0;
%! group_defs(4).geometry.ymin = -0.5 * b;
%! group_defs(4).geometry.ymax = 0.5 * b;
%! group_defs(4).geometry.zmin = -0.5 * c;
%! group_defs(4).geometry.zmax = 0.5 * c;
%! endif
%! fprintf(stderr, "assembling matrices ...\n");
%! mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%! mesh1.material_data.E = 210000e6;
%! mesh1.material_data.nu = 0.3;
%! mesh1.material_data.rho = 7850;
%! mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%! mesh2 = mesh1;
%! mesh2.nodes(:, 1) += a;
%! for i=1:numel(mesh2.groups.tria6)
%!   mesh2.groups.tria6(i).id += 100;
%! endfor
%! data(1).mesh = mesh1;
%! data(2).mesh = mesh2;
%! [mesh] = fem_post_mesh_merge(data);
%! mesh.elements.sfncon6.slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 2)).nodes(:);
%! mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 101)).elements, :);
%! mesh.elements.sfncon6.maxdist = sqrt(eps) * max(abs([a,b,c]));
%! mesh.elements.sfncon6.constraint = FEM_CT_SLIDING;
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 102)).nodes, :) = true;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mtot] = fem_ass_matrix(mesh, ...
%!                         dof_map, ...
%!                         [FEM_MAT_STIFFNESS, ...
%!                          FEM_MAT_MASS, ...
%!                          FEM_SCA_TOT_MASS], ...
%!                         load_case);
%! assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%! fprintf(stderr, "eigenanalysis ...\n");
%! [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);

%! for i=1:numel(sol_eig.f)
%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%! endfor
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect


%!test
%! ## TEST10
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! num_modes = 6;
%! shift = 0;
%! scale_eig = 1.5e-3;
%! E = [210000e6; 70000e6];
%! nu = [0.3, 0.3];
%! rho = [7850; 2700];
%! ri = [4e-3; 2.5e-3];
%! ro = [5e-3; ri(1)];
%! h = 12e-3;
%! scale_def = 5e-3;
%! toldist = 1e-3;
%! mesh_size = linspace(1.25e-3, 0.4e-3, 10)(end);
%! num_nodes = zeros(1, numel(mesh_size));
%! assert(numel(ri), numel(ro));
%! f = nan(numel(mesh_size), num_modes);
%! for m=1:numel(mesh_size)
%! clear data;
%! clear mesh;
%! clear mat_ass;
%! clear sol_eig;
%! for i=1:numel(ri)
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri(i));
%! fprintf(fd, "ro = %g;\n", ro(i));
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(4) = {2,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(8) = {6,1};\n");
%! fputs(fd, "Line Loop(5) = {1,4,5,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%! fprintf(fd, "Physical Volume(\"volume\",%d) = {tmp[1]};\n", (i - 1) * 100 + 1);
%! fprintf(fd, "Physical Surface(\"bottom\",%d) = {tmp[2]};\n", (i - 1) * 100 + 1);
%! fprintf(fd, "Physical Surface(\"outside\",%d) = {tmp[3]};\n", (i - 1) * 100 + 2);
%! fprintf(fd, "Physical Surface(\"inside\",%d) = {tmp[5]};\n", (i - 1) * 100 + 3);
%! fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[4]};\n", (i - 1) * 100 + 4);
%! fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[0]};\n", (i - 1) * 100 + 5); ##x
%! fprintf(fd, "Physical Surface(\"right\",%d) = {6};\n", (i - 1) * 100 + 6);  ##y
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-clmin", sprintf("%g", 0.75 * mesh_size(m)), "-clmax", sprintf("%g", 1.25 *mesh_size(m)), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(data(i).mesh.nodes));
%! unlink([filename, ".msh"]);
%! endfor
%! for i=1:numel(data)
%! data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%! data(i).mesh.material_data.rho = rho(i);
%! data(i).mesh.material_data.C = fem_pre_mat_isotropic(E(i), nu(i));
%! endfor
%! mesh = fem_post_mesh_merge(data);
%! mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([mesh.groups.tria6.id]==3)).elements, :);
%! mesh.elements.sfncon6.slave = mesh.groups.tria6(find([mesh.groups.tria6.id]==102)).nodes(:);
%! mesh.elements.sfncon6.maxdist = toldist * ri(1);
%! mesh.elements.sfncon6.constraint = FEM_CT_SLIDING;
%! group_id = [[mesh.groups.tria6].id];
%! node_constr1 = [mesh.groups.tria6(find((group_id == 1)|(group_id == 101))).nodes];
%! node_constr5 = [mesh.groups.tria6(find(((group_id == 5)|(group_id == 105)))).nodes];
%! node_constr6 = [mesh.groups.tria6(find(((group_id == 6)|(group_id == 106)))).nodes];
%! node_constr = [node_constr1, node_constr5, node_constr6];
%! mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%! idx_joint = int32(0);
%! for i=1:numel(node_constr1)
%!   mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_constr1(i);
%! endfor
%! for i=1:numel(node_constr5)
%!   mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_constr5(i);
%! endfor
%! for i=1:numel(node_constr6)
%!   mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_constr6(i);
%! endfor
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! fprintf(stderr, "assembling matrices ...\n");
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! num_nodes(m) = rows(mesh.nodes);
%! [mat_ass.K, ...
%!  mat_ass.M, ...
%!  mat_ass.mtot] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS, ...
%!                                  FEM_MAT_MASS, ...
%!                                  FEM_SCA_TOT_MASS], ...
%!                                 load_case);
%! fprintf(stderr, "eigenanalysis ...\n");
%! sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, shift);
%! f(m, :) = sol_eig.f;
%! fprintf(stderr, "plotting ...\n");
%! figure("visible", "off");
%! for i=1:numel(data)
%! fem_post_sol_plot(data(i).mesh);
%! endfor
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("undeformed mesh");

%! fref = [52248
%!         65901
%!         105250
%!         113070
%!         132630
%!         201780];
%! assert(sol_eig.f(:), fref, 0.1e-2 * max(fref));

%! for i=1:numel(sol_eig.f)
%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%! endfor
%! endfor
%! figure("visible", "off");
%! hold on;
%! for i=1:columns(f)
%! plot(num_nodes.', f(:, i), sprintf("-;mode %d;%d", i, i));
%! endfor
%! xlabel("nodes [1]");
%! ylabel("f [Hz]");
%! grid on;
%! grid minor on;
%! title("natural frequencies versus mesh size");
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ## TEST11
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! num_modes = 6;
%! shift = 0;
%! scale_eig = 1.5e-3;
%! E = [210000e6; 70000e6];
%! nu = [0.3, 0.3];
%! rho = [7850; 2700];
%! ri = [8e-3; 5e-3];
%! ro = [10e-3; ri(1)];
%! h = 12e-3;
%! scale_def = 5e-3;
%! toldist = 1e-3;
%! tolf = 8e-2;
%! mesh_size = linspace(0.6e-3, 0.3e-3, 5)(1);
%! num_nodes = zeros(1, numel(mesh_size));
%! assert(numel(ri), numel(ro));
%! orange = [2,1];
%! rrange = [false,true];
%! f = nan(numel(mesh_size), num_modes, numel(rrange), numel(orange));
%! for o=1:numel(orange)
%! for r=1:numel(rrange)
%!   if (orange(o) == 2 && rrange(r))
%!     continue;
%!   endif
%!   for m=1:numel(mesh_size)
%!     clear data;
%!     clear mesh;
%!     clear mat_ass;
%!     clear sol_eig;
%!     for i=1:numel(ri)
%!       [fd, msg] = fopen([filename, ".geo"], "wt");
%!       if fd == -1
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "ri = %g;\n", ri(i));
%!       fprintf(fd, "ro = %g;\n", ro(i));
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!       fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!       fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!       fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(4) = {2,5};\n");
%!       fputs(fd, "Line(5) = {5,6};\n");
%!       fputs(fd, "Line(8) = {6,1};\n");
%!       fputs(fd, "Line Loop(5) = {1,4,5,8};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       if (orange(o) == 1 || i == 1)
%!       if (rrange(r) || i == 1)
%!         fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{%d}; Recombine; };\n", ceil(0.5 * pi * ri(1) / mesh_size(m)));
%!         fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!       else
%!         fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{%d};};\n", ceil(0.5 * pi * ri(1) / mesh_size(m)));
%!       endif
%!       else
%!         fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6};};\n");
%!       endif
%!       fprintf(fd, "Physical Volume(\"volume\",%d) = {tmp[1]};\n", (i - 1) * 100 + 1);
%!       fprintf(fd, "Physical Surface(\"bottom\",%d) = {tmp[2]};\n", (i - 1) * 100 + 1);
%!       fprintf(fd, "Physical Surface(\"outside\",%d) = {tmp[3]};\n", (i - 1) * 100 + 2);
%!       fprintf(fd, "Physical Surface(\"inside\",%d) = {tmp[5]};\n", (i - 1) * 100 + 3);
%!       fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[4]};\n", (i - 1) * 100 + 4);
%!       fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[0]};\n", (i - 1) * 100 + 5); ##x
%!       fprintf(fd, "Physical Surface(\"right\",%d) = {6};\n", (i - 1) * 100 + 6);  ##y
%!       fclose(fd);
%!       fprintf(stderr, "meshing ...\n");
%!       unlink([filename, ".msh"]);
%!       if (orange(o) == 1 || i == 1)
%!         optargs = {"-order", "1", "-optimize_threshold", "0.95"};
%!       else
%!         optargs = {"-order", "2", "-optimize_ho"};
%!       endif
%!       pid = spawn("gmsh", {"-format", "msh2", "-3", optargs{:}, "-clmin", sprintf("%g", 0.75 * mesh_size(m) * orange(o)), "-clmax", sprintf("%g", 1.25 *mesh_size(m) * orange(o)), [filename, ".geo"]});
%!       status = spawn_wait(pid);
%!       if status ~= 0
%!         warning("gmsh failed with status %d", status);
%!       endif
%!       unlink([filename, ".geo"]);
%!       fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!       data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!       fprintf(stderr, "%d nodes\n", rows(data(i).mesh.nodes));
%!       unlink([filename, ".msh"]);
%!     endfor
%!     for i=1:numel(data)
%!       if (isfield(data(i).mesh.elements, "iso8"))
%!         data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!       else
%!         data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!       endif
%!       data(i).mesh.material_data.rho = rho(i);
%!       data(i).mesh.material_data.C = fem_pre_mat_isotropic(E(i), nu(i));
%!     endfor
%!     mesh = fem_post_mesh_merge(data);
%!     mesh.elements.sfncon4.master = mesh.elements.iso4(mesh.groups.iso4(find([mesh.groups.iso4.id]==3)).elements, :);
%!     if (orange(o) == 1)
%!       mesh.elements.sfncon4.slave = mesh.groups.iso4(find([mesh.groups.iso4.id]==102)).nodes(:);
%!     else
%!       mesh.elements.sfncon4.slave = mesh.groups.tria6(find([mesh.groups.tria6.id]==102)).nodes(:);
%!     endif
%!     mesh.elements.sfncon4.maxdist = toldist * ri(1);
%!     mesh.elements.sfncon4.constraint = FEM_CT_SLIDING;
%!     if (orange(o) == 1)
%!     group_id = [[mesh.groups.iso4].id];
%!     node_constr1 = [mesh.groups.iso4(find((group_id == 1)|(group_id == 101))).nodes];
%!     node_constr5 = [mesh.groups.iso4(find(((group_id == 5)|(group_id == 105)))).nodes];
%!     node_constr6 = [mesh.groups.iso4(find(((group_id == 6)|(group_id == 106)))).nodes];
%!     else
%!     group_id4 = [[mesh.groups.iso4].id];
%!     group_id6 = [[mesh.groups.tria6].id];
%!     node_constr1 = [mesh.groups.iso4(find((group_id4 == 1))).nodes,   mesh.groups.tria6(find((group_id6 == 101))).nodes];
%!     node_constr5 = [mesh.groups.iso4(find(((group_id4 == 5)))).nodes, mesh.groups.tria6(find(((group_id6 == 105)))).nodes];
%!     node_constr6 = [mesh.groups.iso4(find(((group_id4 == 6)))).nodes, mesh.groups.tria6(find(((group_id6 == 106)))).nodes];
%!     endif
%!     node_constr = [node_constr1, node_constr5, node_constr6];
%!     mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!     idx_joint = int32(0);
%!     for i=1:numel(node_constr1)
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_constr1(i);
%!     endfor
%!     for i=1:numel(node_constr5)
%!       mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_constr5(i);
%!     endfor
%!     for i=1:numel(node_constr6)
%!       mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_constr6(i);
%!     endfor
%!     load_case.locked_dof = false(rows(mesh.nodes), 6);
%!     fprintf(stderr, "assembling matrices ...\n");
%!     dof_map = fem_ass_dof_map(mesh, load_case);
%!     num_nodes(m) = rows(mesh.nodes);
%!     [mat_ass.K, ...
%!      mat_ass.M, ...
%!      mat_ass.mtot] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_STIFFNESS, ...
%!                                      FEM_MAT_MASS, ...
%!                                      FEM_SCA_TOT_MASS], ...
%!                                     load_case);
%!     fprintf(stderr, "eigenanalysis ...\n");
%!     sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, shift);
%!     fem_post_mesh_export([filename, ".msh"], mesh);
%!     f(m, :, r, o) = sol_eig.f;
%!     fprintf(stderr, "plotting ...\n");
%!     figure("visible", "off");
%!     for i=1:numel(data)
%!       fem_post_sol_plot(data(i).mesh);
%!     endfor
%!     xlabel("x [m]");
%!     ylabel("y [m]");
%!     zlabel("z [m]");
%!     grid on;
%!     grid minor on;
%!     title("undeformed mesh");

%!     for i=1:numel(sol_eig.f)
%!       figure("visible", "off");
%!       hold on;
%!       fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!       view(30,30);
%!       xlabel('x [m]');
%!       ylabel('y [m]');
%!       zlabel('z [m]');
%!       grid on;
%!       grid minor on;
%!       title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!     endfor
%!   endfor
%!   figure("visible", "off");
%!   hold on;
%!   for i=1:columns(f)
%!     plot(num_nodes.', f(:, i, r, o), sprintf("-;mode %d;%d", i, i));
%!   endfor
%!   xlabel("nodes [1]");
%!   ylabel("f [Hz]");
%!   grid on;
%!   grid minor on;
%!   title("natural frequencies versus mesh size");
%! endfor
%! endfor
%! fref = [26023  59514  91469  1.0372e+05  1.1294e+05  1.154e+05];
%! for o=1:numel(orange)
%! for r=1:numel(rrange)
%! for i=1:rows(f)
%! if (~all(isnan(f(i, :, r, o))))
%! assert(all(abs(f(i, :, r, o) ./ fref - 1) < tolf));
%! endif
%! endfor
%! endfor
%! endfor
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ## TEST12
%! close all;
%! a = 10e-3;
%! b = 20e-3;
%! c = 10e-3;
%! h = 2e-3;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! px = -10e6;
%! py = -2e6;
%! pz = -3e6;
%! tolstat = 4e-2;
%! scale = 30e-3;
%! maxdist = 1e-2 * max([a,b,c]);
%! eliminate = false;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! for i=1:4
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! unwind_protect
%! switch (i)
%! case {1,2}
%! order = 1;
%! optflags = {"-optimize_threshold", "0.95"};
%! otherwise
%! order = 2;
%! optflags = {"-optimize_ho", "-ho_min", "0.5", "-ho_max", "1.5"};
%! endswitch
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h=%g;\n", h * order);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0.0,0.0,c} {\n");
%! switch (i)
%! case 1
%! fprintf(fd, "  Surface{6}; Layers{%d}; Recombine;\n", ceil(c/h));
%! otherwise
%! fputs(fd, "  Surface{6};\n");
%! endswitch
%! fputs(fd, "};\n");
%! switch (i)
%! case 1
%! fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%! endswitch
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fprintf(fd, "Physical Surface(\"right\",%d) = {tmp[2]};\n", 100 * i + 1);
%! fprintf(fd, "Physical Surface(\"rear\",%d) = {tmp[3]};\n", 100 * i + 2);
%! fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[4]};\n", 100 * i + 3);
%! fprintf(fd, "Physical Surface(\"front\",%d) = {tmp[5]};\n", 100 * i + 4);
%! fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[0]};\n", 100 * i + 5);
%! fprintf(fd, "Physical Surface(\"bottom\",%d) = {6};\n", 100 * i + 6);
%! unwind_protect_cleanup
%! fclose(fd);
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", sprintf("%d", order), optflags{:}, [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(data(i).mesh.nodes));
%! unlink([filename, ".msh"]);
%! endfor
%! data(2).mesh.nodes(:, 1) += a;
%! data(3).mesh.nodes(:, 1) += a;
%! data(3).mesh.nodes(:, 2) += b;
%! data(4).mesh.nodes(:, 2) += b;
%! for i=1:numel(data)
%!   if (isfield(data(i).mesh.elements, "iso8"))
%!     data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!   else
%!     data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!   endif
%!   data(i).mesh.material_data.rho = rho;
%!   data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! endfor
%! mesh = fem_post_mesh_merge(data);
%! constr = FEM_CT_FIXED;
%! mesh.elements.sfncon4(1).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 103)).elements, :);
%! mesh.elements.sfncon4(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 401)).nodes(:);
%! mesh.elements.sfncon4(1).maxdist = maxdist;
%! mesh.elements.sfncon4(1).constraint = constr;
%! mesh.elements.sfncon4(2).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 102)).elements, :);
%! mesh.elements.sfncon4(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 204)).nodes(:);
%! mesh.elements.sfncon4(2).maxdist = maxdist;
%! mesh.elements.sfncon4(2).constraint = constr;
%! mesh.elements.sfncon6(1).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 402)).elements, :);
%! mesh.elements.sfncon6(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 304)).nodes(:);
%! mesh.elements.sfncon6(1).maxdist = maxdist;
%! mesh.elements.sfncon6(1).constraint = constr;
%! mesh.elements.sfncon6(2).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 301)).elements, :);
%! mesh.elements.sfncon6(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 203)).nodes(:);
%! mesh.elements.sfncon6(2).maxdist = maxdist;
%! mesh.elements.sfncon6(2).constraint = constr;
%! slave_nodes = [];
%! for i=1:numel(mesh.elements.sfncon4)
%!   idx_slave = [];
%!   for j=1:numel(slave_nodes)
%!     idx_slave = [idx_slave; find(mesh.elements.sfncon4(i).slave == slave_nodes(j))];
%!   endfor
%!   if (numel(idx_slave))
%!     mesh.elements.sfncon4(i).slave(idx_slave) = 0;
%!     mesh.elements.sfncon4(i).slave = mesh.elements.sfncon4(i).slave(find(mesh.elements.sfncon4(i).slave));
%!   endif
%!   slave_nodes = [slave_nodes; mesh.elements.sfncon4(i).slave];
%! endfor

%! for i=1:numel(mesh.elements.sfncon6)
%!   idx_slave = [];
%!   for j=1:numel(slave_nodes)
%!     idx_slave = [idx_slave; find(mesh.elements.sfncon6(i).slave == slave_nodes(j))];
%!   endfor
%!   if (numel(idx_slave))
%!     mesh.elements.sfncon6(i).slave(idx_slave) = 0;
%!     mesh.elements.sfncon6(i).slave = mesh.elements.sfncon6(i).slave(find(mesh.elements.sfncon6(i).slave));
%!   endif
%!   slave_nodes = [slave_nodes; mesh.elements.sfncon6(i).slave];
%! endfor
%! group_id4 = [[mesh.groups.iso4].id];
%! group_id6 = [[mesh.groups.tria6].id];
%! node_bottom = [[mesh.groups.iso4(find(mod(group_id4, 100) == 6))].nodes,   [mesh.groups.tria6(find(mod(group_id6,100)==6))].nodes];
%! node_front = [[mesh.groups.iso4(find(group_id4 == 104))].nodes, [mesh.groups.tria6(find(group_id6 == 404))].nodes];
%! node_right = [mesh.groups.iso4(find((group_id4 == 101) | (group_id4 == 201))).nodes];
%! node_constr = [node_bottom, node_front, node_right];
%! mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%! idx_joint = int32(0);
%! for i=1:numel(node_bottom)
%!   if (~numel(find(slave_nodes == node_bottom(i))))
%!     mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!     mesh.elements.joints(idx_joint).nodes = node_bottom(i);
%!   endif
%! endfor
%! for i=1:numel(node_front)
%!   if (~numel(find(slave_nodes == node_front(i))))
%!   mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_front(i);
%!   endif
%! endfor
%! for i=1:numel(node_right)
%!   if (~numel(find(slave_nodes == node_right(i))))
%!   mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_right(i);
%!   endif
%! endfor

%! mesh.elements.joints = mesh.elements.joints(1:idx_joint);

%! iso4_top = [[mesh.groups.iso4(find(mod(group_id4, 100) == 5))].elements];
%! tria6_top = [[mesh.groups.tria6(find(mod(group_id6, 100) == 5))].elements];
%! iso4_rear = mesh.groups.iso4(find(group_id4 == 202)).elements;
%! tria6_rear = mesh.groups.tria6(find(group_id6 == 302)).elements;
%! tria6_left = [[mesh.groups.tria6(find((group_id6 == 303) | (group_id6 == 403)))].elements];
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.pressure.tria6.elements = mesh.elements.tria6([tria6_left, tria6_rear, tria6_top], :);
%! load_case.pressure.tria6.p = [repmat(py, numel(tria6_left), 6); repmat(px, numel(tria6_rear), 6); repmat(pz, numel(tria6_top), 6)];
%! load_case.pressure.iso4.elements = mesh.elements.iso4([iso4_rear, iso4_top], :);
%! load_case.pressure.iso4.p = [repmat(px, numel(iso4_rear), 4); repmat(pz, numel(iso4_top), 4)];

%! fprintf(stderr, "assembling matrices ...\n");
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_STIFFNESS, ...
%!                                      FEM_VEC_LOAD_CONSISTENT], ...
%!                                     load_case);
%! if (eliminate)
%! fprintf(stderr, "eliminating constraint equations ...\n");
%! [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%! fprintf(stderr, "linear static solution ...\n");
%! opt_ls.refine_max_iter = int32(100);
%! if (fem_sol_check_func("mumps"))
%! opt_ls.matrix_type = MUMPS_MAT_SYM;
%! opt_ls.verbose = MUMPS_VER_ERR;
%! Kfact = fem_fact_mumps(Kred, opt_ls);
%! else
%! Kfact = fem_fact_umfpack(Kred, opt_ls);
%! endif
%! Ured = Kfact \ Rred;
%! fprintf(stderr, "expanding static solution ...\n");
%! sol_stat.def = fem_post_def_nodal(mesh, dof_map, Tred * Ured);
%! else
%! fprintf(stderr, "linear static solution ...\n");
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! endif
%! sol_stat.F = fem_post_def_nodal(mesh, dof_map, mat_ass.R);
%! Ftot = sum(sol_stat.F, 1);
%! figure("visible", "off");
%! fem_post_sol_plot(mesh);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("undeformed mesh");
%! figure("visible", "off");
%! hold on;
%! fem_post_sol_plot(mesh, sol_stat, scale/max(norm(sol_stat.def, "rows")), 1);
%! view(30,30);
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title("deformed mesh");
%! sigma_a = [-px; -py; -pz; zeros(3, 1)];
%! epsilon_a = mesh.material_data(1).C \ sigma_a;
%! U_a = zeros(rows(mesh.nodes), 3);
%! for i=1:3
%!   U_a(:, i) = epsilon_a(i) * mesh.nodes(:, i);
%! endfor
%! Ftot_a = [-2 * b * c * px, -2 * a * c * py, -4 * a * b * pz, zeros(1, 3)];
%! fprintf(stderr, "max(err)=%g\n", max(max(abs(sol_stat.def(:, 1:3) - U_a))) / max(max(abs(U_a))));
%! fprintf(stderr, "mean(err)=%g\n", mean(mean(abs(sol_stat.def(:, 1:3) - U_a))) / max(max(abs(U_a))));
%! assert(sol_stat.def(:, 1:3), U_a, tolstat * max(max(abs(U_a))));
%! assert(Ftot, Ftot_a, sqrt(eps) * norm(Ftot_a));
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST13
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 60e-3;
%! c = 10e-3;
%! d = 0e-3;
%! h = 2e-3;
%! b = h;
%! Fx = 0;
%! Fz = -1000;
%! My = 0;
%! num_iso = 10;
%! do_post_pro = false;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0, -0.5 * b, -0.5 * c, h};\n");
%! fputs(fd, "Point(2) = {  a, -0.5 * b, -0.5 * c, h};\n");
%! fputs(fd, "Point(3) = {  a,  0.5 * b, -0.5 * c, h};\n");
%! fputs(fd, "Point(4) = {0.0,  0.5 * b, -0.5 * c, h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fputs(fd, "};\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! master_node_idx = int32(rows(mesh.nodes) + 1);
%! mesh.nodes(master_node_idx, 1:3) = [a + d, 0, 0];
%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, master_node_idx);
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);

%! fprintf(stderr, "assembling matrices ...\n");
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%! load_case.loaded_nodes = [master_node_idx];
%! load_case.loads = [Fx, 0, Fz, 0, My, 0];
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.dm] = fem_ass_matrix(mesh, ...
%!                               dof_map, ...
%!                               [FEM_MAT_STIFFNESS, ...
%!                                FEM_VEC_LOAD_CONSISTENT, ...
%!                                FEM_SCA_TOT_MASS], ...
%!                               load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_VEC_STRESS_CAUCH], ...
%!                                  load_case, ...
%!                                  sol_stat);
%! x = linspace(0, a, 100);
%! z = linspace(-0.5 * c, 0.5 * c, 50);
%! [xx, zz] = meshgrid(x, z);
%! xtauel = mesh.nodes(:, 1)(mesh.elements.tet10);
%! ytauel = mesh.nodes(:, 1)(mesh.elements.tet10);
%! ztauel = mesh.nodes(:, 3)(mesh.elements.tet10);
%! tauxxel = sol_stat.stress.taum.tet10(:, :, 1);
%! tauxzel = sol_stat.stress.taum.tet10(:, :, 6);
%! tauxx = griddata(xtauel(:), ztauel(:), tauxxel(:), xx, zz);
%! tauxz = griddata(xtauel(:), ztauel(:), tauxzel(:), xx, zz);
%! Iy = b * c^3 / 12;
%! tauxx_a = -Fz / Iy * (a - xx) .* zz + Fx / (b * c) + My / Iy * zz;
%! tauxz_a = 3 / 2 * Fz / (b * c) * (1 - (zz / (0.5 * c)).^2);
%! scale_tauxx = linspace(min(min(tauxx_a)), max(max(tauxx_a)), num_iso + 1);
%! figure("visible", "off");
%! subplot(2, 1, 1);
%! contourf(xx, zz, tauxx_a, scale_tauxx);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxx [Pa]");
%! grid on;
%! grid minor on;
%! subplot(2, 1, 2);
%! contourf(xx, zz, tauxx, scale_tauxx);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxx [Pa]");
%! grid on;
%! grid minor on;
%! scale_tauxz = linspace(min(min(tauxz_a)), max(max(tauxz_a)), num_iso + 1);
%! figure("visible", "off");
%! subplot(2, 1, 1);
%! contourf(xx, zz, tauxz_a, scale_tauxz);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxz [Pa]");
%! grid on;
%! grid minor on;
%! subplot(2, 1, 2);
%! contourf(xx, zz, tauxz, scale_tauxz);
%! daspect([1,1,1]);
%! colormap jet;
%! colorbar;
%! xlabel("x [m]");
%! ylabel("z [m]");
%! title("stress component tauxz [Pa]");
%! grid on;
%! grid minor on;

%! figure_list();
%! if (do_post_pro)
%!   opts.scale_def = 0.3 * a / max(max(abs(sol_stat.def)));
%!   fem_post_sol_external(mesh, sol_stat, opts);
%! endif
%! idx_x = find((xx(:) > 0.1 * a) & (xx(:) < 0.9 * a));
%! assert(tauxx(:)(idx_x), tauxx_a(:)(idx_x), 0.3e-2 * max(max(max(abs(sol_stat.stress.taum.tet10)))));
%! assert(tauxz(:)(idx_x), tauxz_a(:)(idx_x), 0.3e-2 * max(max(max(abs(sol_stat.stress.taum.tet10)))));
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST14
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! N = 3;
%! a = 10e-3;
%! b = 5e-3;
%! c = 5e-3;
%! h = 5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h = %g;\n", h);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {4,3};\n");
%! fputs(fd, "Line(2) = {3,2};\n");
%! fputs(fd, "Line(3) = {2,1};\n");
%! fputs(fd, "Line(4) = {1,4};\n");
%! fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%! fputs(fd, "  Surface{6};\n");
%! fprintf(fd, "  Layers{%d}; Recombine;\n", ceil(c / h));
%! fputs(fd, "};\n");
%! fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh1.nodes));
%! unlink([filename, ".msh"]);
%! fprintf(stderr, "assembling matrices ...\n");
%! mesh1.materials.iso8 = ones(rows(mesh1.elements.iso8), 1, "int32");
%! mesh1.material_data.E = 210000e6;
%! mesh1.material_data.nu = 0.3;
%! mesh1.material_data.rho = 7850;
%! mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%! data.mesh = mesh1;
%! data = repmat(data, 1, 3);
%! for i=2:3
%!  data(i).mesh.nodes(:, 1) += a * (i - 1);
%!  for j=1:numel(data(i).mesh.groups.iso4)
%!    data(i).mesh.groups.iso4(j).id += 100 * (i - 1);
%!    data(i).mesh.groups.iso4(j).name = sprintf("%s[%d]", data(i).mesh.groups.iso4(j).name, data(i).mesh.groups.iso4(j).id);
%!  endfor
%!  for j=1:numel(data(i).mesh.groups.iso8)
%!    data(i).mesh.groups.iso8(j).id += 100 * (i - 1);
%!    data(i).mesh.groups.iso8(j).name = sprintf("%s[%d]", data(i).mesh.groups.iso8(j).name, data(i).mesh.groups.iso8(j).id);
%!  endfor
%! endfor
%! [mesh] = fem_post_mesh_merge(data);
%! mesh.nodes(:, 2) -= 0.5 * b;
%! mesh.nodes(:, 3) -= 0.5 * c;
%! for i=1:N
%!   if (i > 1)
%!     unwind_protect
%!       fem_post_mesh_export([filename, "_in.msh"], data(i - 1).mesh);
%!       pid = spawn("gmsh", {"-refine", "-format", "msh2", "-o", [filename, "_out.msh"], [filename, "_in.msh"]});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         error("gmsh failed with status %d", status);
%!       endif
%!       data(i).mesh = fem_pre_mesh_import([filename, "_out.msh"]);
%!       data(i).mesh.material_data = mesh.material_data;
%!       data(i).mesh.materials.iso8 = zeros(rows(data(i).mesh.elements.iso8), 1, "int32");
%!       for j=1:numel(data(i).mesh.groups.iso8)
%!         data(i).mesh.materials.iso8(data(i).mesh.groups.iso8(j).elements) = j;
%!       endfor
%!     unwind_protect_cleanup
%!       unlink([filename, "_in.msh"]);
%!       unlink([filename, "_out.msh"]);
%!     end_unwind_protect
%!   else
%!     data(i).mesh = mesh;
%!   endif
%!   fprintf(stderr, "%d: %d nodes\n", i, rows(data(i).mesh.nodes));
%!   data(i).h = h / i;
%!   node_idx_rbe3 = int32(rows(data(i).mesh.nodes) + 1);
%!   data(i).mesh.nodes(node_idx_rbe3, 1:3) = [3 * a, 0, 0];
%!   grp_idx = find([data(i).mesh.groups.iso4.id] == 202);
%!   data(i).mesh.elements.rbe3.nodes = [node_idx_rbe3, data(i).mesh.groups.iso4(grp_idx).nodes];
%!   data(i).mesh.elements.rbe3.weight = ones(numel(data(i).mesh.elements.rbe3.nodes) - 1, 1);
%!   for j=1:2
%!     grp_idx_slave = find([[data(i).mesh.groups.iso4].id] == (j - 1) * 100 + 2);
%!     grp_idx_master = data(i).mesh.groups.iso4(find([[data(i).mesh.groups.iso4].id] == j * 100 + 1)).elements;
%!     data(i).mesh.elements.sfncon4(j).slave = data(i).mesh.groups.iso4(grp_idx_slave).nodes(:);
%!     data(i).mesh.elements.sfncon4(j).master = data(i).mesh.elements.iso4(grp_idx_master, :);
%!     data(i).mesh.elements.sfncon4(j).maxdist = sqrt(eps) * max(abs([a,b,c]));
%!   endfor
%!
%!   data(i).load_case.locked_dof = false(rows(data(i).mesh.nodes), 6);
%!   data(i).load_case.locked_dof(data(i).mesh.groups.iso4(find([[data(i).mesh.groups.iso4].id] == 1)).nodes, :) = true;
%!   data(i).load_case.loaded_nodes = node_idx_rbe3;
%!   data(i).load_case.loads = [0,-1000, 0, 0, 0, 0];
%!   fprintf(stderr, "%d: assembling matrices ...\n", i);
%!   data(i).dof_map = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%!   [data(i).mat_ass.K, ...
%!    data(i).mat_ass.R, ...
%!    data(i).mat_ass.mtot] = fem_ass_matrix(data(i).mesh, ...
%!                                           data(i).dof_map, ...
%!                                           [FEM_MAT_STIFFNESS, ...
%!                                            FEM_VEC_LOAD_CONSISTENT, ...
%!                                            FEM_SCA_TOT_MASS], ...
%!                                           data(i).load_case);
%!   assert(data(i).mat_ass.mtot, ...
%!          a * b * c * sum([data(i).mesh.material_data.rho]), ...
%!          sqrt(eps) * a * b * c * sum([data(i).mesh.material_data.rho]));
%!   fprintf(stderr, "static analysis ...\n");
%!   [data(i).sol_stat, data(i).sol_stat.U] = fem_sol_static(data(i).mesh, data(i).dof_map, data(i).mat_ass);
%!   data(i).sol_stat.stress = fem_ass_matrix(data(i).mesh, ...
%!                                            data(i).dof_map, ...
%!                                            [FEM_VEC_STRESS_CAUCH], ...
%!                                             data(i).load_case, ...
%!                                             data(i).sol_stat);
%!   data(i).W = data(i).sol_stat.U.' * data(i).mat_ass.K * data(i).sol_stat.U;
%! endfor
%! figure("visible", "off");
%! loglog([data.h], [data.W], "-x;W(h) [J];1");
%! xlabel("h [m]");
%! ylabel("W [J]");
%! grid on;
%! grid minor on;
%! title("strain energy");
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test
%! ### TEST15
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! h1 = 10e-3;
%! h2 = 2e-3;
%! w = 4e-3;
%! l = 30e-3;
%! h = 2e-3;
%! do_post_pro = false;
%! scale_def = 20e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "h1=%g;\n", h1);
%! fprintf(fd, "h2=%g;\n", h2)
%! fprintf(fd, "w=%g;\n", w);
%! fprintf(fd, "l=%g;\n", l);
%! fprintf(fd, "h=%g;\n", h);
%! fputs(fd, "Point(1) = {0.0, -0.5 * w, -0.5 * h1 - h2, h};\n");
%! fputs(fd, "Point(2) = {0.0, -0.5 * w, -0.5 * h1, h};\n");
%! fputs(fd, "Point(3) = {0.0,  0.5 * w, -0.5 * h1, h};\n");
%! fputs(fd, "Point(4) = {0.0,  0.5 * w, -0.5 * h1 - h2, h};\n");
%! fputs(fd, "Point(5) = {0.0, -0.5 * w, 0.5 * h1, h};\n");
%! fputs(fd, "Point(6) = {0.0, -0.5 * w, 0.5 * h1 + h2, h};\n");
%! fputs(fd, "Point(7) = {0.0,  0.5 * w, 0.5 * h1 + h2, h};\n");
%! fputs(fd, "Point(8) = {0.0,  0.5 * w, 0.5 * h1, h};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,1};\n");
%! fputs(fd, "Line(5) = {2,5};\n");
%! fputs(fd, "Line(6) = {5,8};\n");
%! fputs(fd, "Line(7) = {8,3};\n");
%! fputs(fd, "Line(8) = {3,2};\n");
%! fputs(fd, "Line(9) = {5,6};\n");
%! fputs(fd, "Line(10) = {6,7};\n");
%! fputs(fd, "Line(11) = {7,8};\n");
%! fputs(fd, "Line(12) = {8,5};\n");
%! fputs(fd, "Line Loop(13) = {1,2,3,4};\n");
%! fputs(fd, "Line Loop(14) = {5,6,7,8};\n");
%! fputs(fd, "Line Loop(15) = {9,10,11,12};\n");
%! fputs(fd, "Plane Surface(16) = {13};\n");
%! fputs(fd, "Plane Surface(17) = {14};\n");
%! fputs(fd, "Plane Surface(18) = {15};\n");
%! fputs(fd, "v1[] = Extrude {l,0.0,0.0} {\n");
%! fputs(fd, "  Surface{16};\n");
%! fputs(fd, "};\n");
%! fputs(fd, "v2[] = Extrude {l,0.0,0.0} {\n");
%! fputs(fd, "  Surface{17};\n");
%! fputs(fd, "};\n");
%! fputs(fd, "v3[] = Extrude {l,0.0,0.0} {\n");
%! fputs(fd, "  Surface{18};\n");
%! fputs(fd, "};\n");
%! fputs(fd, "v = newv;\n");
%! fputs(fd, "BooleanUnion(v) = {Volume{v1[1]};}{ Volume{v2[1],v3[1]};};\n");
%! fputs(fd, "Coherence;\n");
%! fputs(fd, "Physical Volume(\"volume1\",1) = {v1[1]};\n");
%! fputs(fd, "Physical Volume(\"volume2\",2) = {v2[1]};\n");
%! fputs(fd, "Physical Volume(\"volume3\",3) = {v3[1]};\n");
%! fputs(fd, "Physical Surface(\"surface1\",1) = {v1[0]};\n");
%! fputs(fd, "Physical Surface(\"surface2\",2) = {v2[0]};\n");
%! fputs(fd, "Physical Surface(\"surface3\",3) = {v3[0]};\n");
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! mesh.material_data(1).E = 80000e6;
%! mesh.material_data(1).nu = 0.4;
%! mesh.material_data(1).rho = 1700;
%! mesh.material_data(1).alpha = 1e-8;
%! mesh.material_data(1).beta = 1e-6;
%! mesh.material_data(2).E = 10000e6;
%! mesh.material_data(2).nu = 0.4;
%! mesh.material_data(2).rho = 500;
%! mesh.material_data(2).alpha = 0;
%! mesh.material_data(2).beta = 0;
%! for i=1:numel(mesh.material_data)
%!   mesh.material_data(i).C = fem_pre_mat_isotropic(mesh.material_data(i).E, mesh.material_data(i).nu);
%! endfor
%! mesh.materials.tet10 = zeros(rows(mesh.elements.tet10),1);
%! mesh.materials.tet10(mesh.groups.tet10(1).elements) = 1;
%! mesh.materials.tet10(mesh.groups.tet10(2).elements) = 2;
%! mesh.materials.tet10(mesh.groups.tet10(3).elements) = 1;
%! load_case.locked_dof = false(size(mesh.nodes));
%! load_case.locked_dof(find(mesh.nodes(:, 1)==0), 1:3) = true;
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.M, ...
%!  mat_ass.D, ...
%!  mat_ass.K, ...
%!  mat_ass.dm, ...
%!  mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_MAT_MASS, ...
%!                                      FEM_MAT_DAMPING, ...
%!                                      FEM_MAT_STIFFNESS, ...
%!                                      FEM_SCA_TOT_MASS], ...
%!                                     load_case);
%! [sol_eig, err] = fem_sol_modal(mesh, dof_map, mat_ass, 10);
%! [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%! opts.scale_def = scale_def / max(max(max(abs(sol_eig.def))));
%! rho1 = mesh.material_data(1).rho;
%! rho2 = mesh.material_data(2).rho;
%! m1 = rho1 * h2 * w * l;
%! m2 = rho2 * h1 * w * l;
%! m = 2 * m1 + m2;
%! assert(mat_ass.dm, m, sqrt(eps) * m);
%! opts.skin_only = true;
%! if (do_post_pro)
%!   fem_post_sol_external(mesh, sol_eig, opts);
%! endif
%! unwind_protect_cleanup
%!   unlink([filename, ".msh"]);
%!   unlink([filename, ".geo"]);
%! end_unwind_protect

%!test
%! ## TEST 16
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ## K.J.Bathe page 328 4.20a
%! mesh_size = 1.5;
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
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);

%! load_case.pressure.tria6.elements = elno_p1;
%! load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT, ...
%!                               FEM_VEC_LOAD_LUMPED], ...
%!                              load_case);
%!
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%! grp_id_displacement = find([[mesh.groups.tria6].id] == 4);
%! elem_id_displacement = mesh.groups.tria6(grp_id_displacement).elements;
%! elno_id_displacement = mesh.elements.tria6(elem_id_displacement, :);
%! delta = mean(sol_stat.def(elno_id_displacement, 3));
%! grp_id_stress = find([[mesh.groups.tria6].id] == 5);
%! elem_id_stress = mesh.groups.tria6(grp_id_stress).elements;
%! elno_id_stress = mesh.elements.tria6(elem_id_stress, :);
%! taum = zeros(6, numel(elno_id_stress));
%! taum_n = zeros(1, numel(elno_id_stress));
%! for i=1:numel(elno_id_stress)
%!   [ridx, cidx] = find(mesh.elements.tet10 == elno_id_stress(i));
%!   for j=1:numel(ridx)
%!     taum(:, i) += reshape(sol_stat.stress.taum.tet10(ridx(j), cidx(j), :), 6, 1);
%!     ++taum_n(i);
%!   endfor
%! endfor
%! taum *= diag(1 ./ taum_n);
%! sigma1_max = 0;
%! for i=1:columns(taum)
%!   TAU = [taum(1, i), taum(4, i), taum(6, i);
%!          taum(4, i), taum(2, i), taum(5, i);
%!          taum(6, i), taum(5, i), taum(3, i)];
%!   sigma1_max = max(sigma1_max, max(eig(TAU)));
%! endfor
%! fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%! fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%! fprintf(stderr, "delta=%.3f [mm]\n", delta);
%! ## K.J.Bathe page 329 4.20b
%! sigma1_max_ref = 0.6056;
%! delta_ref = -1.669;
%! fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%! fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%! assert(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%! assert(delta, delta_ref, 0.04 * abs(delta_ref));
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!demo
%! ## DEMO 1
%! ## K.J.Bathe 2002, page 328 4.20a
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! animate = true;
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! mesh_size = 1;
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
%! fclose(fd);
%! fprintf(stderr, "meshing ...\n");
%! unlink([filename, ".msh"]);
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! unlink([filename, ".msh"]);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%! elno_p1 = mesh.elements.tria6(elem_id_p1, :);

%! load_case.pressure.tria6.elements = elno_p1;
%! load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = E;
%! mesh.material_data.nu = nu;
%! mesh.material_data.rho = rho;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT, ...
%!                               FEM_VEC_LOAD_LUMPED], ...
%!                              load_case);
%!
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%! grp_id_displacement = find([[mesh.groups.tria6].id] == 4);
%! elem_id_displacement = mesh.groups.tria6(grp_id_displacement).elements;
%! elno_id_displacement = mesh.elements.tria6(elem_id_displacement, :);
%! delta = mean(sol_stat.def(elno_id_displacement, 3));
%! grp_id_stress = find([[mesh.groups.tria6].id] == 5);
%! elem_id_stress = mesh.groups.tria6(grp_id_stress).elements;
%! elno_id_stress = mesh.elements.tria6(elem_id_stress, :);
%! taum = zeros(6, numel(elno_id_stress));
%! taum_n = zeros(1, numel(elno_id_stress));
%! for i=1:numel(elno_id_stress)
%!   [ridx, cidx] = find(mesh.elements.tet10 == elno_id_stress(i));
%!   for j=1:numel(ridx)
%!     taum(:, i) += reshape(sol_stat.stress.taum.tet10(ridx(j), cidx(j), :), 6, 1);
%!     ++taum_n(i);
%!   endfor
%! endfor
%! taum *= diag(1 ./ taum_n);
%! sigma1_max = 0;
%! for i=1:columns(taum)
%!   TAU = [taum(1, i), taum(4, i), taum(6, i);
%!          taum(4, i), taum(2, i), taum(5, i);
%!          taum(6, i), taum(5, i), taum(3, i)];
%!   sigma1_max = max(sigma1_max, max(eig(TAU)));
%! endfor
%! fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%! fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%! fprintf(stderr, "delta=%.3f [mm]\n", delta);
%! ## K.J.Bathe page 329 4.20b
%! sigma1_max_ref = 0.6056;
%! delta_ref = -1.669;
%! fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%! fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%! if (animate)
%!   opt_anim.scale_def = 10;
%!   opt_anim.animation_delay = 1;
%!   opt_anim.print_and_exit = true;
%!   opt_anim.print_to_file = filename;
%!   opt_anim.rotation_angle = [90, 0, 0] * pi / 180;
%!   opt_anim.skin_only = true;
%!   opt_anim.show_element = true;
%!   unwind_protect
%!     fem_post_sol_external(mesh, sol_stat, opt_anim);
%!     [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title("Gmsh - deformed mesh / continuous stress tensor");
%!   unwind_protect_cleanup
%!     unlink([opt_anim.print_to_file, "_001.jpg"]);
%!   end_unwind_protect
%! endif
%! assert(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%! assert(delta, delta_ref, 0.04 * abs(delta_ref));
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!demo
%! ## DEMO2
%! close all;
%! a = 5e-3;
%! b = 5e-3;
%! c = 10e-3;
%! h = 0.5e-3;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! gamma = 10 * pi / 180;
%! tolstat = 4e-2;
%! scale = 5e-3;
%! maxdist = 1e-2 * max([a,b,c]);
%! eliminate = false;
%! animate = true;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! for i=1:4
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! unwind_protect
%! switch (i)
%! case {1,2}
%! order = 1;
%! optflags = {"-optimize_threshold", "0.95"};
%! otherwise
%! order = 2;
%! optflags = {"-optimize_ho", "-ho_min", "0.5", "-ho_max", "1.5"};
%! endswitch
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "a=%g;\n", a);
%! fprintf(fd, "b=%g;\n", b);
%! fprintf(fd, "c=%g;\n", c);
%! fprintf(fd, "h=%g;\n", h * order);
%! fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%! fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%! fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%! fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {0.0,0.0,c} {\n");
%! switch (i)
%! case 1
%! fprintf(fd, "  Surface{6}; Layers{%d}; Recombine;\n", ceil(c/h));
%! otherwise
%! fputs(fd, "  Surface{6};\n");
%! endswitch
%! fputs(fd, "};\n");
%! switch (i)
%! case 1
%! fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%! endswitch
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fprintf(fd, "Physical Surface(\"right\",%d) = {tmp[2]};\n", 100 * i + 1);
%! fprintf(fd, "Physical Surface(\"rear\",%d) = {tmp[3]};\n", 100 * i + 2);
%! fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[4]};\n", 100 * i + 3);
%! fprintf(fd, "Physical Surface(\"front\",%d) = {tmp[5]};\n", 100 * i + 4);
%! fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[0]};\n", 100 * i + 5);
%! fprintf(fd, "Physical Surface(\"bottom\",%d) = {6};\n", 100 * i + 6);
%! unwind_protect_cleanup
%! fclose(fd);
%! end_unwind_protect
%! fprintf(stderr, "meshing ...\n");
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", sprintf("%d", order), optflags{:}, [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%! data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(data(i).mesh.nodes));
%! unlink([filename, ".msh"]);
%! endfor
%! data(2).mesh.nodes(:, 1) += a;
%! data(3).mesh.nodes(:, 1) += a;
%! data(3).mesh.nodes(:, 2) += b;
%! data(4).mesh.nodes(:, 2) += b;
%! for i=1:numel(data)
%!   if (isfield(data(i).mesh.elements, "iso8"))
%!     data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!   else
%!     data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!   endif
%!   data(i).mesh.material_data.rho = rho;
%!   data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! endfor
%! mesh = fem_post_mesh_merge(data);
%! mesh.nodes(:, 3) -= 0.5 * c;
%! constr = FEM_CT_FIXED;
%! mesh.elements.sfncon4(1).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 103)).elements, :);
%! mesh.elements.sfncon4(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 401)).nodes(:);
%! mesh.elements.sfncon4(1).maxdist = maxdist;
%! mesh.elements.sfncon4(1).constraint = constr;
%! mesh.elements.sfncon4(2).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 102)).elements, :);
%! mesh.elements.sfncon4(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 204)).nodes(:);
%! mesh.elements.sfncon4(2).maxdist = maxdist;
%! mesh.elements.sfncon4(2).constraint = constr;
%! mesh.elements.sfncon6(1).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 402)).elements, :);
%! mesh.elements.sfncon6(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 304)).nodes(:);
%! mesh.elements.sfncon6(1).maxdist = maxdist;
%! mesh.elements.sfncon6(1).constraint = constr;
%! mesh.elements.sfncon6(2).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 301)).elements, :);
%! mesh.elements.sfncon6(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 203)).nodes(:);
%! mesh.elements.sfncon6(2).maxdist = maxdist;
%! mesh.elements.sfncon6(2).constraint = constr;
%! slave_nodes = [];
%! for i=1:numel(mesh.elements.sfncon4)
%!   idx_slave = [];
%!   for j=1:numel(slave_nodes)
%!     idx_slave = [idx_slave; find(mesh.elements.sfncon4(i).slave == slave_nodes(j))];
%!   endfor
%!   if (numel(idx_slave))
%!     mesh.elements.sfncon4(i).slave(idx_slave) = 0;
%!     mesh.elements.sfncon4(i).slave = mesh.elements.sfncon4(i).slave(find(mesh.elements.sfncon4(i).slave));
%!   endif
%!   slave_nodes = [slave_nodes; mesh.elements.sfncon4(i).slave];
%! endfor

%! for i=1:numel(mesh.elements.sfncon6)
%!   idx_slave = [];
%!   for j=1:numel(slave_nodes)
%!     idx_slave = [idx_slave; find(mesh.elements.sfncon6(i).slave == slave_nodes(j))];
%!   endfor
%!   if (numel(idx_slave))
%!     mesh.elements.sfncon6(i).slave(idx_slave) = 0;
%!     mesh.elements.sfncon6(i).slave = mesh.elements.sfncon6(i).slave(find(mesh.elements.sfncon6(i).slave));
%!   endif
%!   slave_nodes = [slave_nodes; mesh.elements.sfncon6(i).slave];
%! endfor
%! group_id4 = [[mesh.groups.iso4].id];
%! group_id6 = [[mesh.groups.tria6].id];

%! node_front = [[mesh.groups.iso4(find(group_id4 == 104))].nodes, [mesh.groups.tria6(find(group_id6 == 404))].nodes];
%! node_rear = [[mesh.groups.iso4(find(group_id4 == 202))].nodes, [mesh.groups.tria6(find(group_id6 == 302))].nodes];
%! node_right = [mesh.groups.iso4(find((group_id4 == 101) | (group_id4 == 201))).nodes, ...
%!               mesh.groups.tria6(find((group_id6==303)|(group_id6==403))).nodes];
%! node_bottom = [[mesh.groups.iso4(find(mod(group_id4, 100) == 6))].nodes, ...
%!                [mesh.groups.tria6(find(mod(group_id6,100) == 6))].nodes];
%! node_top = [[mesh.groups.iso4(find(mod(group_id4, 100) == 5))].nodes, ...
%!             [mesh.groups.tria6(find(mod(group_id6,100) == 5))].nodes];
%! node_constr = [node_front, node_rear, node_right, node_bottom, node_top];
%! mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%! load_case.joints = repmat(struct("U",[]), 1, numel(node_constr));
%! idx_joint = int32(0);
%! for i=1:numel(node_front)
%!   if (~numel(find(slave_nodes == node_front(i))))
%!   mesh.elements.joints(++idx_joint).C = [[1, 0, 0; 0, 0, 1],  zeros(2, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_front(i);
%!   load_case.joints(idx_joint).U = [gamma * mesh.nodes(node_front(i), 3); 0];
%!   endif
%! endfor
%! for i=1:numel(node_rear)
%!   if (~numel(find(slave_nodes == node_rear(i))))
%!     mesh.elements.joints(++idx_joint).C = [[1, 0, 0; 0, 0, 1], zeros(2, 3)];
%!     mesh.elements.joints(idx_joint).nodes = node_rear(i);
%!     load_case.joints(idx_joint).U = [gamma * mesh.nodes(node_rear(i), 3); 0];
%!   endif
%! endfor

%! for i=1:numel(node_right)
%!   if (~numel(find(slave_nodes == node_right(i))))
%!   mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!   mesh.elements.joints(idx_joint).nodes = node_right(i);
%!   load_case.joints(idx_joint).U = 0;
%!   endif
%! endfor

%! for i=1:numel(node_bottom)
%!   xi = mesh.nodes(node_bottom(i), 1);
%!   if (~(numel(find(slave_nodes == node_bottom(i))) || xi == 0 || xi == 2 * a))
%!     mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!     mesh.elements.joints(idx_joint).nodes = node_bottom(i);
%!     load_case.joints(idx_joint).U = 0;
%!   endif
%! endfor

%! for i=1:numel(node_top)
%!   xi = mesh.nodes(node_top(i), 1);
%!   if (~(numel(find(slave_nodes == node_top(i))) || xi == 0 || xi == 2 * a))
%!     mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!     mesh.elements.joints(idx_joint).nodes = node_top(i);
%!     load_case.joints(idx_joint).U = 0;
%!   endif
%! endfor
%!
%! mesh.elements.joints = mesh.elements.joints(1:idx_joint);
%! load_case.joints = load_case.joints(1:idx_joint);
%! load_case.locked_dof = false(size(mesh.nodes));
%! fprintf(stderr, "assembling matrices ...\n");
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! if (eliminate)
%!   K_mat_type = FEM_MAT_STIFFNESS;
%! else
%!   K_mat_type = FEM_MAT_STIFFNESS;
%! endif
%! [mat_ass.K, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [K_mat_type, ...
%!                                      FEM_VEC_LOAD_CONSISTENT], ...
%!                                     load_case);
%! if (eliminate)
%! fprintf(stderr, "eliminating constraint equations ...\n");
%! [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%! fprintf(stderr, "linear static solution ...\n");
%! opt_ls.refine_max_iter = int32(100);
%! if (fem_sol_check_func("mumps"))
%! opt_ls.matrix_type = MUMPS_MAT_SYM;
%! opt_ls.verbose = MUMPS_VER_ERR;
%! Kfact = fem_fact_mumps(Kred, opt_ls);
%! else
%! Kfact = fem_fact_umfpack(Kred, opt_ls);
%! endif
%! Ured = Kfact \ Rred;
%! fprintf(stderr, "expanding static solution ...\n");
%! sol_stat.def = fem_post_def_nodal(mesh, dof_map, Tred * Ured);
%! else
%! fprintf(stderr, "linear static solution ...\n");
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! endif

%! [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%! if (animate)
%!   opt_anim.scale_def = 1;
%!   opt_anim.animation_delay = 1;
%!   opt_anim.print_and_exit = true;
%!   opt_anim.print_to_file = filename;
%!   opt_anim.rotation_angle = [286, 2, 205] * pi / 180;
%!   opt_anim.skin_only = true;
%!   opt_anim.show_element = true;
%!   unwind_protect
%!     fem_post_sol_external(mesh, sol_stat, opt_anim);
%!     [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title("Gmsh - deformed mesh / continuous stress tensor");
%!   unwind_protect_cleanup
%!     unlink([opt_anim.print_to_file, "_001.jpg"]);
%!   end_unwind_protect
%! endif
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
