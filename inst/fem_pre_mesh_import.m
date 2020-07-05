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

function [mesh, load_case] = fem_pre_mesh_import(filename, format, options)
  if (nargin < 1 || nargin > 3 || nargout > 2)
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

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "promote_elem"))
    options.promote_elem = {"tet4", "tria3", "prism6"};
  endif

  switch (format)
    case "gmsh"
      mesh = fem_load_mesh_gmsh(filename, format, options);
    case "apdl"
      mesh = fem_load_mesh_apdl(filename, format);
    case "eossp"
      [mesh, load_case] = fem_load_mesh_eossp(filename, format);
    otherwise
      error("mesh format \"%s\" not supported", format);
  endswitch
endfunction

function mesh = fem_load_mesh_gmsh(filename, format, options)
  mesh.nodes = zeros(0, 6);
  mesh.elements = struct();
  mesh.groups = struct();

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(filename, "rt");

    if (fd == -1)
      error("failed to open file \"%s\": %s", filename, msg);
    endif

    nodes = zeros(0, 4);
    elements = zeros(0, 15, "int32");
    p_name = {};
    p_dim = [];
    p_id = [];

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

          eltype = struct("name", {"tet10", "tria6", "tet4", "prism6", "tria3", "iso8", "iso4", "iso20", "quad8"}, ...
                          "promote", {-1, -1, 6, 6, 7, -1, -1, -1, -1}, ...
                          "id", {11, 9, 4, 6, 2, 5, 3, 17, 16}, ...
                          "dim", {3, 2, 3, 3, 2, 3, 2, 3, 2}, ...
                          "norder", {[1:8, 10, 9], ...
                                     1:6, ...
                                     [4, 4, 4, 4, 1:3, 3], ...
                                     [6,4,4,5,3,1,1,2], ...
                                     [1:3, 3], ...
                                     [5:8, 1:4], ...
				     1:4, ...
				     [5:8, 1:4, 17, 19, 20, 18, 9, 12, 14, 10, 11, 13, 15, 16], ...
				     1:8}, ...
                         "nordernonp", {[],[],[],[],[1:3],[],[],[],[]});

          for i=1:numel(eltype)
            if (eltype(i).promote > 0)
              switch (eltype(i).name)
                case options.promote_elem
                otherwise
                  if (~isempty(eltype(i).nordernonp))
                    eltype(i).promote = -1;
                    eltype(i).norder = eltype(i).nordernonp;
                  endif
              endswitch
            endif
          endfor

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
                                                   numel(eltype(k).norder) * numel(mshgrp(end).elements)));

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

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(filename, "rt");

    if (fd == -1)
      error("failed to open file \"%s\": %s", filename, msg);
    endif

    nodes = zeros(0, 4);
    elements = zeros(0, 15, "int32");
    elemtype = -1;

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

  fd = -1;

  unwind_protect
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
    if (fd ~= -1)
      fclose(fd);
    endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 30e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 3.5e-3;
%!   p = 25e6;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   if (do_rotate)
%!     fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!   endif
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
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
%!     group_defs(1).elem_type = "tria6";
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
%!     group_defs(2).elem_type = "tria6";
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
%!     group_defs(3).elem_type = "tria6";
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
%!     group_defs(4).elem_type = "tria6";
%!     groups = fem_pre_mesh_groups_create(mesh, group_defs, sqrt(eps));
%!     assert(numel(groups.tria6), 4);
%!     assert([groups.tria6.id], [group_defs.id]);
%!     assert(groups.tria6(1).nodes, mesh.groups.tria6(1).nodes);
%!     assert(groups.tria6(2).nodes, mesh.groups.tria6(1).nodes);
%!     assert(groups.tria6(3).nodes, mesh.groups.tria6(2).nodes);
%!     assert(groups.tria6(4).nodes, mesh.groups.tria6(2).nodes);
%!     assert(groups.tria6(1).elements, mesh.groups.tria6(1).elements);
%!     assert(groups.tria6(2).elements, mesh.groups.tria6(1).elements);
%!     assert(groups.tria6(3).elements, mesh.groups.tria6(2).elements);
%!     assert(groups.tria6(4).elements, mesh.groups.tria6(2).elements);
%!   endif
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   load_case.pressure.tria6.elements = mesh.elements.tria6(mesh.groups.tria6(find([mesh.groups.tria6.id] == 2)).elements, :);
%!   Xp = mesh.nodes(load_case.pressure.tria6.elements, 1:3);
%!   xp = reshape(Xp(:, 1), rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%!   yp = reshape(Xp(:, 2), rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%!   zp = reshape(Xp(:, 3), rows(load_case.pressure.tria6.elements), columns(load_case.pressure.tria6.elements));
%!   load_case.pressure.tria6.p = p / 2 * (yp / b + zp / c);
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.Mlumped, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_MAT_MASS_LUMPED, ...
%!                            FEM_VEC_LOAD_CONSISTENT, ...
%!                            FEM_VEC_LOAD_LUMPED, ...
%!                            FEM_SCA_TOT_MASS], load_case);
%!   assert(mtot, a * b * c * mesh.material_data.rho, sqrt(eps) * a * b * c * mesh.material_data.rho);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   [sol_eig_lumped] = fem_sol_modal(mesh, dof_map, setfield(mat_ass, "M", mat_ass.Mlumped), number_of_modes);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%!   X = mesh.nodes(unique(load_case.pressure.tria6.elements), 1:3).';
%!   dof_idx = dof_map.ndof(unique(load_case.pressure.tria6.elements), 1:3);
%!   F_con = full(mat_ass.R(dof_idx)).';
%!   F_lumped = full(mat_ass.Rlumped(dof_idx)).';
%!   M_con = cross(X, F_con);
%!   M_lumped = cross(X, F_lumped);
%!   Ftot_con = sum(F_con, 2);
%!   Mtot_con = sum(M_con, 2);
%!   Ftot_lumped = sum(F_lumped, 2);
%!   Mtot_lumped = sum(M_lumped, 2);
%!   Fx_an = -(b * c * p) / 2; ## grind(integrate(integrate(-1/2*p*(y/b+z/c),z,0,c),y,0,b));
%!   Mz_an = (7 * b^2 * c * p) / 24; ## grind(integrate(integrate(1/2*p*(y/b+z/c)*y,z,0,c),y,0,b));
%!   My_an = -(7 * b * c^2 * p) / 24; ## grind(integrate(integrate(-1/2*p*(y/b+z/c)*z,z,0,c),y,0,b));
%!   F_an = [Fx_an; 0; 0];
%!   M_an = [0; My_an; Mz_an];
%!   assert(Ftot_con, F_an, eps^0.9 * norm(F_an));
%!   assert(Ftot_lumped, F_an, eps^0.9 * norm(F_an));
%!   assert(Mtot_con, M_an, eps^0.9 * norm(M_an));
%!   assert(Mtot_lumped, M_an, 5e-3 * norm(M_an));
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
%!   assert(all(f_lumped <= f));
%!   assert(f, f_ref, tol * max(f_ref));

%!   figure("visible","off");
%!   hold on;
%!   fem_post_sol_plot(mesh);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title('undeformed mesh');
%!   opts_plot.elem_types = {"tria6", "tet10"};
%!   opts_plot.elem_groups.tria6 = [mesh.groups.tria6.id];
%!   opts_plot.elem_groups.tet10 = [mesh.groups.tet10.id];
%!   for i=1:min(number_of_modes, length(sol_eig.f))
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, :, i), "rows")), i, opts_plot);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("%d. eigenmode: %gHz",i,sol_eig.f(i)));
%!   endfor

%!   figure("visible", "off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat, scale_eig/max(norm(sol_stat.def, "rows")), 1, opts_plot);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title("deformed mesh consistent load vector");

%!   figure("visible", "off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat_lumped, scale_eig/max(norm(sol_stat_lumped.def, "rows")), 1, opts_plot);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title("deformed mesh lumped load vector");

%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST2
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   ri = 8e-3;
%!   ro = 10e-3;
%!   h = 12e-3;
%!   c = 2e-3;
%!   b = h - 2 * c;
%!   p1 = 25.79e6;
%!   p2 = 7.83e6;
%!   p3 = 1.3758e6;
%!   scale_def = 5e-3;
%!   mesh_size = 3e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "ri = %g;\n", ri);
%!   fprintf(fd, "ro = %g;\n", ro);
%!   fprintf(fd, "h = %g;\n", h);
%!   fprintf(fd, "c = %g;\n", c);
%!   fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!   fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!   fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!   fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!   fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!   fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!   fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!   fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Line(7) = {7,8};\n");
%!   fputs(fd, "Line(8) = {8,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!   fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%!   grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!   grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%!   grp_id_p3 = find([[mesh.groups.tria6].id] == 4);
%!   elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%!   elem_id_p2 = mesh.groups.tria6(grp_id_p2).elements;
%!   elem_id_p3 = mesh.groups.tria6(grp_id_p3).elements;
%!   elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%!   elno_p2 = mesh.elements.tria6(elem_id_p2, :);
%!   elno_p3 = mesh.elements.tria6(elem_id_p3, :);

%!   load_case.pressure.tria6.elements = [elno_p1; elno_p2; elno_p3];
%!   load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1));
%!                                 repmat(p2, rows(elno_p2), columns(elno_p2));
%!                                 repmat(p3, rows(elno_p3), columns(elno_p3))];

%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT, ...
%!                                       FEM_VEC_LOAD_LUMPED], ...
%!                                      load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat_lumped, scale_def / max(norm(sol_stat_lumped.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - lumped pressure load");

%!   X1 = mesh.nodes(unique(elno_p1), 1:3).';
%!   X2 = mesh.nodes(unique(elno_p2), 1:3).';
%!   X3 = mesh.nodes(unique(elno_p3), 1:3).';
%!   dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%!   dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%!   dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%!   F1_con = full(mat_ass.R(dof1)).';
%!   F2_con = full(mat_ass.R(dof2)).';
%!   F3_con = full(mat_ass.R(dof3)).';
%!   M1_con = cross(X1, F1_con);
%!   M2_con = cross(X2, F2_con);
%!   M3_con = cross(X3, F3_con);
%!   F1_lumped = full(mat_ass.Rlumped(dof1)).';
%!   F2_lumped = full(mat_ass.Rlumped(dof2)).';
%!   F3_lumped = full(mat_ass.Rlumped(dof3)).';
%!   M1_lumped = cross(X1, F1_lumped);
%!   M2_lumped = cross(X2, F2_lumped);
%!   M3_lumped = cross(X3, F3_lumped);
%!   Ftot1_con = sum(F1_con, 2);
%!   Ftot2_con = sum(F2_con, 2);
%!   Ftot3_con = sum(F3_con, 2);
%!   Mtot1_con = sum(M1_con, 2);
%!   Mtot2_con = sum(M2_con, 2);
%!   Mtot3_con = sum(M3_con, 2);

%!   Ftot1_lumped = sum(F1_lumped, 2);
%!   Ftot2_lumped = sum(F2_lumped, 2);
%!   Ftot3_lumped = sum(F3_lumped, 2);
%!   Mtot1_lumped = sum(M1_lumped, 2);
%!   Mtot2_lumped = sum(M2_lumped, 2);
%!   Mtot3_lumped = sum(M3_lumped, 2);

%!   ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%!   F1_an = [ri * b * p1;
%!            ri * b * p1;
%!            0];

%!   M1_an = [-ri * b * p1 * (c + b/2);
%!            ri * b * p1 * (c + b/2);
%!            0];

%!   F2_an = [-ro * b * p2;
%!            -ro * b * p2;
%!            0];

%!   M2_an = [ ro * b * p2 * (c + b/2);
%!             -ro * b * p2 * (c + b/2);
%!             0];

%!   F3_an = [0;
%!            0;
%!            -p3 * (ro^2 - ri^2) * pi / 4];

%!   M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            0];

%!   assert(Ftot1_con, F1_an, eps^0.9 * norm(F1_an));
%!   assert(Ftot2_con, F2_an, eps^0.9 * norm(F2_an));
%!   assert(Ftot1_lumped, F1_an, eps^0.9 * norm(F1_an));
%!   assert(Ftot2_lumped, F2_an, eps^0.9 * norm(F2_an));

%!   assert(Mtot1_con, M1_an, eps^0.9 * norm(M1_an));
%!   assert(Mtot2_con, M2_an, eps^0.9 * norm(M2_an));
%!   assert(Mtot1_lumped, M1_an, eps^0.2 * norm(M1_an));
%!   assert(Mtot2_lumped, M2_an, eps^0.2 * norm(M2_an));

%!   assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%!   assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%!   assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%!   assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif

%!   a = 30e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   d = -5e-3;
%!   e = 35e-3;
%!   h = 4e-3;
%!   alpha = 1e-6;
%!   beta = 1e-4;
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"modal\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"interfaces\",2) = {tmp[2]};\n");
%!   fputs(fd, "SetFactory(\"Built-in\");\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   cms_opt.modes.number = int32(6);
%!   cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!   cms_opt.nodes.interfaces.number = int32(rows(mesh.nodes) + 2);
%!   cms_opt.invariants = false;
%!   cms_opt.algorithm = "unsymmetric";
%!   mesh.nodes(cms_opt.nodes.modal.number, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%!   mesh.nodes([cms_opt.nodes.interfaces.number], :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                    [1, 2], ...
%!                                                    [cms_opt.nodes.modal.number, ...
%!                                                     cms_opt.nodes.interfaces.number]);
%!   for i=1:numel(mesh.elements.rbe3)
%!     assert(sum(mesh.elements.rbe3(i).weight), b * c, sqrt(eps) * (b * c));
%!   endfor
%!   [mesh, mat_ass_cms, dof_map_cms, sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%!   mat_ass_cms.Dred = alpha * mat_ass_cms.Mred + beta * mat_ass_cms.Kred;
%!   
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
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST4
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   ri = 8e-3;
%!   ro = 10e-3;
%!   h = 12e-3;
%!   c = 2e-3;
%!   b = h - 2 * c;
%!   p1 = 25.79e6;
%!   p2 = 7.83e6;
%!   p3 = 1.3758e6;
%!   scale_def = 5e-3;
%!   mesh_size = 1e-3;
%!   enable_linear_dist = false;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "ri = %g;\n", ri);
%!   fprintf(fd, "ro = %g;\n", ro);
%!   fprintf(fd, "h = %g;\n", h);
%!   fprintf(fd, "c = %g;\n", c);
%!   fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!   fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!   fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!   fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!   fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!   fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!   fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!   fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Line(7) = {7,8};\n");
%!   fputs(fd, "Line(8) = {8,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!   fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%!   grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!   grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%!   grp_id_p3 = find([[mesh.groups.tria6].id] == 4);
%!   elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%!   elem_id_p2 = mesh.groups.tria6(grp_id_p2).elements;
%!   elem_id_p3 = mesh.groups.tria6(grp_id_p3).elements;
%!   elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%!   elno_p2 = mesh.elements.tria6(elem_id_p2, :);
%!   elno_p3 = mesh.elements.tria6(elem_id_p3, :);
%!   x1 = mesh.nodes(:, 1)(elno_p1);
%!   y1 = mesh.nodes(:, 2)(elno_p1);
%!   z1 = mesh.nodes(:, 3)(elno_p1);
%!   Phi1 = atan2(y1, x1);

%!   x2 = mesh.nodes(:, 1)(elno_p2);
%!   y2 = mesh.nodes(:, 2)(elno_p2);
%!   z2 = mesh.nodes(:, 3)(elno_p2);
%!   Phi2 = atan2(y2, x2);

%!   load_case.pressure.tria6.elements = [elno_p1; elno_p2; elno_p3];
%!   load_case.pressure.tria6.p = [p1 * Phi1 / (pi / 2) .* z1 / h;
%!                                 p2 * Phi2 / (pi / 2) .* z2 / h;
%!                                 repmat(p3, rows(elno_p3), columns(elno_p3))];
%!   if (enable_linear_dist)
%!     p_mid = load_case.pressure.tria6.p(:, 1:3);
%!     load_case.pressure.tria6.p(:, 4:6) = [0.5 * (p_mid(:, 1) + p_mid(:, 2)), 0.5 * (p_mid(:, 2) + p_mid(:, 3)), 0.5 * (p_mid(:, 1) + p_mid(:, 3))];
%!   endif
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT, ...
%!                                       FEM_VEC_LOAD_LUMPED], ...
%!                                      load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat_lumped, scale_def / max(norm(sol_stat_lumped.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - lumped pressure load");

%!   X1 = mesh.nodes(unique(elno_p1), 1:3).';
%!   X2 = mesh.nodes(unique(elno_p2), 1:3).';
%!   X3 = mesh.nodes(unique(elno_p3), 1:3).';
%!   dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%!   dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%!   dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%!   F1_con = full(mat_ass.R(dof1)).';
%!   F2_con = full(mat_ass.R(dof2)).';
%!   F3_con = full(mat_ass.R(dof3)).';
%!   M1_con = cross(X1, F1_con);
%!   M2_con = cross(X2, F2_con);
%!   M3_con = cross(X3, F3_con);
%!   F1_lumped = full(mat_ass.Rlumped(dof1)).';
%!   F2_lumped = full(mat_ass.Rlumped(dof2)).';
%!   F3_lumped = full(mat_ass.Rlumped(dof3)).';
%!   M1_lumped = cross(X1, F1_lumped);
%!   M2_lumped = cross(X2, F2_lumped);
%!   M3_lumped = cross(X3, F3_lumped);
%!   Ftot1_con = sum(F1_con, 2);
%!   Ftot2_con = sum(F2_con, 2);
%!   Ftot3_con = sum(F3_con, 2);
%!   Mtot1_con = sum(M1_con, 2);
%!   Mtot2_con = sum(M2_con, 2);
%!   Mtot3_con = sum(M3_con, 2);

%!   Ftot1_lumped = sum(F1_lumped, 2);
%!   Ftot2_lumped = sum(F2_lumped, 2);
%!   Ftot3_lumped = sum(F3_lumped, 2);
%!   Mtot1_lumped = sum(M1_lumped, 2);
%!   Mtot2_lumped = sum(M2_lumped, 2);
%!   Mtot3_lumped = sum(M3_lumped, 2);

%!   F1_an = [(2*(pi/2-1)*((h^2-2*c*h+c^2)/2-c^2/2)*p1*ri)/(pi*h);
%!            (2*((h^2-2*c*h+c^2)/2-c^2/2)*p1*ri)/(pi*h);
%!            0];

%!   F2_an = [-(2*(pi/2-1)*((h^2-2*c*h+c^2)/2-c^2/2)*p2*ro)/(pi*h);
%!            -(2*((h^2-2*c*h+c^2)/2-c^2/2)*p2*ro)/(pi*h);
%!            0];

%!   M1_an = [-(2*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p1*ri)/(pi*h);
%!            (2*(pi/2-1)*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p1*ri)/(pi*h);
%!            0];

%!   M2_an = [(2*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p2*ro)/(pi*h);
%!            -(2*(pi/2-1)*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p2*ro)/(pi*h);
%!            0];

%!   ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%!   F3_an = [0;
%!            0;
%!            -p3 * (ro^2 - ri^2) * pi / 4];

%!   M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            0];

%!   assert(Ftot1_con, F1_an, 1e-4 * norm(F1_an));
%!   assert(Ftot2_con, F2_an, 1e-4 * norm(F2_an));
%!   assert(Ftot1_lumped, F1_an, 2e-3 * norm(F1_an));
%!   assert(Ftot2_lumped, F2_an, 2e-3 * norm(F2_an));

%!   assert(Mtot1_con, M1_an, 1e-4 * norm(M1_an));
%!   assert(Mtot2_con, M2_an, 1e-4 * norm(M2_an));
%!   assert(Mtot1_lumped, M1_an, 5e-3 * norm(M1_an));
%!   assert(Mtot2_lumped, M2_an, 5e-3 * norm(M2_an));

%!   assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%!   assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%!   assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%!   assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));

%!   fprintf(stderr, "err(F1_con)=%e\n", norm(Ftot1_con - F1_an) / norm(F1_an));
%!   fprintf(stderr, "err(F2_con)=%e\n", norm(Ftot2_con - F2_an) / norm(F2_an));
%!   fprintf(stderr, "err(F3_con)=%e\n", norm(Ftot3_con - F3_an) / norm(F3_an));
%!   fprintf(stderr, "err(M1_con)=%e\n", norm(Mtot1_con - M1_an) / norm(M1_an));
%!   fprintf(stderr, "err(M2_con)=%e\n", norm(Mtot2_con - M2_an) / norm(M2_an));
%!   fprintf(stderr, "err(M3_con)=%e\n", norm(Mtot3_con - M3_an) / norm(M3_an));
%!   fprintf(stderr, "err(F1_lumped)=%e\n", norm(Ftot1_lumped - F1_an) / norm(F1_an));
%!   fprintf(stderr, "err(F2_lumped)=%e\n", norm(Ftot2_lumped - F2_an) / norm(F2_an));
%!   fprintf(stderr, "err(F3_lumped)=%e\n", norm(Ftot3_lumped - F3_an) / norm(F3_an));
%!   fprintf(stderr, "err(M1_lumped)=%e\n", norm(Mtot1_lumped - M1_an) / norm(M1_an));
%!   fprintf(stderr, "err(M2_lumped)=%e\n", norm(Mtot2_lumped - M2_an) / norm(M2_an));
%!   fprintf(stderr, "err(M3_lumped)=%e\n", norm(Mtot3_lumped - M3_an) / norm(M3_an));

%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST5
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   ri = 8e-3;
%!   ro = 10e-3;
%!   h = 12e-3;
%!   c = 2e-3;
%!   b = h - 2 * c;
%!   p1 = 25.79e6;
%!   p2 = 7.83e6;
%!   p3 = 1.3758e6;
%!   scale_def = 100;
%!   mesh_size = 2e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "ri = %g;\n", ri);
%!   fprintf(fd, "ro = %g;\n", ro);
%!   fprintf(fd, "h = %g;\n", h);
%!   fprintf(fd, "c = %g;\n", c);
%!   fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!   fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!   fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!   fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!   fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!   fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!   fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!   fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Line(7) = {7,8};\n");
%!   fputs(fd, "Line(8) = {8,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!   fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%!   grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!   grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%!   grp_id_p3 = find([[mesh.groups.tria6].id] == 4);
%!   elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%!   elem_id_p2 = mesh.groups.tria6(grp_id_p2).elements;
%!   elem_id_p3 = mesh.groups.tria6(grp_id_p3).elements;
%!   elno_p1 = mesh.elements.tria6(elem_id_p1, :);
%!   elno_p2 = mesh.elements.tria6(elem_id_p2, :);
%!   elno_p3 = mesh.elements.tria6(elem_id_p3, :);
%!   x1 = mesh.nodes(:, 1)(elno_p1);
%!   y1 = mesh.nodes(:, 2)(elno_p1);
%!   z1 = mesh.nodes(:, 3)(elno_p1);
%!   Phi1 = atan2(y1, x1);

%!   x2 = mesh.nodes(:, 1)(elno_p2);
%!   y2 = mesh.nodes(:, 2)(elno_p2);
%!   z2 = mesh.nodes(:, 3)(elno_p2);
%!   Phi2 = atan2(y2, x2);

%!   load_case.pressure.tria6.elements = [elno_p1; elno_p2; elno_p3];
%!   load_case.pressure.tria6.p = [p1 * cos(Phi1) .* sin(pi * z1 / h);
%!                                 p2 * cos(Phi2) .* sin(pi * z2 / h);
%!                                 repmat(p3, rows(elno_p3), columns(elno_p3))];

%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT, ...
%!                                       FEM_VEC_LOAD_LUMPED], ...
%!                                      load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat_lumped, scale_def);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - lumped pressure load");

%!   X1 = mesh.nodes(unique(elno_p1), 1:3).';
%!   X2 = mesh.nodes(unique(elno_p2), 1:3).';
%!   X3 = mesh.nodes(unique(elno_p3), 1:3).';
%!   dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%!   dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%!   dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%!   F1_con = full(mat_ass.R(dof1)).';
%!   F2_con = full(mat_ass.R(dof2)).';
%!   F3_con = full(mat_ass.R(dof3)).';
%!   M1_con = cross(X1, F1_con);
%!   M2_con = cross(X2, F2_con);
%!   M3_con = cross(X3, F3_con);
%!   F1_lumped = full(mat_ass.Rlumped(dof1)).';
%!   F2_lumped = full(mat_ass.Rlumped(dof2)).';
%!   F3_lumped = full(mat_ass.Rlumped(dof3)).';
%!   M1_lumped = cross(X1, F1_lumped);
%!   M2_lumped = cross(X2, F2_lumped);
%!   M3_lumped = cross(X3, F3_lumped);
%!   Ftot1_con = sum(F1_con, 2);
%!   Ftot2_con = sum(F2_con, 2);
%!   Ftot3_con = sum(F3_con, 2);
%!   Mtot1_con = sum(M1_con, 2);
%!   Mtot2_con = sum(M2_con, 2);
%!   Mtot3_con = sum(M3_con, 2);

%!   Ftot1_lumped = sum(F1_lumped, 2);
%!   Ftot2_lumped = sum(F2_lumped, 2);
%!   Ftot3_lumped = sum(F3_lumped, 2);
%!   Mtot1_lumped = sum(M1_lumped, 2);
%!   Mtot2_lumped = sum(M2_lumped, 2);
%!   Mtot3_lumped = sum(M3_lumped, 2);

%!   F1_an = [(pi*((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p1*ri)/4;
%!            (((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p1*ri)/2;
%!            0];
%!   F2_an = [-(pi*((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p2*ro)/4;
%!            -(((cos((pi*c)/h)*h)/pi-(h*cos((pi*h-pi*c)/h))/pi)*p2*ro)/2;
%!            0];
%!   M1_an = [-(((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p1*ri)/2;
%!            (pi*((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p1*ri)/4;
%!            0];
%!   M2_an = [(((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p2*ro)/2;
%!            -(pi*((h^2*sin((pi*h-pi*c)/h)+(pi*c*h-pi*h^2)*cos((pi*h-pi*c)/h))/pi^2-(sin((pi*c)/h)*h^2-pi*c*cos((pi*c)/h)*h)/pi^2)*p2*ro)/4;
%!            0];
%!   ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%!   F3_an = [0;
%!            0;
%!            -p3 * (ro^2 - ri^2) * pi / 4];

%!   M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            0];

%!   fprintf(stderr, "err(F1_con)=%e\n", norm(Ftot1_con - F1_an) / norm(F1_an));
%!   fprintf(stderr, "err(F2_con)=%e\n", norm(Ftot2_con - F2_an) / norm(F2_an));
%!   fprintf(stderr, "err(F3_con)=%e\n", norm(Ftot3_con - F3_an) / norm(F3_an));
%!   fprintf(stderr, "err(M1_con)=%e\n", norm(Mtot1_con - M1_an) / norm(M1_an));
%!   fprintf(stderr, "err(M2_con)=%e\n", norm(Mtot2_con - M2_an) / norm(M2_an));
%!   fprintf(stderr, "err(M3_con)=%e\n", norm(Mtot3_con - M3_an) / norm(M3_an));
%!   fprintf(stderr, "err(F1_lumped)=%e\n", norm(Ftot1_lumped - F1_an) / norm(F1_an));
%!   fprintf(stderr, "err(F2_lumped)=%e\n", norm(Ftot2_lumped - F2_an) / norm(F2_an));
%!   fprintf(stderr, "err(F3_lumped)=%e\n", norm(Ftot3_lumped - F3_an) / norm(F3_an));
%!   fprintf(stderr, "err(M1_lumped)=%e\n", norm(Mtot1_lumped - M1_an) / norm(M1_an));
%!   fprintf(stderr, "err(M2_lumped)=%e\n", norm(Mtot2_lumped - M2_an) / norm(M2_an));
%!   fprintf(stderr, "err(M3_lumped)=%e\n", norm(Mtot3_lumped - M3_an) / norm(M3_an));

%!   assert(Ftot1_con, F1_an, eps^0.3 * norm(F1_an));
%!   assert(Ftot2_con, F2_an, eps^0.3 * norm(F2_an));
%!   assert(Ftot1_lumped, F1_an, 2e-2 * norm(F1_an));
%!   assert(Ftot2_lumped, F2_an, 2e-2 * norm(F2_an));

%!   assert(Mtot1_con, M1_an, eps^0.3 * norm(M1_an));
%!   assert(Mtot2_con, M2_an, eps^0.3 * norm(M2_an));
%!   assert(Mtot1_lumped, M1_an, 1e-2 * norm(M1_an));
%!   assert(Mtot2_lumped, M2_an, 1e-2 * norm(M2_an));

%!   assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%!   assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%!   assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%!   assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));

%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 15e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 3.5e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   if (do_rotate)
%!     fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!   endif
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!   mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   fprintf(stderr, "%d nodes\n", rows(mesh1.nodes));
%!   unlink([filename, ".msh"]);
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
%!   endif
%!   mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%!   mesh1.material_data.E = 210000e6;
%!   mesh1.material_data.nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%!   mesh2 = mesh1;
%!   mesh2.nodes(:, 1) += a;
%!   for i=1:numel(mesh2.groups.tria6)
%!     mesh2.groups.tria6(i).id += 100;
%!   endfor
%!   data(1).mesh = mesh1;
%!   data(2).mesh = mesh2;
%!   [mesh] = fem_post_mesh_merge(data);
%!   mesh.elements.sfncon6.slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 2)).nodes(:);
%!   mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 101)).elements, :);
%!   mesh.elements.sfncon6.maxdist = sqrt(eps) * max(abs([a,b,c]));
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], load_case);
%!   assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_SCA_STRESS_VMIS], ...
%!                                     load_case, ...
%!                                     sol_eig);

%!   for i=1:numel(sol_eig.f)
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!   endfor

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
%!   for i=1:length(sol_eig)
%!     fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig(i).f);
%!   endfor
%!   assert(sol_eig.f(:), f_ref, tol * max(f_ref));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 10e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 3.5e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   if (do_rotate)
%!     fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!   endif
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
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
%!   endif
%!   mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%!   mesh1.material_data.E = 210000e6;
%!   mesh1.material_data.nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%!   data.mesh = mesh1;
%!   data = repmat(data, 1, 3);
%!   for i=2:3
%!     data(i).mesh.nodes(:, 1) += a * (i - 1);
%!     for j=1:numel(data(i).mesh.groups.tria6)
%!       data(i).mesh.groups.tria6(j).id += 100 * (i - 1);
%!     endfor
%!   endfor
%!   [mesh] = fem_post_mesh_merge(data);

%!   for i=1:2
%!     grp_idx_slave = find([[mesh.groups.tria6].id] == (i - 1) * 100 + 2);
%!     grp_idx_master = mesh.groups.tria6(find([[mesh.groups.tria6].id] == i * 100 + 1)).elements;
%!     mesh.elements.sfncon6(i).slave = mesh.groups.tria6(grp_idx_slave).nodes(:);
%!     mesh.elements.sfncon6(i).master = mesh.elements.tria6(grp_idx_master, :);
%!     mesh.elements.sfncon6(i).maxdist = sqrt(eps) * max(abs([a,b,c]));
%!   endfor
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], ...
%!                           load_case);
%!   assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   for i=1:numel(sol_eig.f)
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3,i), "rows")), i);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!   endfor
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
%!   for i=1:length(sol_eig.f)
%!     fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig.f(i));
%!   endfor
%!   assert(sol_eig.f(:), f_ref, tol * max(f_ref));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 10e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 3.5e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   if (do_rotate)
%!     fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!   endif
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
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
%!   endif
%!   mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%!   mesh1.material_data.E = 210000e6;
%!   mesh1.material_data.nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%!   data.mesh = mesh1;
%!   data = repmat(data, 1, 3);
%!   for i=2:3
%!     data(i).mesh.nodes(:, 1) += a * (i - 1);
%!     for j=1:numel(data(i).mesh.groups.tria6)
%!       data(i).mesh.groups.tria6(j).id += 100 * (i - 1);
%!     endfor
%!   endfor
%!   [mesh] = fem_post_mesh_merge(data);

%!   for i=1:2
%!     grp_idx_slave = find([[mesh.groups.tria6].id] == i * 100 + 1);
%!     grp_idx_master = mesh.groups.tria6(find([[mesh.groups.tria6].id] == (i - 1) * 100 + 2)).elements;
%!     mesh.elements.sfncon6(i).slave = mesh.groups.tria6(grp_idx_slave).nodes(:);
%!     mesh.elements.sfncon6(i).master = mesh.elements.tria6(grp_idx_master, :);
%!     mesh.elements.sfncon6(i).maxdist = sqrt(eps) * max(abs([a,b,c]));
%!   endfor
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], ...
%!                           load_case);
%!   assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);

%!   for i=1:numel(sol_eig.f)
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!   endfor
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
%!   for i=1:length(sol_eig.f)
%!     fprintf(stderr, "mode %d f=%.0f\n", i, sol_eig.f(i));
%!   endfor
%!   assert(sol_eig.f(:), f_ref, tol * max(f_ref));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 15e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 3.5e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   if (do_rotate)
%!     fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!   endif
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
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
%!   endif
%!   mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%!   mesh1.material_data.E = 210000e6;
%!   mesh1.material_data.nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%!   mesh2 = mesh1;
%!   mesh2.nodes(:, 1) += a;
%!   for i=1:numel(mesh2.groups.tria6)
%!     mesh2.groups.tria6(i).id += 100;
%!   endfor
%!   data(1).mesh = mesh1;
%!   data(2).mesh = mesh2;
%!   [mesh] = fem_post_mesh_merge(data);
%!   mesh.elements.sfncon6.slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 2)).nodes(:);
%!   mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 101)).elements, :);
%!   mesh.elements.sfncon6.maxdist = sqrt(eps) * max(abs([a,b,c]));
%!   mesh.elements.sfncon6.constraint = FEM_CT_SLIDING;
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 102)).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], ...
%!                           load_case);
%!   assert(mtot, a * b * c * sum([mesh.material_data.rho]), sqrt(eps) * a * b * c * sum([mesh.material_data.rho]));
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);

%!   for i=1:numel(sol_eig.f)
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!   endfor
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect


%!test
%! ## TEST10
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   num_modes = 6;
%!   shift = 0;
%!   scale_eig = 1.5e-3;
%!   E = [210000e6; 70000e6];
%!   nu = [0.3, 0.3];
%!   rho = [7850; 2700];
%!   ri = [4e-3; 2.5e-3];
%!   ro = [5e-3; ri(1)];
%!   h = 12e-3;
%!   scale_def = 5e-3;
%!   toldist = 1e-3;
%!   mesh_size = linspace(1.25e-3, 0.4e-3, 10)(end);
%!   num_nodes = zeros(1, numel(mesh_size));
%!   assert(numel(ri), numel(ro));
%!   f = nan(numel(mesh_size), num_modes);
%!   for m=1:numel(mesh_size)
%!     clear data;
%!     clear mesh;
%!     clear mat_ass;
%!     clear sol_eig;
%!     for i=1:numel(ri)
%!       fd = -1;
%!       unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "wt");
%!       if (fd == -1)
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
%!       fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; };\n");
%!       fprintf(fd, "Physical Volume(\"volume\",%d) = {tmp[1]};\n", (i - 1) * 100 + 1);
%!       fprintf(fd, "Physical Surface(\"bottom\",%d) = {tmp[2]};\n", (i - 1) * 100 + 1);
%!       fprintf(fd, "Physical Surface(\"outside\",%d) = {tmp[3]};\n", (i - 1) * 100 + 2);
%!       fprintf(fd, "Physical Surface(\"inside\",%d) = {tmp[5]};\n", (i - 1) * 100 + 3);
%!       fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[4]};\n", (i - 1) * 100 + 4);
%!       fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[0]};\n", (i - 1) * 100 + 5); ##x
%!       fprintf(fd, "Physical Surface(\"right\",%d) = {6};\n", (i - 1) * 100 + 6);  ##y
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       unlink([filename, ".msh"]);
%!       pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size(m)), "-clmax", sprintf("%g", 1.25 *mesh_size(m)), [filename, ".geo"]});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         warning("gmsh failed with status %d", status);
%!       endif
%!       unlink([filename, ".geo"]);
%!       data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!       unlink([filename, ".msh"]);
%!     endfor
%!     for i=1:numel(data)
%!       data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!       data(i).mesh.material_data.rho = rho(i);
%!       data(i).mesh.material_data.C = fem_pre_mat_isotropic(E(i), nu(i));
%!     endfor
%!     mesh = fem_post_mesh_merge(data);
%!     mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(find([mesh.groups.tria6.id]==3)).elements, :);
%!     mesh.elements.sfncon6.slave = mesh.groups.tria6(find([mesh.groups.tria6.id]==102)).nodes(:);
%!     mesh.elements.sfncon6.maxdist = toldist * ri(1);
%!     mesh.elements.sfncon6.constraint = FEM_CT_SLIDING;
%!     group_id = [[mesh.groups.tria6].id];
%!     node_constr1 = [mesh.groups.tria6(find((group_id == 1)|(group_id == 101))).nodes];
%!     node_constr5 = [mesh.groups.tria6(find(((group_id == 5)|(group_id == 105)))).nodes];
%!     node_constr6 = [mesh.groups.tria6(find(((group_id == 6)|(group_id == 106)))).nodes];
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
%!     sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, shift);
%!     f(m, :) = sol_eig.f;
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

%!     fref = [52248
%!             65901
%!             105250
%!             113070
%!             132630
%!             201780];
%!     assert(sol_eig.f(:), fref, 0.1e-2 * max(fref));

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
%!     plot(num_nodes.', f(:, i), sprintf("-;mode %d;%d", i, i));
%!   endfor
%!   xlabel("nodes [1]");
%!   ylabel("f [Hz]");
%!   grid on;
%!   grid minor on;
%!   title("natural frequencies versus mesh size");
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST11
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   num_modes = 6;
%!   shift = 0;
%!   scale_eig = 1.5e-3;
%!   E = [210000e6; 70000e6];
%!   nu = [0.3, 0.3];
%!   rho = [7850; 2700];
%!   ri = [8e-3; 5e-3];
%!   ro = [10e-3; ri(1)];
%!   h = 12e-3;
%!   scale_def = 5e-3;
%!   toldist = 1e-3;
%!   tolf = 8e-2;
%!   mesh_size = linspace(0.6e-3, 0.3e-3, 5)(1);
%!   num_nodes = zeros(1, numel(mesh_size));
%!   assert(numel(ri), numel(ro));
%!   orange = [2,1];
%!   rrange = [false,true];
%!   f = nan(numel(mesh_size), num_modes, numel(rrange), numel(orange));
%!   for o=1:numel(orange)
%!     for r=1:numel(rrange)
%!       if (orange(o) == 2 && rrange(r))
%!         continue;
%!       endif
%!       for m=1:numel(mesh_size)
%!         clear data;
%!         clear mesh;
%!         clear mat_ass;
%!         clear sol_eig;
%!         for i=1:numel(ri)
%!           fd = -1;
%!           unwind_protect
%!           [fd, msg] = fopen([filename, ".geo"], "wt");
%!           if fd == -1
%!             error("failed to open file \"%s.geo\"", filename);
%!           endif
%!           fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!           fprintf(fd, "ri = %g;\n", ri(i));
%!           fprintf(fd, "ro = %g;\n", ro(i));
%!           fprintf(fd, "h = %g;\n", h);
%!           fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!           fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!           fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!           fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!           fputs(fd, "Line(1) = {1,2};\n");
%!           fputs(fd, "Line(4) = {2,5};\n");
%!           fputs(fd, "Line(5) = {5,6};\n");
%!           fputs(fd, "Line(8) = {6,1};\n");
%!           fputs(fd, "Line Loop(5) = {1,4,5,8};\n");
%!           fputs(fd, "Plane Surface(6) = {5};\n");
%!           if (orange(o) == 1 || i == 1)
%!             if (rrange(r) || i == 1)
%!               fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{%d}; Recombine; };\n", ceil(0.5 * pi * ri(1) / mesh_size(m)));
%!               fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!             else
%!               fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{%d};};\n", ceil(0.5 * pi * ri(1) / mesh_size(m)));
%!             endif
%!           else
%!             fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6};};\n");
%!           endif
%!           fprintf(fd, "Physical Volume(\"volume\",%d) = {tmp[1]};\n", (i - 1) * 100 + 1);
%!           fprintf(fd, "Physical Surface(\"bottom\",%d) = {tmp[2]};\n", (i - 1) * 100 + 1);
%!           fprintf(fd, "Physical Surface(\"outside\",%d) = {tmp[3]};\n", (i - 1) * 100 + 2);
%!           fprintf(fd, "Physical Surface(\"inside\",%d) = {tmp[5]};\n", (i - 1) * 100 + 3);
%!           fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[4]};\n", (i - 1) * 100 + 4);
%!           fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[0]};\n", (i - 1) * 100 + 5); ##x
%!           fprintf(fd, "Physical Surface(\"right\",%d) = {6};\n", (i - 1) * 100 + 6);  ##y
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!           end_unwind_protect
%!           unlink([filename, ".msh"]);
%!           if (orange(o) == 1 || i == 1)
%!             optargs = {"-order", "1"};
%!           else
%!             optargs = {"-order", "2"};
%!           endif
%!           pid = spawn("gmsh", {"-format", "msh2", "-3", optargs{:}, "-clmin", sprintf("%g", 0.75 * mesh_size(m) * orange(o)), "-clmax", sprintf("%g", 1.25 *mesh_size(m) * orange(o)), [filename, ".geo"]});
%!           status = spawn_wait(pid);
%!           if status ~= 0
%!             warning("gmsh failed with status %d", status);
%!           endif
%!           unlink([filename, ".geo"]);
%!           data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!           unlink([filename, ".msh"]);
%!         endfor
%!         for i=1:numel(data)
%!           if (isfield(data(i).mesh.elements, "iso8"))
%!             data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!           else
%!             data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!           endif
%!           data(i).mesh.material_data.rho = rho(i);
%!           data(i).mesh.material_data.C = fem_pre_mat_isotropic(E(i), nu(i));
%!         endfor
%!         mesh = fem_post_mesh_merge(data);
%!         mesh.elements.sfncon4.master = mesh.elements.iso4(mesh.groups.iso4(find([mesh.groups.iso4.id]==3)).elements, :);
%!         if (orange(o) == 1)
%!           mesh.elements.sfncon4.slave = mesh.groups.iso4(find([mesh.groups.iso4.id]==102)).nodes(:);
%!         else
%!           mesh.elements.sfncon4.slave = mesh.groups.tria6(find([mesh.groups.tria6.id]==102)).nodes(:);
%!         endif
%!         mesh.elements.sfncon4.maxdist = toldist * ri(1);
%!         mesh.elements.sfncon4.constraint = FEM_CT_SLIDING;
%!         if (orange(o) == 1)
%!           group_id = [[mesh.groups.iso4].id];
%!           node_constr1 = [mesh.groups.iso4(find((group_id == 1)|(group_id == 101))).nodes];
%!           node_constr5 = [mesh.groups.iso4(find(((group_id == 5)|(group_id == 105)))).nodes];
%!           node_constr6 = [mesh.groups.iso4(find(((group_id == 6)|(group_id == 106)))).nodes];
%!         else
%!           group_id4 = [[mesh.groups.iso4].id];
%!           group_id6 = [[mesh.groups.tria6].id];
%!           node_constr1 = [mesh.groups.iso4(find((group_id4 == 1))).nodes,   mesh.groups.tria6(find((group_id6 == 101))).nodes];
%!           node_constr5 = [mesh.groups.iso4(find(((group_id4 == 5)))).nodes, mesh.groups.tria6(find(((group_id6 == 105)))).nodes];
%!           node_constr6 = [mesh.groups.iso4(find(((group_id4 == 6)))).nodes, mesh.groups.tria6(find(((group_id6 == 106)))).nodes];
%!         endif
%!         node_constr = [node_constr1, node_constr5, node_constr6];
%!         mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!         idx_joint = int32(0);
%!         for i=1:numel(node_constr1)
%!           mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!           mesh.elements.joints(idx_joint).nodes = node_constr1(i);
%!         endfor
%!         for i=1:numel(node_constr5)
%!           mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!           mesh.elements.joints(idx_joint).nodes = node_constr5(i);
%!         endfor
%!         for i=1:numel(node_constr6)
%!           mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!           mesh.elements.joints(idx_joint).nodes = node_constr6(i);
%!         endfor
%!         load_case.locked_dof = false(rows(mesh.nodes), 6);
%!         dof_map = fem_ass_dof_map(mesh, load_case);
%!         num_nodes(m) = rows(mesh.nodes);
%!         [mat_ass.K, ...
%!          mat_ass.M, ...
%!          mat_ass.mtot] = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_MAT_STIFFNESS, ...
%!                                          FEM_MAT_MASS, ...
%!                                          FEM_SCA_TOT_MASS], ...
%!                                         load_case);
%!         sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes, shift);
%!         fem_post_mesh_export([filename, ".msh"], mesh);
%!         f(m, :, r, o) = sol_eig.f;
%!         figure("visible", "off");
%!         for i=1:numel(data)
%!           fem_post_sol_plot(data(i).mesh);
%!         endfor
%!         xlabel("x [m]");
%!         ylabel("y [m]");
%!         zlabel("z [m]");
%!         grid on;
%!         grid minor on;
%!         title("undeformed mesh");

%!         for i=1:numel(sol_eig.f)
%!           figure("visible", "off");
%!           hold on;
%!           fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, 1:3, i), "rows")), i);
%!           view(30,30);
%!           xlabel('x [m]');
%!           ylabel('y [m]');
%!           zlabel('z [m]');
%!           grid on;
%!           grid minor on;
%!           title(sprintf("mode %d: %.1fHz", i, sol_eig.f(i)));
%!         endfor
%!       endfor
%!       figure("visible", "off");
%!       hold on;
%!       for i=1:columns(f)
%!         plot(num_nodes.', f(:, i, r, o), sprintf("-;mode %d;%d", i, i));
%!       endfor
%!       xlabel("nodes [1]");
%!       ylabel("f [Hz]");
%!       grid on;
%!       grid minor on;
%!       title("natural frequencies versus mesh size");
%!     endfor
%!   endfor
%!   fref = [26023  59514  91469  1.0372e+05  1.1294e+05  1.154e+05];
%!   for o=1:numel(orange)
%!     for r=1:numel(rrange)
%!       for i=1:rows(f)
%!         if (~all(isnan(f(i, :, r, o))))
%!           assert(all(abs(f(i, :, r, o) ./ fref - 1) < tolf));
%!         endif
%!       endfor
%!     endfor
%!   endfor
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   for i=1:4
%!     fd = -1;
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "wt");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       switch (i)
%!         case {1,2}
%!           order = 1;
%!         otherwise
%!           order = 2;
%!       endswitch
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "a=%g;\n", a);
%!       fprintf(fd, "b=%g;\n", b);
%!       fprintf(fd, "c=%g;\n", c);
%!       fprintf(fd, "h=%g;\n", h * order);
%!       fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!       fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!       fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!       fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,1};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp[] = Extrude {0.0,0.0,c} {\n");
%!       switch (i)
%!         case 1
%!           fprintf(fd, "  Surface{6}; Layers{%d}; Recombine;\n", ceil(c/h));
%!         otherwise
%!           fputs(fd, "  Surface{6};\n");
%!       endswitch
%!       fputs(fd, "};\n");
%!       switch (i)
%!         case 1
%!           fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!       endswitch
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!       fprintf(fd, "Physical Surface(\"right\",%d) = {tmp[2]};\n", 100 * i + 1);
%!       fprintf(fd, "Physical Surface(\"rear\",%d) = {tmp[3]};\n", 100 * i + 2);
%!       fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[4]};\n", 100 * i + 3);
%!       fprintf(fd, "Physical Surface(\"front\",%d) = {tmp[5]};\n", 100 * i + 4);
%!       fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[0]};\n", 100 * i + 5);
%!       fprintf(fd, "Physical Surface(\"bottom\",%d) = {6};\n", 100 * i + 6);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", sprintf("%d", order), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     unlink([filename, ".geo"]);
%!     data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!     unlink([filename, ".msh"]);
%!   endfor
%!   data(2).mesh.nodes(:, 1) += a;
%!   data(3).mesh.nodes(:, 1) += a;
%!   data(3).mesh.nodes(:, 2) += b;
%!   data(4).mesh.nodes(:, 2) += b;
%!   for i=1:numel(data)
%!     if (isfield(data(i).mesh.elements, "iso8"))
%!       data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!     else
%!       data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!     endif
%!     data(i).mesh.material_data.rho = rho;
%!     data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   endfor
%!   mesh = fem_post_mesh_merge(data);
%!   constr = FEM_CT_FIXED;
%!   mesh.elements.sfncon4(1).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 103)).elements, :);
%!   mesh.elements.sfncon4(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 401)).nodes(:);
%!   mesh.elements.sfncon4(1).maxdist = maxdist;
%!   mesh.elements.sfncon4(1).constraint = constr;
%!   mesh.elements.sfncon4(2).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 102)).elements, :);
%!   mesh.elements.sfncon4(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 204)).nodes(:);
%!   mesh.elements.sfncon4(2).maxdist = maxdist;
%!   mesh.elements.sfncon4(2).constraint = constr;
%!   mesh.elements.sfncon6(1).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 402)).elements, :);
%!   mesh.elements.sfncon6(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 304)).nodes(:);
%!   mesh.elements.sfncon6(1).maxdist = maxdist;
%!   mesh.elements.sfncon6(1).constraint = constr;
%!   mesh.elements.sfncon6(2).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 301)).elements, :);
%!   mesh.elements.sfncon6(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 203)).nodes(:);
%!   mesh.elements.sfncon6(2).maxdist = maxdist;
%!   mesh.elements.sfncon6(2).constraint = constr;
%!   slave_nodes = [];
%!   for i=1:numel(mesh.elements.sfncon4)
%!     idx_slave = [];
%!     for j=1:numel(slave_nodes)
%!       idx_slave = [idx_slave; find(mesh.elements.sfncon4(i).slave == slave_nodes(j))];
%!     endfor
%!     if (numel(idx_slave))
%!       mesh.elements.sfncon4(i).slave(idx_slave) = 0;
%!       mesh.elements.sfncon4(i).slave = mesh.elements.sfncon4(i).slave(find(mesh.elements.sfncon4(i).slave));
%!     endif
%!     slave_nodes = [slave_nodes; mesh.elements.sfncon4(i).slave];
%!   endfor

%!   for i=1:numel(mesh.elements.sfncon6)
%!     idx_slave = [];
%!     for j=1:numel(slave_nodes)
%!       idx_slave = [idx_slave; find(mesh.elements.sfncon6(i).slave == slave_nodes(j))];
%!     endfor
%!     if (numel(idx_slave))
%!       mesh.elements.sfncon6(i).slave(idx_slave) = 0;
%!       mesh.elements.sfncon6(i).slave = mesh.elements.sfncon6(i).slave(find(mesh.elements.sfncon6(i).slave));
%!     endif
%!     slave_nodes = [slave_nodes; mesh.elements.sfncon6(i).slave];
%!   endfor
%!   group_id4 = [[mesh.groups.iso4].id];
%!   group_id6 = [[mesh.groups.tria6].id];
%!   node_bottom = [[mesh.groups.iso4(find(mod(group_id4, 100) == 6))].nodes,   [mesh.groups.tria6(find(mod(group_id6,100)==6))].nodes];
%!   node_front = [[mesh.groups.iso4(find(group_id4 == 104))].nodes, [mesh.groups.tria6(find(group_id6 == 404))].nodes];
%!   node_right = [mesh.groups.iso4(find((group_id4 == 101) | (group_id4 == 201))).nodes];
%!   node_constr = [node_bottom, node_front, node_right];
%!   mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!   idx_joint = int32(0);
%!   for i=1:numel(node_bottom)
%!     if (~numel(find(slave_nodes == node_bottom(i))))
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_bottom(i);
%!     endif
%!   endfor
%!   for i=1:numel(node_front)
%!     if (~numel(find(slave_nodes == node_front(i))))
%!       mesh.elements.joints(++idx_joint).C = [1,0,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_front(i);
%!     endif
%!   endfor
%!   for i=1:numel(node_right)
%!     if (~numel(find(slave_nodes == node_right(i))))
%!       mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_right(i);
%!     endif
%!   endfor

%!   mesh.elements.joints = mesh.elements.joints(1:idx_joint);

%!   iso4_top = [[mesh.groups.iso4(find(mod(group_id4, 100) == 5))].elements];
%!   tria6_top = [[mesh.groups.tria6(find(mod(group_id6, 100) == 5))].elements];
%!   iso4_rear = mesh.groups.iso4(find(group_id4 == 202)).elements;
%!   tria6_rear = mesh.groups.tria6(find(group_id6 == 302)).elements;
%!   tria6_left = [[mesh.groups.tria6(find((group_id6 == 303) | (group_id6 == 403)))].elements];
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.pressure.tria6.elements = mesh.elements.tria6([tria6_left, tria6_rear, tria6_top], :);
%!   load_case.pressure.tria6.p = [repmat(py, numel(tria6_left), 6); repmat(px, numel(tria6_rear), 6); repmat(pz, numel(tria6_top), 6)];
%!   load_case.pressure.iso4.elements = mesh.elements.iso4([iso4_rear, iso4_top], :);
%!   load_case.pressure.iso4.p = [repmat(px, numel(iso4_rear), 4); repmat(pz, numel(iso4_top), 4)];
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_STIFFNESS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case);
%!   if (eliminate)
%!     [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%!     opt_ls.refine_max_iter = int32(100);
%!     Kfact = fem_sol_factor(Kred, opt_ls);
%!     Ured = Kfact \ Rred;
%!     sol_stat.def = fem_post_def_nodal(mesh, dof_map, Tred * Ured);
%!   else
%!     sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   endif
%!   sol_stat.F = fem_post_def_nodal(mesh, dof_map, mat_ass.R);
%!   Ftot = sum(sol_stat.F, 1);
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh);
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("undeformed mesh");
%!   figure("visible", "off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat, scale/max(norm(sol_stat.def, "rows")), 1);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title("deformed mesh");
%!   sigma_a = [-px; -py; -pz; zeros(3, 1)];
%!   epsilon_a = mesh.material_data(1).C \ sigma_a;
%!   U_a = zeros(rows(mesh.nodes), 3);
%!   for i=1:3
%!     U_a(:, i) = epsilon_a(i) * mesh.nodes(:, i);
%!   endfor
%!   Ftot_a = [-2 * b * c * px, -2 * a * c * py, -4 * a * b * pz, zeros(1, 3)];
%!   fprintf(stderr, "max(err)=%g\n", max(max(abs(sol_stat.def(:, 1:3) - U_a))) / max(max(abs(U_a))));
%!   fprintf(stderr, "mean(err)=%g\n", mean(mean(abs(sol_stat.def(:, 1:3) - U_a))) / max(max(abs(U_a))));
%!   assert(sol_stat.def(:, 1:3), U_a, tolstat * max(max(abs(U_a))));
%!   assert(Ftot, Ftot_a, sqrt(eps) * norm(Ftot_a));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ### TEST13
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 60e-3;
%!   c = 10e-3;
%!   d = 0e-3;
%!   h = 2e-3;
%!   b = h;
%!   Fx = 0;
%!   Fz = -1000;
%!   My = 0;
%!   num_iso = 10;
%!   do_post_pro = false;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0, -0.5 * b, -0.5 * c, h};\n");
%!   fputs(fd, "Point(2) = {  a, -0.5 * b, -0.5 * c, h};\n");
%!   fputs(fd, "Point(3) = {  a,  0.5 * b, -0.5 * c, h};\n");
%!   fputs(fd, "Point(4) = {0.0,  0.5 * b, -0.5 * c, h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   master_node_idx = int32(rows(mesh.nodes) + 1);
%!   mesh.nodes(master_node_idx, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, master_node_idx);
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   load_case.loaded_nodes = [master_node_idx];
%!   load_case.loads = [Fx, 0, Fz, 0, My, 0];
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.dm] = fem_ass_matrix(mesh, ...
%!                                 dof_map, ...
%!                                 [FEM_MAT_STIFFNESS, ...
%!                                  FEM_VEC_LOAD_CONSISTENT, ...
%!                                  FEM_SCA_TOT_MASS], ...
%!                                 load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   x = linspace(0, a, 100);
%!   z = linspace(-0.5 * c, 0.5 * c, 50);
%!   [xx, zz] = meshgrid(x, z);
%!   xtauel = mesh.nodes(:, 1)(mesh.elements.tet10);
%!   ytauel = mesh.nodes(:, 1)(mesh.elements.tet10);
%!   ztauel = mesh.nodes(:, 3)(mesh.elements.tet10);
%!   tauxxel = sol_stat.stress.taum.tet10(:, :, 1);
%!   tauxzel = sol_stat.stress.taum.tet10(:, :, 6);
%!   tauxx = griddata(xtauel(:), ztauel(:), tauxxel(:), xx, zz);
%!   tauxz = griddata(xtauel(:), ztauel(:), tauxzel(:), xx, zz);
%!   Iy = b * c^3 / 12;
%!   tauxx_a = -Fz / Iy * (a - xx) .* zz + Fx / (b * c) + My / Iy * zz;
%!   tauxz_a = 3 / 2 * Fz / (b * c) * (1 - (zz / (0.5 * c)).^2);
%!   scale_tauxx = linspace(min(min(tauxx_a)), max(max(tauxx_a)), num_iso + 1);
%!   figure("visible", "off");
%!   subplot(2, 1, 1);
%!   contourf(xx, zz, tauxx_a, scale_tauxx);
%!   daspect([1,1,1]);
%!   colormap jet;
%!   colorbar;
%!   xlabel("x [m]");
%!   ylabel("z [m]");
%!   title("stress component tauxx [Pa]");
%!   grid on;
%!   grid minor on;
%!   subplot(2, 1, 2);
%!   contourf(xx, zz, tauxx, scale_tauxx);
%!   daspect([1,1,1]);
%!   colormap jet;
%!   colorbar;
%!   xlabel("x [m]");
%!   ylabel("z [m]");
%!   title("stress component tauxx [Pa]");
%!   grid on;
%!   grid minor on;
%!   scale_tauxz = linspace(min(min(tauxz_a)), max(max(tauxz_a)), num_iso + 1);
%!   figure("visible", "off");
%!   subplot(2, 1, 1);
%!   contourf(xx, zz, tauxz_a, scale_tauxz);
%!   daspect([1,1,1]);
%!   colormap jet;
%!   colorbar;
%!   xlabel("x [m]");
%!   ylabel("z [m]");
%!   title("stress component tauxz [Pa]");
%!   grid on;
%!   grid minor on;
%!   subplot(2, 1, 2);
%!   contourf(xx, zz, tauxz, scale_tauxz);
%!   daspect([1,1,1]);
%!   colormap jet;
%!   colorbar;
%!   xlabel("x [m]");
%!   ylabel("z [m]");
%!   title("stress component tauxz [Pa]");
%!   grid on;
%!   grid minor on;

%!   figure_list();
%!   if (do_post_pro)
%!     opts.scale_def = 0.3 * a / max(max(abs(sol_stat.def)));
%!     fem_post_sol_external(mesh, sol_stat, opts);
%!   endif
%!   idx_x = find((xx(:) > 0.1 * a) & (xx(:) < 0.9 * a));
%!   assert(tauxx(:)(idx_x), tauxx_a(:)(idx_x), 0.3e-2 * max(max(max(abs(sol_stat.stress.taum.tet10)))));
%!   assert(tauxz(:)(idx_x), tauxz_a(:)(idx_x), 0.3e-2 * max(max(max(abs(sol_stat.stress.taum.tet10)))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ### TEST14
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   N = 3;
%!   a = 10e-3;
%!   b = 5e-3;
%!   c = 5e-3;
%!   h = 5e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6};\n");
%!   fprintf(fd, "  Layers{%d}; Recombine;\n", ceil(c / h));
%!   fputs(fd, "};\n");
%!   fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh1 = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   mesh1.materials.iso8 = ones(rows(mesh1.elements.iso8), 1, "int32");
%!   mesh1.material_data.E = 210000e6;
%!   mesh1.material_data.nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(mesh1.material_data.E, mesh1.material_data.nu);
%!   data.mesh = mesh1;
%!   data = repmat(data, 1, 3);
%!   for i=2:3
%!     data(i).mesh.nodes(:, 1) += a * (i - 1);
%!     for j=1:numel(data(i).mesh.groups.iso4)
%!       data(i).mesh.groups.iso4(j).id += 100 * (i - 1);
%!       data(i).mesh.groups.iso4(j).name = sprintf("%s[%d]", data(i).mesh.groups.iso4(j).name, data(i).mesh.groups.iso4(j).id);
%!     endfor
%!     for j=1:numel(data(i).mesh.groups.iso8)
%!       data(i).mesh.groups.iso8(j).id += 100 * (i - 1);
%!       data(i).mesh.groups.iso8(j).name = sprintf("%s[%d]", data(i).mesh.groups.iso8(j).name, data(i).mesh.groups.iso8(j).id);
%!     endfor
%!   endfor
%!   [mesh] = fem_post_mesh_merge(data);
%!   mesh.nodes(:, 2) -= 0.5 * b;
%!   mesh.nodes(:, 3) -= 0.5 * c;
%!   for i=1:N
%!     if (i > 1)
%!       unwind_protect
%!         fem_post_mesh_export([filename, "_in.msh"], data(i - 1).mesh);
%!         pid = spawn("gmsh", {"-refine", "-format", "msh2", "-o", [filename, "_out.msh"], [filename, "_in.msh"]});
%!         status = spawn_wait(pid);
%!         if (status ~= 0)
%!           error("gmsh failed with status %d", status);
%!         endif
%!         data(i).mesh = fem_pre_mesh_import([filename, "_out.msh"]);
%!         data(i).mesh.material_data = mesh.material_data;
%!         data(i).mesh.materials.iso8 = zeros(rows(data(i).mesh.elements.iso8), 1, "int32");
%!         for j=1:numel(data(i).mesh.groups.iso8)
%!           data(i).mesh.materials.iso8(data(i).mesh.groups.iso8(j).elements) = j;
%!         endfor
%!       unwind_protect_cleanup
%!         unlink([filename, "_in.msh"]);
%!         unlink([filename, "_out.msh"]);
%!       end_unwind_protect
%!     else
%!       data(i).mesh = mesh;
%!     endif
%!     data(i).h = h / i;
%!     node_idx_rbe3 = int32(rows(data(i).mesh.nodes) + 1);
%!     data(i).mesh.nodes(node_idx_rbe3, 1:3) = [3 * a, 0, 0];
%!     grp_idx = find([data(i).mesh.groups.iso4.id] == 202);
%!     data(i).mesh.elements.rbe3.nodes = [node_idx_rbe3, data(i).mesh.groups.iso4(grp_idx).nodes];
%!     data(i).mesh.elements.rbe3.weight = ones(numel(data(i).mesh.elements.rbe3.nodes) - 1, 1);
%!     for j=1:2
%!       grp_idx_slave = find([[data(i).mesh.groups.iso4].id] == (j - 1) * 100 + 2);
%!       grp_idx_master = data(i).mesh.groups.iso4(find([[data(i).mesh.groups.iso4].id] == j * 100 + 1)).elements;
%!       data(i).mesh.elements.sfncon4(j).slave = data(i).mesh.groups.iso4(grp_idx_slave).nodes(:);
%!       data(i).mesh.elements.sfncon4(j).master = data(i).mesh.elements.iso4(grp_idx_master, :);
%!       data(i).mesh.elements.sfncon4(j).maxdist = sqrt(eps) * max(abs([a,b,c]));
%!     endfor

%!     data(i).load_case.locked_dof = false(rows(data(i).mesh.nodes), 6);
%!     data(i).load_case.locked_dof(data(i).mesh.groups.iso4(find([[data(i).mesh.groups.iso4].id] == 1)).nodes, :) = true;
%!     data(i).load_case.loaded_nodes = node_idx_rbe3;
%!     data(i).load_case.loads = [0,-1000, 0, 0, 0, 0];
%!     data(i).dof_map = fem_ass_dof_map(data(i).mesh, data(i).load_case);
%!     [data(i).mat_ass.K, ...
%!      data(i).mat_ass.R, ...
%!      data(i).mat_ass.mtot] = fem_ass_matrix(data(i).mesh, ...
%!                                             data(i).dof_map, ...
%!                                             [FEM_MAT_STIFFNESS, ...
%!                                              FEM_VEC_LOAD_CONSISTENT, ...
%!                                              FEM_SCA_TOT_MASS], ...
%!                                             data(i).load_case);
%!     assert(data(i).mat_ass.mtot, ...
%!            a * b * c * sum([data(i).mesh.material_data.rho]), ...
%!            sqrt(eps) * a * b * c * sum([data(i).mesh.material_data.rho]));
%!     [data(i).sol_stat, data(i).sol_stat.U] = fem_sol_static(data(i).mesh, data(i).dof_map, data(i).mat_ass);
%!     data(i).sol_stat.stress = fem_ass_matrix(data(i).mesh, ...
%!                                              data(i).dof_map, ...
%!                                              [FEM_VEC_STRESS_CAUCH], ...
%!                                              data(i).load_case, ...
%!                                              data(i).sol_stat);
%!     data(i).W = data(i).sol_stat.U.' * data(i).mat_ass.K * data(i).sol_stat.U;
%!   endfor
%!   figure("visible", "off");
%!   loglog([data.h], [data.W], "-x;W(h) [J];1");
%!   xlabel("h [m]");
%!   ylabel("W [J]");
%!   grid on;
%!   grid minor on;
%!   title("strain energy");
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ### TEST 15
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   h1 = 10e-3;
%!   h2 = 2e-3;
%!   w = 4e-3;
%!   l = 30e-3;
%!   h = 2e-3;
%!   do_post_pro = false;
%!   scale_def = 20e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "h1=%g;\n", h1);
%!   fprintf(fd, "h2=%g;\n", h2)
%!   fprintf(fd, "w=%g;\n", w);
%!   fprintf(fd, "l=%g;\n", l);
%!   fprintf(fd, "h=%g;\n", h);
%!   fputs(fd, "Point(1) = {0.0, -0.5 * w, -0.5 * h1 - h2, h};\n");
%!   fputs(fd, "Point(2) = {0.0, -0.5 * w, -0.5 * h1, h};\n");
%!   fputs(fd, "Point(3) = {0.0,  0.5 * w, -0.5 * h1, h};\n");
%!   fputs(fd, "Point(4) = {0.0,  0.5 * w, -0.5 * h1 - h2, h};\n");
%!   fputs(fd, "Point(5) = {0.0, -0.5 * w, 0.5 * h1, h};\n");
%!   fputs(fd, "Point(6) = {0.0, -0.5 * w, 0.5 * h1 + h2, h};\n");
%!   fputs(fd, "Point(7) = {0.0,  0.5 * w, 0.5 * h1 + h2, h};\n");
%!   fputs(fd, "Point(8) = {0.0,  0.5 * w, 0.5 * h1, h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line(5) = {2,5};\n");
%!   fputs(fd, "Line(6) = {5,8};\n");
%!   fputs(fd, "Line(7) = {8,3};\n");
%!   fputs(fd, "Line(8) = {3,2};\n");
%!   fputs(fd, "Line(9) = {5,6};\n");
%!   fputs(fd, "Line(10) = {6,7};\n");
%!   fputs(fd, "Line(11) = {7,8};\n");
%!   fputs(fd, "Line(12) = {8,5};\n");
%!   fputs(fd, "Line Loop(13) = {1,2,3,4};\n");
%!   fputs(fd, "Line Loop(14) = {5,6,7,8};\n");
%!   fputs(fd, "Line Loop(15) = {9,10,11,12};\n");
%!   fputs(fd, "Plane Surface(16) = {13};\n");
%!   fputs(fd, "Plane Surface(17) = {14};\n");
%!   fputs(fd, "Plane Surface(18) = {15};\n");
%!   fputs(fd, "v1[] = Extrude {l,0.0,0.0} {\n");
%!   fputs(fd, "  Surface{16};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "v2[] = Extrude {l,0.0,0.0} {\n");
%!   fputs(fd, "  Surface{17};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "v3[] = Extrude {l,0.0,0.0} {\n");
%!   fputs(fd, "  Surface{18};\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "v = newv;\n");
%!   fputs(fd, "BooleanUnion(v) = {Volume{v1[1]};}{ Volume{v2[1],v3[1]};};\n");
%!   fputs(fd, "Coherence;\n");
%!   fputs(fd, "Physical Volume(\"volume1\",1) = {v1[1]};\n");
%!   fputs(fd, "Physical Volume(\"volume2\",2) = {v2[1]};\n");
%!   fputs(fd, "Physical Volume(\"volume3\",3) = {v3[1]};\n");
%!   fputs(fd, "Physical Surface(\"surface1\",1) = {v1[0]};\n");
%!   fputs(fd, "Physical Surface(\"surface2\",2) = {v2[0]};\n");
%!   fputs(fd, "Physical Surface(\"surface3\",3) = {v3[0]};\n");
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
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   mesh.material_data(1).E = 80000e6;
%!   mesh.material_data(1).nu = 0.4;
%!   mesh.material_data(1).rho = 1700;
%!   mesh.material_data(1).alpha = 1e-8;
%!   mesh.material_data(1).beta = 1e-6;
%!   mesh.material_data(2).E = 10000e6;
%!   mesh.material_data(2).nu = 0.4;
%!   mesh.material_data(2).rho = 500;
%!   mesh.material_data(2).alpha = 0;
%!   mesh.material_data(2).beta = 0;
%!   for i=1:numel(mesh.material_data)
%!     mesh.material_data(i).C = fem_pre_mat_isotropic(mesh.material_data(i).E, mesh.material_data(i).nu);
%!   endfor
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10),1);
%!   mesh.materials.tet10(mesh.groups.tet10(1).elements) = 1;
%!   mesh.materials.tet10(mesh.groups.tet10(2).elements) = 2;
%!   mesh.materials.tet10(mesh.groups.tet10(3).elements) = 1;
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(find(mesh.nodes(:, 1)==0), 1:3) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.M, ...
%!    mat_ass.D, ...
%!    mat_ass.K, ...
%!    mat_ass.dm, ...
%!    mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [FEM_MAT_MASS, ...
%!                                        FEM_MAT_DAMPING, ...
%!                                        FEM_MAT_STIFFNESS, ...
%!                                        FEM_SCA_TOT_MASS], ...
%!                                       load_case);
%!   [sol_eig, err] = fem_sol_modal(mesh, dof_map, mat_ass, 10);
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig);
%!   opts.scale_def = scale_def / max(max(max(abs(sol_eig.def))));
%!   rho1 = mesh.material_data(1).rho;
%!   rho2 = mesh.material_data(2).rho;
%!   m1 = rho1 * h2 * w * l;
%!   m2 = rho2 * h1 * w * l;
%!   m = 2 * m1 + m2;
%!   assert(mat_ass.dm, m, sqrt(eps) * m);
%!   opts.skin_only = true;
%!   if (do_post_pro)
%!     fem_post_sol_external(mesh, sol_eig, opts);
%!   endif
%! unwind_protect_cleanup
%!   unlink([filename, ".msh"]);
%!   unlink([filename, ".geo"]);
%! end_unwind_protect

%!test
%! ## TEST 16
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   ## K.J.Bathe page 328 4.20a
%!   mesh_size = 1.5;
%!   p1 = 0.006;
%!   E = 55;
%!   nu = 0.3;
%!   rho = 1000e-12;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%!   fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Circle(3) = {3,11,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Circle(7) = {7,12,8};\n");
%!   fputs(fd, "Line(8) = {8,9};\n");
%!   fputs(fd, "Line(9) = {9,10};\n");
%!   fputs(fd, "Line(10) = {10,1};\n");
%!   fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%!   fputs(fd, "Plane Surface(14) = {11};\n");
%!   fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; };\n", mesh_size);
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%!   fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%!   fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%!   fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%!   grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!   elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%!   elno_p1 = mesh.elements.tria6(elem_id_p1, :);

%!   load_case.pressure.tria6.elements = elno_p1;
%!   load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_VEC_LOAD_LUMPED], ...
%!                                load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   grp_id_displacement = find([[mesh.groups.tria6].id] == 4);
%!   elem_id_displacement = mesh.groups.tria6(grp_id_displacement).elements;
%!   elno_id_displacement = mesh.elements.tria6(elem_id_displacement, :);
%!   delta = mean(sol_stat.def(elno_id_displacement, 3));
%!   grp_id_stress = find([[mesh.groups.tria6].id] == 5);
%!   elem_id_stress = mesh.groups.tria6(grp_id_stress).elements;
%!   elno_id_stress = mesh.elements.tria6(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(mesh.elements.tet10 == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(sol_stat.stress.taum.tet10(ridx(j), cidx(j), :), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   sigma1_max = 0;
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max = max(sigma1_max, max(eig(TAU)));
%!   endfor
%!   fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   ## K.J.Bathe page 329 4.20b
%!   sigma1_max_ref = 0.6056;
%!   delta_ref = -1.669;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   assert(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%!   assert(delta, delta_ref, 0.04 * abs(delta_ref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO 1
%! ## K.J.Bathe 2002, page 328 4.20a
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   animate = true;
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   mesh_size = 1;
%!   p1 = 0.006;
%!   E = 55;
%!   nu = 0.3;
%!   rho = 1000e-12;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%!   fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Circle(3) = {3,11,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Circle(7) = {7,12,8};\n");
%!   fputs(fd, "Line(8) = {8,9};\n");
%!   fputs(fd, "Line(9) = {9,10};\n");
%!   fputs(fd, "Line(10) = {10,1};\n");
%!   fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%!   fputs(fd, "Plane Surface(14) = {11};\n");
%!   fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; };\n", mesh_size);
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%!   fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%!   fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%!   fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.tria6].id] == 2);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, 1:3) = true;
%!   grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!   elem_id_p1 = mesh.groups.tria6(grp_id_p1).elements;
%!   elno_p1 = mesh.elements.tria6(elem_id_p1, :);

%!   load_case.pressure.tria6.elements = elno_p1;
%!   load_case.pressure.tria6.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_VEC_LOAD_LUMPED], ...
%!                                load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   grp_id_displacement = find([[mesh.groups.tria6].id] == 4);
%!   elem_id_displacement = mesh.groups.tria6(grp_id_displacement).elements;
%!   elno_id_displacement = mesh.elements.tria6(elem_id_displacement, :);
%!   delta = mean(sol_stat.def(elno_id_displacement, 3));
%!   grp_id_stress = find([[mesh.groups.tria6].id] == 5);
%!   elem_id_stress = mesh.groups.tria6(grp_id_stress).elements;
%!   elno_id_stress = mesh.elements.tria6(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(mesh.elements.tet10 == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(sol_stat.stress.taum.tet10(ridx(j), cidx(j), :), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   sigma1_max = 0;
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max = max(sigma1_max, max(eig(TAU)));
%!   endfor
%!   fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   ## K.J.Bathe page 329 4.20b
%!   sigma1_max_ref = 0.6056;
%!   delta_ref = -1.669;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   if (animate)
%!     opt_anim.scale_def = 10;
%!     opt_anim.animation_delay = 1;
%!     opt_anim.print_and_exit = true;
%!     opt_anim.print_to_file = filename;
%!     opt_anim.rotation_angle = [90, 0, 0] * pi / 180;
%!     opt_anim.skin_only = true;
%!     opt_anim.show_element = true;
%!     unwind_protect
%!       fem_post_sol_external(mesh, sol_stat, opt_anim);
%!       [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title("Gmsh - deformed mesh / continuous stress tensor");
%!     unwind_protect_cleanup
%!       unlink([opt_anim.print_to_file, "_001.jpg"]);
%!     end_unwind_protect
%!   endif
%!   assert(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%!   assert(delta, delta_ref, 0.04 * abs(delta_ref));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
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
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   for i=1:4
%!     fd = -1;
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "wt");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       switch (i)
%!         case {1,2}
%!           order = 1;
%!         otherwise
%!           order = 2;
%!       endswitch
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "a=%g;\n", a);
%!       fprintf(fd, "b=%g;\n", b);
%!       fprintf(fd, "c=%g;\n", c);
%!       fprintf(fd, "h=%g;\n", h * order);
%!       fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!       fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!       fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!       fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,1};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp[] = Extrude {0.0,0.0,c} {\n");
%!       switch (i)
%!         case 1
%!           fprintf(fd, "  Surface{6}; Layers{%d}; Recombine;\n", ceil(c/h));
%!         otherwise
%!           fputs(fd, "  Surface{6};\n");
%!       endswitch
%!       fputs(fd, "};\n");
%!       switch (i)
%!         case 1
%!           fprintf(fd, "Recombine Surface{6,tmp[0]};\n");
%!       endswitch
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!       fprintf(fd, "Physical Surface(\"right\",%d) = {tmp[2]};\n", 100 * i + 1);
%!       fprintf(fd, "Physical Surface(\"rear\",%d) = {tmp[3]};\n", 100 * i + 2);
%!       fprintf(fd, "Physical Surface(\"left\",%d) = {tmp[4]};\n", 100 * i + 3);
%!       fprintf(fd, "Physical Surface(\"front\",%d) = {tmp[5]};\n", 100 * i + 4);
%!       fprintf(fd, "Physical Surface(\"top\",%d) = {tmp[0]};\n", 100 * i + 5);
%!       fprintf(fd, "Physical Surface(\"bottom\",%d) = {6};\n", 100 * i + 6);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", sprintf("%d", order), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     unlink([filename, ".geo"]);
%!     data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!     unlink([filename, ".msh"]);
%!   endfor
%!   data(2).mesh.nodes(:, 1) += a;
%!   data(3).mesh.nodes(:, 1) += a;
%!   data(3).mesh.nodes(:, 2) += b;
%!   data(4).mesh.nodes(:, 2) += b;
%!   for i=1:numel(data)
%!     if (isfield(data(i).mesh.elements, "iso8"))
%!       data(i).mesh.materials.iso8 = ones(rows(data(i).mesh.elements.iso8), 1, "int32");
%!     else
%!       data(i).mesh.materials.tet10 = ones(rows(data(i).mesh.elements.tet10), 1, "int32");
%!     endif
%!     data(i).mesh.material_data.rho = rho;
%!     data(i).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   endfor
%!   mesh = fem_post_mesh_merge(data);
%!   mesh.nodes(:, 3) -= 0.5 * c;
%!   constr = FEM_CT_FIXED;
%!   mesh.elements.sfncon4(1).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 103)).elements, :);
%!   mesh.elements.sfncon4(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 401)).nodes(:);
%!   mesh.elements.sfncon4(1).maxdist = maxdist;
%!   mesh.elements.sfncon4(1).constraint = constr;
%!   mesh.elements.sfncon4(2).master = mesh.elements.iso4(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 102)).elements, :);
%!   mesh.elements.sfncon4(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 204)).nodes(:);
%!   mesh.elements.sfncon4(2).maxdist = maxdist;
%!   mesh.elements.sfncon4(2).constraint = constr;
%!   mesh.elements.sfncon6(1).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 402)).elements, :);
%!   mesh.elements.sfncon6(1).slave = mesh.groups.tria6(find([[mesh.groups.tria6].id] == 304)).nodes(:);
%!   mesh.elements.sfncon6(1).maxdist = maxdist;
%!   mesh.elements.sfncon6(1).constraint = constr;
%!   mesh.elements.sfncon6(2).master = mesh.elements.tria6(mesh.groups.tria6(find([[mesh.groups.tria6].id] == 301)).elements, :);
%!   mesh.elements.sfncon6(2).slave = mesh.groups.iso4(find([[mesh.groups.iso4].id] == 203)).nodes(:);
%!   mesh.elements.sfncon6(2).maxdist = maxdist;
%!   mesh.elements.sfncon6(2).constraint = constr;
%!   slave_nodes = [];
%!   for i=1:numel(mesh.elements.sfncon4)
%!     idx_slave = [];
%!     for j=1:numel(slave_nodes)
%!       idx_slave = [idx_slave; find(mesh.elements.sfncon4(i).slave == slave_nodes(j))];
%!     endfor
%!     if (numel(idx_slave))
%!       mesh.elements.sfncon4(i).slave(idx_slave) = 0;
%!       mesh.elements.sfncon4(i).slave = mesh.elements.sfncon4(i).slave(find(mesh.elements.sfncon4(i).slave));
%!     endif
%!     slave_nodes = [slave_nodes; mesh.elements.sfncon4(i).slave];
%!   endfor

%!   for i=1:numel(mesh.elements.sfncon6)
%!     idx_slave = [];
%!     for j=1:numel(slave_nodes)
%!       idx_slave = [idx_slave; find(mesh.elements.sfncon6(i).slave == slave_nodes(j))];
%!     endfor
%!     if (numel(idx_slave))
%!       mesh.elements.sfncon6(i).slave(idx_slave) = 0;
%!       mesh.elements.sfncon6(i).slave = mesh.elements.sfncon6(i).slave(find(mesh.elements.sfncon6(i).slave));
%!     endif
%!     slave_nodes = [slave_nodes; mesh.elements.sfncon6(i).slave];
%!   endfor
%!   group_id4 = [[mesh.groups.iso4].id];
%!   group_id6 = [[mesh.groups.tria6].id];

%!   node_front = [[mesh.groups.iso4(find(group_id4 == 104))].nodes, [mesh.groups.tria6(find(group_id6 == 404))].nodes];
%!   node_rear = [[mesh.groups.iso4(find(group_id4 == 202))].nodes, [mesh.groups.tria6(find(group_id6 == 302))].nodes];
%!   node_right = [mesh.groups.iso4(find((group_id4 == 101) | (group_id4 == 201))).nodes, ...
%!                 mesh.groups.tria6(find((group_id6==303)|(group_id6==403))).nodes];
%!   node_bottom = [[mesh.groups.iso4(find(mod(group_id4, 100) == 6))].nodes, ...
%!                  [mesh.groups.tria6(find(mod(group_id6,100) == 6))].nodes];
%!   node_top = [[mesh.groups.iso4(find(mod(group_id4, 100) == 5))].nodes, ...
%!               [mesh.groups.tria6(find(mod(group_id6,100) == 5))].nodes];
%!   node_constr = [node_front, node_rear, node_right, node_bottom, node_top];
%!   mesh.elements.joints = repmat(struct("nodes",[],"C",[]), 1, numel(node_constr));
%!   load_case.joints = repmat(struct("U",[]), 1, numel(node_constr));
%!   idx_joint = int32(0);
%!   for i=1:numel(node_front)
%!     if (~numel(find(slave_nodes == node_front(i))))
%!       mesh.elements.joints(++idx_joint).C = [[1, 0, 0; 0, 0, 1],  zeros(2, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_front(i);
%!       load_case.joints(idx_joint).U = [gamma * mesh.nodes(node_front(i), 3); 0];
%!     endif
%!   endfor
%!   for i=1:numel(node_rear)
%!     if (~numel(find(slave_nodes == node_rear(i))))
%!       mesh.elements.joints(++idx_joint).C = [[1, 0, 0; 0, 0, 1], zeros(2, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_rear(i);
%!       load_case.joints(idx_joint).U = [gamma * mesh.nodes(node_rear(i), 3); 0];
%!     endif
%!   endfor

%!   for i=1:numel(node_right)
%!     if (~numel(find(slave_nodes == node_right(i))))
%!       mesh.elements.joints(++idx_joint).C = [0,1,0, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_right(i);
%!       load_case.joints(idx_joint).U = 0;
%!     endif
%!   endfor

%!   for i=1:numel(node_bottom)
%!     xi = mesh.nodes(node_bottom(i), 1);
%!     if (~(numel(find(slave_nodes == node_bottom(i))) || xi == 0 || xi == 2 * a))
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_bottom(i);
%!       load_case.joints(idx_joint).U = 0;
%!     endif
%!   endfor

%!   for i=1:numel(node_top)
%!     xi = mesh.nodes(node_top(i), 1);
%!     if (~(numel(find(slave_nodes == node_top(i))) || xi == 0 || xi == 2 * a))
%!       mesh.elements.joints(++idx_joint).C = [0,0,1, zeros(1, 3)];
%!       mesh.elements.joints(idx_joint).nodes = node_top(i);
%!       load_case.joints(idx_joint).U = 0;
%!     endif
%!   endfor

%!   mesh.elements.joints = mesh.elements.joints(1:idx_joint);
%!   load_case.joints = load_case.joints(1:idx_joint);
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   if (eliminate)
%!     K_mat_type = FEM_MAT_STIFFNESS;
%!   else
%!     K_mat_type = FEM_MAT_STIFFNESS;
%!   endif
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info] = fem_ass_matrix(mesh, ...
%!                                       dof_map, ...
%!                                       [K_mat_type, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case);
%!   if (eliminate)
%!     [Tred, Kred, Mred, Rred] = fem_cms_constr_elim(mesh, dof_map, mat_ass);
%!     opt_ls.refine_max_iter = int32(100);
%!     Kfact = fem_sol_factor(Kred, opt_ls);
%!     Ured = Kfact \ Rred;
%!     sol_stat.def = fem_post_def_nodal(mesh, dof_map, Tred * Ured);
%!   else
%!     sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   endif

%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   if (animate)
%!     opt_anim.scale_def = 1;
%!     opt_anim.animation_delay = 1;
%!     opt_anim.print_and_exit = true;
%!     opt_anim.print_to_file = filename;
%!     opt_anim.rotation_angle = [286, 2, 205] * pi / 180;
%!     opt_anim.skin_only = true;
%!     opt_anim.show_element = true;
%!     unwind_protect
%!       fem_post_sol_external(mesh, sol_stat, opt_anim);
%!       [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title("Gmsh - deformed mesh / continuous stress tensor");
%!     unwind_protect_cleanup
%!       unlink([opt_anim.print_to_file, "_001.jpg"]);
%!     end_unwind_protect
%!   endif
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST17
%! fd = -1;
%! filename = "";
%! unwind_protect
%!   unwind_protect
%!     fn = fullfile(tempdir(), "fem_pre_mesh_import_XXXXXX");
%!     if (ispc())
%!       fn(fn == "\\") = "/";
%!     endif
%!     [fd, filename] = mkstemp(fn);
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", fn);
%!     endif
%!     fputs(fd, "$ Input file from EOSSP\n");
%!     fputs(fd, "XFEP 'HEXAHEDRAL MESH'\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ Nodal coordinates\n");
%!     fputs(fd, "KOOR KNR     1 X      0.00  Y      0.00  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     2 X      0.00  Y    -15.00  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     3 X      8.59  Y    -13.20  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     4 X     13.71  Y     -6.08  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     5 X     15.75  Y      0.12  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     6 X     13.23  Y      7.08  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     7 X     12.08  Y     10.11  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     8 X      0.00  Y     10.11  Z      0.00\n");
%!     fputs(fd, "KOOR KNR     9 X      0.00  Y    -30.00  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    10 X     16.36  Y    -25.15  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    11 X     27.42  Y    -12.16  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    12 X     30.00  Y      0.45  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    13 X     26.45  Y     14.15  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    14 X     23.01  Y     19.25  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    15 X      0.00  Y     19.25  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    22 X     28.78  Y      8.46  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    39 X      0.00  Y     30.00  Z      0.00\n");
%!     fputs(fd, "KOOR KNR    50 X     15.66  Y     25.59  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   103 X     -8.59  Y    -13.20  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   104 X    -13.71  Y     -6.08  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   105 X    -15.75  Y      0.12  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   106 X    -13.23  Y      7.08  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   107 X    -12.08  Y     10.11  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   110 X    -16.36  Y    -25.15  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   111 X    -27.42  Y    -12.16  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   112 X    -30.00  Y      0.45  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   113 X    -26.45  Y     14.15  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   114 X    -23.01  Y     19.25  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   122 X    -28.78  Y      8.46  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   150 X    -15.66  Y     25.59  Z      0.00\n");
%!     fputs(fd, "KOOR KNR   201 X      0.00  Y      0.00  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   202 X      0.00  Y    -15.00  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   203 X      8.59  Y    -13.20  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   204 X     13.71  Y     -6.08  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   205 X     15.75  Y      0.12  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   206 X     13.23  Y      7.08  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   207 X     12.08  Y     10.11  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   208 X      0.00  Y     10.11  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   209 X      0.00  Y    -30.00  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   210 X     16.36  Y    -25.15  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   211 X     27.42  Y    -12.16  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   212 X     30.00  Y      0.45  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   213 X     26.45  Y     14.15  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   214 X     23.01  Y     19.25  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   215 X      0.00  Y     19.25  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   222 X     28.78  Y      8.46  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   239 X      0.00  Y     30.00  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   250 X     15.66  Y     25.59  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   303 X     -8.59  Y    -13.20  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   304 X    -13.71  Y     -6.08  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   305 X    -15.75  Y      0.12  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   306 X    -13.23  Y      7.08  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   307 X    -12.08  Y     10.11  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   310 X    -16.36  Y    -25.15  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   311 X    -27.42  Y    -12.16  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   312 X    -30.00  Y      0.45  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   313 X    -26.45  Y     14.15  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   314 X    -23.01  Y     19.25  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   322 X    -28.78  Y      8.46  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   350 X    -15.66  Y     25.59  Z     -4.17\n");
%!     fputs(fd, "KOOR KNR   401 X      0.00  Y      0.00  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   402 X      0.00  Y    -15.00  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   403 X      8.59  Y    -13.20  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   404 X     13.71  Y     -6.08  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   405 X     15.75  Y      0.12  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   406 X     13.23  Y      7.08  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   407 X     12.08  Y     10.11  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   408 X      0.00  Y     10.11  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   409 X      0.00  Y    -30.00  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   410 X     16.36  Y    -25.15  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   411 X     27.42  Y    -12.16  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   412 X     30.00  Y      0.45  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   413 X     26.45  Y     14.15  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   414 X     23.01  Y     19.25  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   415 X      0.00  Y     19.25  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   422 X     28.78  Y      8.46  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   439 X      0.00  Y     30.00  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   450 X     15.66  Y     25.59  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   503 X     -8.59  Y    -13.20  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   504 X    -13.71  Y     -6.08  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   505 X    -15.75  Y      0.12  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   506 X    -13.23  Y      7.08  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   507 X    -12.08  Y     10.11  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   510 X    -16.36  Y    -25.15  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   511 X    -27.42  Y    -12.16  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   512 X    -30.00  Y      0.45  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   513 X    -26.45  Y     14.15  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   514 X    -23.01  Y     19.25  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   522 X    -28.78  Y      8.46  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   550 X    -15.66  Y     25.59  Z     -8.33\n");
%!     fputs(fd, "KOOR KNR   601 X      0.00  Y      0.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   602 X      0.00  Y    -15.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   603 X      8.59  Y    -13.20  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   604 X     13.71  Y     -6.08  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   605 X     15.75  Y      0.12  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   606 X     13.23  Y      7.08  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   607 X     12.08  Y     10.11  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   608 X      0.00  Y     10.11  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   609 X      0.00  Y    -30.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   610 X     16.36  Y    -25.15  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   611 X     27.42  Y    -12.16  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   612 X     30.00  Y      0.45  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   613 X     26.45  Y     14.15  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   614 X     23.01  Y     19.25  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   615 X      0.00  Y     19.25  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   616 X      0.00  Y    -35.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   617 X     19.08  Y    -29.34  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   618 X     32.00  Y    -14.19  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   619 X     35.34  Y      0.63  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   620 X     33.84  Y      9.95  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   621 X     32.26  Y     19.25  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   622 X     28.78  Y      8.46  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   625 X      0.00  Y     41.75  Z    -14.81\n");
%!     fputs(fd, "KOOR KNR   626 X      0.00  Y     45.50  Z    -17.12\n");
%!     fputs(fd, "KOOR KNR   627 X      2.41  Y     44.87  Z    -16.73\n");
%!     fputs(fd, "KOOR KNR   628 X      3.37  Y     43.40  Z    -15.83\n");
%!     fputs(fd, "KOOR KNR   629 X      3.86  Y     42.54  Z    -15.30\n");
%!     fputs(fd, "KOOR KNR   630 X      3.69  Y     41.08  Z    -14.40\n");
%!     fputs(fd, "KOOR KNR   631 X      2.66  Y     38.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   632 X      0.00  Y     38.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   633 X      0.00  Y     64.25  Z    -28.67\n");
%!     fputs(fd, "KOOR KNR   634 X     13.76  Y     59.55  Z    -25.83\n");
%!     fputs(fd, "KOOR KNR   635 X     20.20  Y     51.66  Z    -21.04\n");
%!     fputs(fd, "KOOR KNR   636 X     22.04  Y     46.26  Z    -17.74\n");
%!     fputs(fd, "KOOR KNR   637 X     22.14  Y     37.75  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   638 X     19.19  Y     30.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   639 X      0.00  Y     30.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   640 X      0.00  Y     69.25  Z    -31.75\n");
%!     fputs(fd, "KOOR KNR   641 X     16.82  Y     63.51  Z    -28.30\n");
%!     fputs(fd, "KOOR KNR   642 X     24.69  Y     53.87  Z    -22.46\n");
%!     fputs(fd, "KOOR KNR   643 X     26.94  Y     47.26  Z    -18.43\n");
%!     fputs(fd, "KOOR KNR   644 X     28.87  Y     37.58  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   645 X     30.32  Y     30.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   650 X     15.66  Y     25.59  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   703 X     -8.59  Y    -13.20  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   704 X    -13.71  Y     -6.08  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   705 X    -15.75  Y      0.12  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   706 X    -13.23  Y      7.08  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   707 X    -12.08  Y     10.11  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   710 X    -16.36  Y    -25.15  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   711 X    -27.42  Y    -12.16  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   712 X    -30.00  Y      0.45  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   713 X    -26.45  Y     14.15  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   714 X    -23.01  Y     19.25  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   717 X    -19.08  Y    -29.34  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   718 X    -32.00  Y    -14.19  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   719 X    -35.34  Y      0.63  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   720 X    -33.84  Y      9.95  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   721 X    -32.26  Y     19.25  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   722 X    -28.78  Y      8.46  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   727 X     -2.41  Y     44.87  Z    -16.73\n");
%!     fputs(fd, "KOOR KNR   728 X     -3.37  Y     43.40  Z    -15.83\n");
%!     fputs(fd, "KOOR KNR   729 X     -3.86  Y     42.54  Z    -15.30\n");
%!     fputs(fd, "KOOR KNR   730 X     -3.69  Y     41.08  Z    -14.40\n");
%!     fputs(fd, "KOOR KNR   731 X     -2.66  Y     38.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   734 X    -13.76  Y     59.55  Z    -25.83\n");
%!     fputs(fd, "KOOR KNR   735 X    -20.20  Y     51.66  Z    -21.04\n");
%!     fputs(fd, "KOOR KNR   736 X    -22.04  Y     46.26  Z    -17.74\n");
%!     fputs(fd, "KOOR KNR   737 X    -22.14  Y     37.75  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   738 X    -19.19  Y     30.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   741 X    -16.82  Y     63.51  Z    -28.30\n");
%!     fputs(fd, "KOOR KNR   742 X    -24.69  Y     53.87  Z    -22.46\n");
%!     fputs(fd, "KOOR KNR   743 X    -26.94  Y     47.26  Z    -18.43\n");
%!     fputs(fd, "KOOR KNR   744 X    -28.87  Y     37.58  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   745 X    -30.32  Y     30.00  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   750 X    -15.66  Y     25.59  Z    -12.50\n");
%!     fputs(fd, "KOOR KNR   801 X      0.00  Y      0.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   802 X      0.00  Y    -15.00  Z    -17.04\n");
%!     fputs(fd, "KOOR KNR   803 X      8.59  Y    -13.20  Z    -17.37\n");
%!     fputs(fd, "KOOR KNR   804 X     13.71  Y     -6.08  Z    -18.68\n");
%!     fputs(fd, "KOOR KNR   805 X     15.75  Y      0.12  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   806 X     13.23  Y      7.08  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   807 X     12.08  Y     10.11  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   808 X      0.00  Y     10.11  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   809 X      0.00  Y    -30.00  Z    -14.25\n");
%!     fputs(fd, "KOOR KNR   810 X     16.36  Y    -25.15  Z    -15.13\n");
%!     fputs(fd, "KOOR KNR   811 X     27.42  Y    -12.16  Z    -17.50\n");
%!     fputs(fd, "KOOR KNR   812 X     30.00  Y      0.45  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   813 X     26.45  Y     14.15  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   814 X     23.01  Y     19.25  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   815 X      0.00  Y     19.25  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   816 X      0.00  Y    -35.00  Z    -13.32\n");
%!     fputs(fd, "KOOR KNR   817 X     19.08  Y    -29.34  Z    -14.34\n");
%!     fputs(fd, "KOOR KNR   818 X     32.00  Y    -14.19  Z    -17.10\n");
%!     fputs(fd, "KOOR KNR   819 X     35.34  Y      0.63  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   820 X     33.84  Y      9.95  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   821 X     32.26  Y     19.25  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   822 X     28.78  Y      8.46  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   825 X      0.00  Y     41.75  Z    -21.37\n");
%!     fputs(fd, "KOOR KNR   826 X      0.00  Y     45.50  Z    -22.91\n");
%!     fputs(fd, "KOOR KNR   827 X      2.41  Y     44.87  Z    -22.65\n");
%!     fputs(fd, "KOOR KNR   828 X      3.37  Y     43.40  Z    -22.05\n");
%!     fputs(fd, "KOOR KNR   829 X      3.86  Y     42.54  Z    -21.70\n");
%!     fputs(fd, "KOOR KNR   830 X      3.69  Y     41.08  Z    -21.10\n");
%!     fputs(fd, "KOOR KNR   831 X      2.66  Y     38.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   832 X      0.00  Y     38.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   833 X      0.00  Y     64.25  Z    -30.61\n");
%!     fputs(fd, "KOOR KNR   834 X     13.76  Y     59.55  Z    -28.72\n");
%!     fputs(fd, "KOOR KNR   835 X     20.20  Y     51.66  Z    -25.53\n");
%!     fputs(fd, "KOOR KNR   836 X     22.04  Y     46.26  Z    -23.32\n");
%!     fputs(fd, "KOOR KNR   837 X     22.14  Y     37.75  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   838 X     19.19  Y     30.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   839 X      0.00  Y     30.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   840 X      0.00  Y     69.25  Z    -32.67\n");
%!     fputs(fd, "KOOR KNR   841 X     16.82  Y     63.51  Z    -30.36\n");
%!     fputs(fd, "KOOR KNR   842 X     24.69  Y     53.87  Z    -26.47\n");
%!     fputs(fd, "KOOR KNR   843 X     26.94  Y     47.26  Z    -23.78\n");
%!     fputs(fd, "KOOR KNR   844 X     28.87  Y     37.58  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   845 X     30.32  Y     30.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   850 X     15.66  Y     25.59  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   903 X     -8.59  Y    -13.20  Z    -17.37\n");
%!     fputs(fd, "KOOR KNR   904 X    -13.71  Y     -6.08  Z    -18.68\n");
%!     fputs(fd, "KOOR KNR   905 X    -15.75  Y      0.12  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   906 X    -13.23  Y      7.08  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   907 X    -12.08  Y     10.11  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   910 X    -16.36  Y    -25.15  Z    -15.13\n");
%!     fputs(fd, "KOOR KNR   911 X    -27.42  Y    -12.16  Z    -17.50\n");
%!     fputs(fd, "KOOR KNR   912 X    -30.00  Y      0.45  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   913 X    -26.45  Y     14.15  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   914 X    -23.01  Y     19.25  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   917 X    -19.08  Y    -29.34  Z    -14.34\n");
%!     fputs(fd, "KOOR KNR   918 X    -32.00  Y    -14.19  Z    -17.10\n");
%!     fputs(fd, "KOOR KNR   919 X    -35.34  Y      0.63  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   920 X    -33.84  Y      9.95  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   921 X    -32.26  Y     19.25  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   922 X    -28.78  Y      8.46  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   927 X     -2.41  Y     44.87  Z    -22.65\n");
%!     fputs(fd, "KOOR KNR   928 X     -3.37  Y     43.40  Z    -22.05\n");
%!     fputs(fd, "KOOR KNR   929 X     -3.86  Y     42.54  Z    -21.70\n");
%!     fputs(fd, "KOOR KNR   930 X     -3.69  Y     41.08  Z    -21.10\n");
%!     fputs(fd, "KOOR KNR   931 X     -2.66  Y     38.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   934 X    -13.76  Y     59.55  Z    -28.72\n");
%!     fputs(fd, "KOOR KNR   935 X    -20.20  Y     51.66  Z    -25.53\n");
%!     fputs(fd, "KOOR KNR   936 X    -22.04  Y     46.26  Z    -23.32\n");
%!     fputs(fd, "KOOR KNR   937 X    -22.14  Y     37.75  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   938 X    -19.19  Y     30.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   941 X    -16.82  Y     63.51  Z    -30.36\n");
%!     fputs(fd, "KOOR KNR   942 X    -24.69  Y     53.87  Z    -26.47\n");
%!     fputs(fd, "KOOR KNR   943 X    -26.94  Y     47.26  Z    -23.78\n");
%!     fputs(fd, "KOOR KNR   944 X    -28.87  Y     37.58  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   945 X    -30.32  Y     30.00  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR   950 X    -15.66  Y     25.59  Z    -19.83\n");
%!     fputs(fd, "KOOR KNR  1001 X      0.00  Y      0.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1002 X      0.00  Y    -15.00  Z    -21.58\n");
%!     fputs(fd, "KOOR KNR  1003 X      8.59  Y    -13.20  Z    -22.24\n");
%!     fputs(fd, "KOOR KNR  1004 X     13.71  Y     -6.08  Z    -24.87\n");
%!     fputs(fd, "KOOR KNR  1005 X     15.75  Y      0.12  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1006 X     13.23  Y      7.08  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1007 X     12.08  Y     10.11  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1008 X      0.00  Y     10.11  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1009 X      0.00  Y    -30.00  Z    -16.00\n");
%!     fputs(fd, "KOOR KNR  1010 X     16.36  Y    -25.15  Z    -17.75\n");
%!     fputs(fd, "KOOR KNR  1011 X     27.42  Y    -12.16  Z    -22.50\n");
%!     fputs(fd, "KOOR KNR  1012 X     30.00  Y      0.45  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1013 X     26.45  Y     14.15  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1014 X     23.01  Y     19.25  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1015 X      0.00  Y     19.25  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1016 X      0.00  Y    -35.00  Z    -14.13\n");
%!     fputs(fd, "KOOR KNR  1017 X     19.08  Y    -29.34  Z    -16.18\n");
%!     fputs(fd, "KOOR KNR  1018 X     32.00  Y    -14.19  Z    -21.70\n");
%!     fputs(fd, "KOOR KNR  1019 X     35.34  Y      0.63  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1020 X     33.84  Y      9.95  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1021 X     32.26  Y     19.25  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1022 X     28.78  Y      8.46  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1025 X      0.00  Y     41.75  Z    -27.94\n");
%!     fputs(fd, "KOOR KNR  1026 X      0.00  Y     45.50  Z    -28.71\n");
%!     fputs(fd, "KOOR KNR  1027 X      2.41  Y     44.87  Z    -28.58\n");
%!     fputs(fd, "KOOR KNR  1028 X      3.37  Y     43.40  Z    -28.28\n");
%!     fputs(fd, "KOOR KNR  1029 X      3.86  Y     42.54  Z    -28.10\n");
%!     fputs(fd, "KOOR KNR  1030 X      3.69  Y     41.08  Z    -27.80\n");
%!     fputs(fd, "KOOR KNR  1031 X      2.66  Y     38.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1032 X      0.00  Y     38.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1033 X      0.00  Y     64.25  Z    -32.56\n");
%!     fputs(fd, "KOOR KNR  1034 X     13.76  Y     59.55  Z    -31.61\n");
%!     fputs(fd, "KOOR KNR  1035 X     20.20  Y     51.66  Z    -30.01\n");
%!     fputs(fd, "KOOR KNR  1036 X     22.04  Y     46.26  Z    -28.91\n");
%!     fputs(fd, "KOOR KNR  1037 X     22.14  Y     37.75  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1038 X     19.19  Y     30.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1039 X      0.00  Y     30.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1040 X      0.00  Y     69.25  Z    -33.58\n");
%!     fputs(fd, "KOOR KNR  1041 X     16.82  Y     63.51  Z    -32.43\n");
%!     fputs(fd, "KOOR KNR  1042 X     24.69  Y     53.87  Z    -30.49\n");
%!     fputs(fd, "KOOR KNR  1043 X     26.94  Y     47.26  Z    -29.14\n");
%!     fputs(fd, "KOOR KNR  1044 X     28.87  Y     37.58  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1045 X     30.32  Y     30.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1050 X     15.66  Y     25.59  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1103 X     -8.59  Y    -13.20  Z    -22.24\n");
%!     fputs(fd, "KOOR KNR  1104 X    -13.71  Y     -6.08  Z    -24.87\n");
%!     fputs(fd, "KOOR KNR  1105 X    -15.75  Y      0.12  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1106 X    -13.23  Y      7.08  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1107 X    -12.08  Y     10.11  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1110 X    -16.36  Y    -25.15  Z    -17.75\n");
%!     fputs(fd, "KOOR KNR  1111 X    -27.42  Y    -12.16  Z    -22.50\n");
%!     fputs(fd, "KOOR KNR  1112 X    -30.00  Y      0.45  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1113 X    -26.45  Y     14.15  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1114 X    -23.01  Y     19.25  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1117 X    -19.08  Y    -29.34  Z    -16.18\n");
%!     fputs(fd, "KOOR KNR  1118 X    -32.00  Y    -14.19  Z    -21.70\n");
%!     fputs(fd, "KOOR KNR  1119 X    -35.34  Y      0.63  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1120 X    -33.84  Y      9.95  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1121 X    -32.26  Y     19.25  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1122 X    -28.78  Y      8.46  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1127 X     -2.41  Y     44.87  Z    -28.58\n");
%!     fputs(fd, "KOOR KNR  1128 X     -3.37  Y     43.40  Z    -28.28\n");
%!     fputs(fd, "KOOR KNR  1129 X     -3.86  Y     42.54  Z    -28.10\n");
%!     fputs(fd, "KOOR KNR  1130 X     -3.69  Y     41.08  Z    -27.80\n");
%!     fputs(fd, "KOOR KNR  1131 X     -2.66  Y     38.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1134 X    -13.76  Y     59.55  Z    -31.61\n");
%!     fputs(fd, "KOOR KNR  1135 X    -20.20  Y     51.66  Z    -30.01\n");
%!     fputs(fd, "KOOR KNR  1136 X    -22.04  Y     46.26  Z    -28.91\n");
%!     fputs(fd, "KOOR KNR  1137 X    -22.14  Y     37.75  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1138 X    -19.19  Y     30.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1141 X    -16.82  Y     63.51  Z    -32.43\n");
%!     fputs(fd, "KOOR KNR  1142 X    -24.69  Y     53.87  Z    -30.49\n");
%!     fputs(fd, "KOOR KNR  1143 X    -26.94  Y     47.26  Z    -29.14\n");
%!     fputs(fd, "KOOR KNR  1144 X    -28.87  Y     37.58  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1145 X    -30.32  Y     30.00  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1150 X    -15.66  Y     25.59  Z    -27.17\n");
%!     fputs(fd, "KOOR KNR  1201 X      0.00  Y      0.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1202 X      0.00  Y    -15.00  Z    -26.12\n");
%!     fputs(fd, "KOOR KNR  1203 X      8.59  Y    -13.20  Z    -27.11\n");
%!     fputs(fd, "KOOR KNR  1204 X     13.71  Y     -6.08  Z    -31.05\n");
%!     fputs(fd, "KOOR KNR  1205 X     15.75  Y      0.12  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1206 X     13.23  Y      7.08  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1207 X     12.08  Y     10.11  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1208 X      0.00  Y     10.11  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1209 X      0.00  Y    -30.00  Z    -17.74\n");
%!     fputs(fd, "KOOR KNR  1210 X     16.36  Y    -25.15  Z    -20.38\n");
%!     fputs(fd, "KOOR KNR  1211 X     27.42  Y    -12.16  Z    -27.50\n");
%!     fputs(fd, "KOOR KNR  1212 X     30.00  Y      0.45  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1213 X     26.45  Y     14.15  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1214 X     23.01  Y     19.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1215 X      0.00  Y     19.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1216 X      0.00  Y    -35.00  Z    -14.95\n");
%!     fputs(fd, "KOOR KNR  1217 X     19.08  Y    -29.34  Z    -18.01\n");
%!     fputs(fd, "KOOR KNR  1218 X     32.00  Y    -14.19  Z    -26.29\n");
%!     fputs(fd, "KOOR KNR  1219 X     35.34  Y      0.63  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1220 X     33.84  Y      9.95  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1221 X     32.26  Y     19.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1222 X     28.78  Y      8.46  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1225 X      0.00  Y     41.75  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1226 X      0.00  Y     45.50  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1227 X      2.41  Y     44.87  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1228 X      3.37  Y     43.40  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1229 X      3.86  Y     42.54  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1230 X      3.69  Y     41.08  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1231 X      2.66  Y     38.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1232 X      0.00  Y     38.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1233 X      0.00  Y     64.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1234 X     13.76  Y     59.55  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1235 X     20.20  Y     51.66  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1236 X     22.04  Y     46.26  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1237 X     22.14  Y     37.75  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1238 X     19.19  Y     30.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1239 X      0.00  Y     30.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1240 X      0.00  Y     69.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1241 X     16.82  Y     63.51  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1242 X     24.69  Y     53.87  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1243 X     26.94  Y     47.26  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1244 X     28.87  Y     37.58  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1245 X     30.32  Y     30.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1250 X     15.66  Y     25.59  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1303 X     -8.59  Y    -13.20  Z    -27.11\n");
%!     fputs(fd, "KOOR KNR  1304 X    -13.71  Y     -6.08  Z    -31.05\n");
%!     fputs(fd, "KOOR KNR  1305 X    -15.75  Y      0.12  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1306 X    -13.23  Y      7.08  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1307 X    -12.08  Y     10.11  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1310 X    -16.36  Y    -25.15  Z    -20.38\n");
%!     fputs(fd, "KOOR KNR  1311 X    -27.42  Y    -12.16  Z    -27.50\n");
%!     fputs(fd, "KOOR KNR  1312 X    -30.00  Y      0.45  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1313 X    -26.45  Y     14.15  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1314 X    -23.01  Y     19.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1317 X    -19.08  Y    -29.34  Z    -18.01\n");
%!     fputs(fd, "KOOR KNR  1318 X    -32.00  Y    -14.19  Z    -26.29\n");
%!     fputs(fd, "KOOR KNR  1319 X    -35.34  Y      0.63  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1320 X    -33.84  Y      9.95  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1321 X    -32.26  Y     19.25  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1322 X    -28.78  Y      8.46  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1327 X     -2.41  Y     44.87  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1328 X     -3.37  Y     43.40  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1329 X     -3.86  Y     42.54  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1330 X     -3.69  Y     41.08  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1331 X     -2.66  Y     38.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1334 X    -13.76  Y     59.55  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1335 X    -20.20  Y     51.66  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1336 X    -22.04  Y     46.26  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1337 X    -22.14  Y     37.75  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1338 X    -19.19  Y     30.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1341 X    -16.82  Y     63.51  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1342 X    -24.69  Y     53.87  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1343 X    -26.94  Y     47.26  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1344 X    -28.87  Y     37.58  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1345 X    -30.32  Y     30.00  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1350 X    -15.66  Y     25.59  Z    -34.50\n");
%!     fputs(fd, "KOOR KNR  1415 X      0.00  Y     19.25  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1425 X      0.00  Y     41.75  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1426 X      0.00  Y     45.50  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1427 X      2.41  Y     44.87  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1428 X      3.37  Y     43.40  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1429 X      3.86  Y     42.54  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1430 X      3.69  Y     41.08  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1431 X      2.66  Y     38.00  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1432 X      0.00  Y     38.00  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1433 X      0.00  Y     64.25  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1434 X     13.76  Y     59.55  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1435 X     20.20  Y     51.66  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1436 X     22.04  Y     46.26  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1437 X     22.14  Y     37.75  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1438 X     19.19  Y     30.00  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1439 X      0.00  Y     30.00  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1450 X     15.66  Y     25.59  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1527 X     -2.41  Y     44.87  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1528 X     -3.37  Y     43.40  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1529 X     -3.86  Y     42.54  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1530 X     -3.69  Y     41.08  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1531 X     -2.66  Y     38.00  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1534 X    -13.76  Y     59.55  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1535 X    -20.20  Y     51.66  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1536 X    -22.04  Y     46.26  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1537 X    -22.14  Y     37.75  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1538 X    -19.19  Y     30.00  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1550 X    -15.66  Y     25.59  Z    -38.17\n");
%!     fputs(fd, "KOOR KNR  1615 X      0.00  Y     19.25  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1625 X      0.00  Y     41.75  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1626 X      0.00  Y     45.50  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1627 X      2.41  Y     44.87  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1628 X      3.37  Y     43.40  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1629 X      3.86  Y     42.54  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1630 X      3.69  Y     41.08  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1631 X      2.66  Y     38.00  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1632 X      0.00  Y     38.00  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1633 X      0.00  Y     64.25  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1634 X     13.76  Y     59.55  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1635 X     20.20  Y     51.66  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1636 X     22.04  Y     46.26  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1637 X     22.14  Y     37.75  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1638 X     19.19  Y     30.00  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1639 X      0.00  Y     30.00  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1650 X     15.66  Y     25.59  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1727 X     -2.41  Y     44.87  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1728 X     -3.37  Y     43.40  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1729 X     -3.86  Y     42.54  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1730 X     -3.69  Y     41.08  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1731 X     -2.66  Y     38.00  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1734 X    -13.76  Y     59.55  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1735 X    -20.20  Y     51.66  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1736 X    -22.04  Y     46.26  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1737 X    -22.14  Y     37.75  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1738 X    -19.19  Y     30.00  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1750 X    -15.66  Y     25.59  Z    -41.83\n");
%!     fputs(fd, "KOOR KNR  1815 X      0.00  Y     19.25  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1825 X      0.00  Y     41.75  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1826 X      0.00  Y     45.50  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1827 X      2.41  Y     44.87  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1828 X      3.37  Y     43.40  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1829 X      3.86  Y     42.54  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1830 X      3.69  Y     41.08  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1831 X      2.66  Y     38.00  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1832 X      0.00  Y     38.00  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1833 X      0.00  Y     64.25  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1834 X     13.76  Y     59.55  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1835 X     20.20  Y     51.66  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1836 X     22.04  Y     46.26  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1837 X     22.14  Y     37.75  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1838 X     19.19  Y     30.00  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1839 X      0.00  Y     30.00  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1850 X     15.66  Y     25.59  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1927 X     -2.41  Y     44.87  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1928 X     -3.37  Y     43.40  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1929 X     -3.86  Y     42.54  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1930 X     -3.69  Y     41.08  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1931 X     -2.66  Y     38.00  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1934 X    -13.76  Y     59.55  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1935 X    -20.20  Y     51.66  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1936 X    -22.04  Y     46.26  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1937 X    -22.14  Y     37.75  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1938 X    -19.19  Y     30.00  Z    -45.50\n");
%!     fputs(fd, "KOOR KNR  1950 X    -15.66  Y     25.59  Z    -45.50\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ MATERIALDATEN\n");
%!     fputs(fd, "MATE MNR    1 E  177500.  NUE  0.290  RHO 0.0000711\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ ELEMENTBESCHREIBUNG\n");
%!     fputs(fd, "STRU BNR  1 ENR   1 TNR 12 MNR  1 KNR1     2 KNR2     3 KNR3     4 KNR4     1 $$\n");
%!     fputs(fd, "                                  KNR5   202 KNR6   203 KNR7   204 KNR8   201\n");
%!     fputs(fd, "STRU BNR  1 ENR   2 TNR 11 MNR  1 KNR1     1 KNR2     4 KNR3     5 KNR4   201 $$\n");
%!     fputs(fd, "                                  KNR5   204 KNR6   205\n");
%!     fputs(fd, "STRU BNR  1 ENR   3 TNR 12 MNR  1 KNR1     1 KNR2     6 KNR3     7 KNR4     8 $$\n");
%!     fputs(fd, "                                  KNR5   201 KNR6   206 KNR7   207 KNR8   208\n");
%!     fputs(fd, "STRU BNR  1 ENR   4 TNR 12 MNR  1 KNR1     9 KNR2    10 KNR3     3 KNR4     2 $$\n");
%!     fputs(fd, "                                  KNR5   209 KNR6   210 KNR7   203 KNR8   202\n");
%!     fputs(fd, "STRU BNR  1 ENR   5 TNR 12 MNR  1 KNR1    10 KNR2    11 KNR3     4 KNR4     3 $$\n");
%!     fputs(fd, "                                  KNR5   210 KNR6   211 KNR7   204 KNR8   203\n");
%!     fputs(fd, "STRU BNR  1 ENR   6 TNR 12 MNR  1 KNR1     4 KNR2    11 KNR3    12 KNR4     5 $$\n");
%!     fputs(fd, "                                  KNR5   204 KNR6   211 KNR7   212 KNR8   205\n");
%!     fputs(fd, "STRU BNR  1 ENR   7 TNR 11 MNR  1 KNR1     1 KNR2     5 KNR3     6 KNR4   201 $$\n");
%!     fputs(fd, "                                  KNR5   205 KNR6   206\n");
%!     fputs(fd, "STRU BNR  1 ENR   8 TNR 11 MNR  1 KNR1     5 KNR2    12 KNR3    22 KNR4   205 $$\n");
%!     fputs(fd, "                                  KNR5   212 KNR6   222\n");
%!     fputs(fd, "STRU BNR  1 ENR   9 TNR 12 MNR  1 KNR1     5 KNR2    22 KNR3    13 KNR4     6 $$\n");
%!     fputs(fd, "                                  KNR5   205 KNR6   222 KNR7   213 KNR8   206\n");
%!     fputs(fd, "STRU BNR  1 ENR  10 TNR 12 MNR  1 KNR1     8 KNR2     7 KNR3    14 KNR4    15 $$\n");
%!     fputs(fd, "                                  KNR5   208 KNR6   207 KNR7   214 KNR8   215\n");
%!     fputs(fd, "STRU BNR  1 ENR  11 TNR 12 MNR  1 KNR1     6 KNR2    13 KNR3    14 KNR4     7 $$\n");
%!     fputs(fd, "                                  KNR5   206 KNR6   213 KNR7   214 KNR8   207\n");
%!     fputs(fd, "STRU BNR  1 ENR  41 TNR 11 MNR  1 KNR1    15 KNR2    50 KNR3    39 KNR4   215 $$\n");
%!     fputs(fd, "                                  KNR5   250 KNR6   239\n");
%!     fputs(fd, "STRU BNR  1 ENR  42 TNR 11 MNR  1 KNR1    15 KNR2    14 KNR3    50 KNR4   215 $$\n");
%!     fputs(fd, "                                  KNR5   214 KNR6   250\n");
%!     fputs(fd, "STRU BNR  1 ENR 101 TNR 12 MNR  1 KNR1     2 KNR2   103 KNR3   104 KNR4     1 $$\n");
%!     fputs(fd, "                                  KNR5   202 KNR6   303 KNR7   304 KNR8   201\n");
%!     fputs(fd, "STRU BNR  1 ENR 102 TNR 11 MNR  1 KNR1     1 KNR2   104 KNR3   105 KNR4   201 $$\n");
%!     fputs(fd, "                                  KNR5   304 KNR6   305\n");
%!     fputs(fd, "STRU BNR  1 ENR 103 TNR 12 MNR  1 KNR1     1 KNR2   106 KNR3   107 KNR4     8 $$\n");
%!     fputs(fd, "                                  KNR5   201 KNR6   306 KNR7   307 KNR8   208\n");
%!     fputs(fd, "STRU BNR  1 ENR 104 TNR 12 MNR  1 KNR1     9 KNR2   110 KNR3   103 KNR4     2 $$\n");
%!     fputs(fd, "                                  KNR5   209 KNR6   310 KNR7   303 KNR8   202\n");
%!     fputs(fd, "STRU BNR  1 ENR 105 TNR 12 MNR  1 KNR1   110 KNR2   111 KNR3   104 KNR4   103 $$\n");
%!     fputs(fd, "                                  KNR5   310 KNR6   311 KNR7   304 KNR8   303\n");
%!     fputs(fd, "STRU BNR  1 ENR 106 TNR 12 MNR  1 KNR1   104 KNR2   111 KNR3   112 KNR4   105 $$\n");
%!     fputs(fd, "                                  KNR5   304 KNR6   311 KNR7   312 KNR8   305\n");
%!     fputs(fd, "STRU BNR  1 ENR 107 TNR 11 MNR  1 KNR1     1 KNR2   105 KNR3   106 KNR4   201 $$\n");
%!     fputs(fd, "                                  KNR5   305 KNR6   306\n");
%!     fputs(fd, "STRU BNR  1 ENR 108 TNR 11 MNR  1 KNR1   105 KNR2   112 KNR3   122 KNR4   305 $$\n");
%!     fputs(fd, "                                  KNR5   312 KNR6   322\n");
%!     fputs(fd, "STRU BNR  1 ENR 109 TNR 12 MNR  1 KNR1   105 KNR2   122 KNR3   113 KNR4   106 $$\n");
%!     fputs(fd, "                                  KNR5   305 KNR6   322 KNR7   313 KNR8   306\n");
%!     fputs(fd, "STRU BNR  1 ENR 110 TNR 12 MNR  1 KNR1     8 KNR2   107 KNR3   114 KNR4    15 $$\n");
%!     fputs(fd, "                                  KNR5   208 KNR6   307 KNR7   314 KNR8   215\n");
%!     fputs(fd, "STRU BNR  1 ENR 111 TNR 12 MNR  1 KNR1   106 KNR2   113 KNR3   114 KNR4   107 $$\n");
%!     fputs(fd, "                                  KNR5   306 KNR6   313 KNR7   314 KNR8   307\n");
%!     fputs(fd, "STRU BNR  1 ENR 141 TNR 11 MNR  1 KNR1    15 KNR2   150 KNR3    39 KNR4   215 $$\n");
%!     fputs(fd, "                                  KNR5   350 KNR6   239\n");
%!     fputs(fd, "STRU BNR  1 ENR 142 TNR 11 MNR  1 KNR1    15 KNR2   114 KNR3   150 KNR4   215 $$\n");
%!     fputs(fd, "                                  KNR5   314 KNR6   350\n");
%!     fputs(fd, "STRU BNR  1 ENR 201 TNR 12 MNR  1 KNR1   202 KNR2   203 KNR3   204 KNR4   201 $$\n");
%!     fputs(fd, "                                  KNR5   402 KNR6   403 KNR7   404 KNR8   401\n");
%!     fputs(fd, "STRU BNR  1 ENR 202 TNR 11 MNR  1 KNR1   201 KNR2   204 KNR3   205 KNR4   401 $$\n");
%!     fputs(fd, "                                  KNR5   404 KNR6   405\n");
%!     fputs(fd, "STRU BNR  1 ENR 203 TNR 12 MNR  1 KNR1   201 KNR2   206 KNR3   207 KNR4   208 $$\n");
%!     fputs(fd, "                                  KNR5   401 KNR6   406 KNR7   407 KNR8   408\n");
%!     fputs(fd, "STRU BNR  1 ENR 204 TNR 12 MNR  1 KNR1   209 KNR2   210 KNR3   203 KNR4   202 $$\n");
%!     fputs(fd, "                                  KNR5   409 KNR6   410 KNR7   403 KNR8   402\n");
%!     fputs(fd, "STRU BNR  1 ENR 205 TNR 12 MNR  1 KNR1   210 KNR2   211 KNR3   204 KNR4   203 $$\n");
%!     fputs(fd, "                                  KNR5   410 KNR6   411 KNR7   404 KNR8   403\n");
%!     fputs(fd, "STRU BNR  1 ENR 206 TNR 12 MNR  1 KNR1   204 KNR2   211 KNR3   212 KNR4   205 $$\n");
%!     fputs(fd, "                                  KNR5   404 KNR6   411 KNR7   412 KNR8   405\n");
%!     fputs(fd, "STRU BNR  1 ENR 207 TNR 11 MNR  1 KNR1   201 KNR2   205 KNR3   206 KNR4   401 $$\n");
%!     fputs(fd, "                                  KNR5   405 KNR6   406\n");
%!     fputs(fd, "STRU BNR  1 ENR 208 TNR 11 MNR  1 KNR1   205 KNR2   212 KNR3   222 KNR4   405 $$\n");
%!     fputs(fd, "                                  KNR5   412 KNR6   422\n");
%!     fputs(fd, "STRU BNR  1 ENR 209 TNR 12 MNR  1 KNR1   205 KNR2   222 KNR3   213 KNR4   206 $$\n");
%!     fputs(fd, "                                  KNR5   405 KNR6   422 KNR7   413 KNR8   406\n");
%!     fputs(fd, "STRU BNR  1 ENR 210 TNR 12 MNR  1 KNR1   208 KNR2   207 KNR3   214 KNR4   215 $$\n");
%!     fputs(fd, "                                  KNR5   408 KNR6   407 KNR7   414 KNR8   415\n");
%!     fputs(fd, "STRU BNR  1 ENR 211 TNR 12 MNR  1 KNR1   206 KNR2   213 KNR3   214 KNR4   207 $$\n");
%!     fputs(fd, "                                  KNR5   406 KNR6   413 KNR7   414 KNR8   407\n");
%!     fputs(fd, "STRU BNR  1 ENR 241 TNR 11 MNR  1 KNR1   215 KNR2   250 KNR3   239 KNR4   415 $$\n");
%!     fputs(fd, "                                  KNR5   450 KNR6   439\n");
%!     fputs(fd, "STRU BNR  1 ENR 242 TNR 11 MNR  1 KNR1   215 KNR2   214 KNR3   250 KNR4   415 $$\n");
%!     fputs(fd, "                                  KNR5   414 KNR6   450\n");
%!     fputs(fd, "STRU BNR  1 ENR 301 TNR 12 MNR  1 KNR1   202 KNR2   303 KNR3   304 KNR4   201 $$\n");
%!     fputs(fd, "                                  KNR5   402 KNR6   503 KNR7   504 KNR8   401\n");
%!     fputs(fd, "STRU BNR  1 ENR 302 TNR 11 MNR  1 KNR1   201 KNR2   304 KNR3   305 KNR4   401 $$\n");
%!     fputs(fd, "                                  KNR5   504 KNR6   505\n");
%!     fputs(fd, "STRU BNR  1 ENR 303 TNR 12 MNR  1 KNR1   201 KNR2   306 KNR3   307 KNR4   208 $$\n");
%!     fputs(fd, "                                  KNR5   401 KNR6   506 KNR7   507 KNR8   408\n");
%!     fputs(fd, "STRU BNR  1 ENR 304 TNR 12 MNR  1 KNR1   209 KNR2   310 KNR3   303 KNR4   202 $$\n");
%!     fputs(fd, "                                  KNR5   409 KNR6   510 KNR7   503 KNR8   402\n");
%!     fputs(fd, "STRU BNR  1 ENR 305 TNR 12 MNR  1 KNR1   310 KNR2   311 KNR3   304 KNR4   303 $$\n");
%!     fputs(fd, "                                  KNR5   510 KNR6   511 KNR7   504 KNR8   503\n");
%!     fputs(fd, "STRU BNR  1 ENR 306 TNR 12 MNR  1 KNR1   304 KNR2   311 KNR3   312 KNR4   305 $$\n");
%!     fputs(fd, "                                  KNR5   504 KNR6   511 KNR7   512 KNR8   505\n");
%!     fputs(fd, "STRU BNR  1 ENR 307 TNR 11 MNR  1 KNR1   201 KNR2   305 KNR3   306 KNR4   401 $$\n");
%!     fputs(fd, "                                  KNR5   505 KNR6   506\n");
%!     fputs(fd, "STRU BNR  1 ENR 308 TNR 11 MNR  1 KNR1   305 KNR2   312 KNR3   322 KNR4   505 $$\n");
%!     fputs(fd, "                                  KNR5   512 KNR6   522\n");
%!     fputs(fd, "STRU BNR  1 ENR 309 TNR 12 MNR  1 KNR1   305 KNR2   322 KNR3   313 KNR4   306 $$\n");
%!     fputs(fd, "                                  KNR5   505 KNR6   522 KNR7   513 KNR8   506\n");
%!     fputs(fd, "STRU BNR  1 ENR 310 TNR 12 MNR  1 KNR1   208 KNR2   307 KNR3   314 KNR4   215 $$\n");
%!     fputs(fd, "                                  KNR5   408 KNR6   507 KNR7   514 KNR8   415\n");
%!     fputs(fd, "STRU BNR  1 ENR 311 TNR 12 MNR  1 KNR1   306 KNR2   313 KNR3   314 KNR4   307 $$\n");
%!     fputs(fd, "                                  KNR5   506 KNR6   513 KNR7   514 KNR8   507\n");
%!     fputs(fd, "STRU BNR  1 ENR 341 TNR 11 MNR  1 KNR1   215 KNR2   350 KNR3   239 KNR4   415 $$\n");
%!     fputs(fd, "                                  KNR5   550 KNR6   439\n");
%!     fputs(fd, "STRU BNR  1 ENR 342 TNR 11 MNR  1 KNR1   215 KNR2   314 KNR3   350 KNR4   415 $$\n");
%!     fputs(fd, "                                  KNR5   514 KNR6   550\n");
%!     fputs(fd, "STRU BNR  1 ENR 401 TNR 12 MNR  1 KNR1   402 KNR2   403 KNR3   404 KNR4   401 $$\n");
%!     fputs(fd, "                                  KNR5   602 KNR6   603 KNR7   604 KNR8   601\n");
%!     fputs(fd, "STRU BNR  1 ENR 402 TNR 11 MNR  1 KNR1   401 KNR2   404 KNR3   405 KNR4   601 $$\n");
%!     fputs(fd, "                                  KNR5   604 KNR6   605\n");
%!     fputs(fd, "STRU BNR  1 ENR 403 TNR 12 MNR  1 KNR1   401 KNR2   406 KNR3   407 KNR4   408 $$\n");
%!     fputs(fd, "                                  KNR5   601 KNR6   606 KNR7   607 KNR8   608\n");
%!     fputs(fd, "STRU BNR  1 ENR 404 TNR 12 MNR  1 KNR1   409 KNR2   410 KNR3   403 KNR4   402 $$\n");
%!     fputs(fd, "                                  KNR5   609 KNR6   610 KNR7   603 KNR8   602\n");
%!     fputs(fd, "STRU BNR  1 ENR 405 TNR 12 MNR  1 KNR1   410 KNR2   411 KNR3   404 KNR4   403 $$\n");
%!     fputs(fd, "                                  KNR5   610 KNR6   611 KNR7   604 KNR8   603\n");
%!     fputs(fd, "STRU BNR  1 ENR 406 TNR 12 MNR  1 KNR1   404 KNR2   411 KNR3   412 KNR4   405 $$\n");
%!     fputs(fd, "                                  KNR5   604 KNR6   611 KNR7   612 KNR8   605\n");
%!     fputs(fd, "STRU BNR  1 ENR 407 TNR 11 MNR  1 KNR1   401 KNR2   405 KNR3   406 KNR4   601 $$\n");
%!     fputs(fd, "                                  KNR5   605 KNR6   606\n");
%!     fputs(fd, "STRU BNR  1 ENR 408 TNR 11 MNR  1 KNR1   405 KNR2   412 KNR3   422 KNR4   605 $$\n");
%!     fputs(fd, "                                  KNR5   612 KNR6   622\n");
%!     fputs(fd, "STRU BNR  1 ENR 409 TNR 12 MNR  1 KNR1   405 KNR2   422 KNR3   413 KNR4   406 $$\n");
%!     fputs(fd, "                                  KNR5   605 KNR6   622 KNR7   613 KNR8   606\n");
%!     fputs(fd, "STRU BNR  1 ENR 410 TNR 12 MNR  1 KNR1   408 KNR2   407 KNR3   414 KNR4   415 $$\n");
%!     fputs(fd, "                                  KNR5   608 KNR6   607 KNR7   614 KNR8   615\n");
%!     fputs(fd, "STRU BNR  1 ENR 411 TNR 12 MNR  1 KNR1   406 KNR2   413 KNR3   414 KNR4   407 $$\n");
%!     fputs(fd, "                                  KNR5   606 KNR6   613 KNR7   614 KNR8   607\n");
%!     fputs(fd, "STRU BNR  1 ENR 441 TNR 11 MNR  1 KNR1   415 KNR2   450 KNR3   439 KNR4   615 $$\n");
%!     fputs(fd, "                                  KNR5   650 KNR6   639\n");
%!     fputs(fd, "STRU BNR  1 ENR 442 TNR 11 MNR  1 KNR1   415 KNR2   414 KNR3   450 KNR4   615 $$\n");
%!     fputs(fd, "                                  KNR5   614 KNR6   650\n");
%!     fputs(fd, "STRU BNR  1 ENR 501 TNR 12 MNR  1 KNR1   402 KNR2   503 KNR3   504 KNR4   401 $$\n");
%!     fputs(fd, "                                  KNR5   602 KNR6   703 KNR7   704 KNR8   601\n");
%!     fputs(fd, "STRU BNR  1 ENR 502 TNR 11 MNR  1 KNR1   401 KNR2   504 KNR3   505 KNR4   601 $$\n");
%!     fputs(fd, "                                  KNR5   704 KNR6   705\n");
%!     fputs(fd, "STRU BNR  1 ENR 503 TNR 12 MNR  1 KNR1   401 KNR2   506 KNR3   507 KNR4   408 $$\n");
%!     fputs(fd, "                                  KNR5   601 KNR6   706 KNR7   707 KNR8   608\n");
%!     fputs(fd, "STRU BNR  1 ENR 504 TNR 12 MNR  1 KNR1   409 KNR2   510 KNR3   503 KNR4   402 $$\n");
%!     fputs(fd, "                                  KNR5   609 KNR6   710 KNR7   703 KNR8   602\n");
%!     fputs(fd, "STRU BNR  1 ENR 505 TNR 12 MNR  1 KNR1   510 KNR2   511 KNR3   504 KNR4   503 $$\n");
%!     fputs(fd, "                                  KNR5   710 KNR6   711 KNR7   704 KNR8   703\n");
%!     fputs(fd, "STRU BNR  1 ENR 506 TNR 12 MNR  1 KNR1   504 KNR2   511 KNR3   512 KNR4   505 $$\n");
%!     fputs(fd, "                                  KNR5   704 KNR6   711 KNR7   712 KNR8   705\n");
%!     fputs(fd, "STRU BNR  1 ENR 507 TNR 11 MNR  1 KNR1   401 KNR2   505 KNR3   506 KNR4   601 $$\n");
%!     fputs(fd, "                                  KNR5   705 KNR6   706\n");
%!     fputs(fd, "STRU BNR  1 ENR 508 TNR 11 MNR  1 KNR1   505 KNR2   512 KNR3   522 KNR4   705 $$\n");
%!     fputs(fd, "                                  KNR5   712 KNR6   722\n");
%!     fputs(fd, "STRU BNR  1 ENR 509 TNR 12 MNR  1 KNR1   505 KNR2   522 KNR3   513 KNR4   506 $$\n");
%!     fputs(fd, "                                  KNR5   705 KNR6   722 KNR7   713 KNR8   706\n");
%!     fputs(fd, "STRU BNR  1 ENR 510 TNR 12 MNR  1 KNR1   408 KNR2   507 KNR3   514 KNR4   415 $$\n");
%!     fputs(fd, "                                  KNR5   608 KNR6   707 KNR7   714 KNR8   615\n");
%!     fputs(fd, "STRU BNR  1 ENR 511 TNR 12 MNR  1 KNR1   506 KNR2   513 KNR3   514 KNR4   507 $$\n");
%!     fputs(fd, "                                  KNR5   706 KNR6   713 KNR7   714 KNR8   707\n");
%!     fputs(fd, "STRU BNR  1 ENR 541 TNR 11 MNR  1 KNR1   415 KNR2   550 KNR3   439 KNR4   615 $$\n");
%!     fputs(fd, "                                  KNR5   750 KNR6   639\n");
%!     fputs(fd, "STRU BNR  1 ENR 542 TNR 11 MNR  1 KNR1   415 KNR2   514 KNR3   550 KNR4   615 $$\n");
%!     fputs(fd, "                                  KNR5   714 KNR6   750\n");
%!     fputs(fd, "STRU BNR  2 ENR   1 TNR 12 MNR  1 KNR1   602 KNR2   603 KNR3   604 KNR4   601 $$\n");
%!     fputs(fd, "                                  KNR5   802 KNR6   803 KNR7   804 KNR8   801\n");
%!     fputs(fd, "STRU BNR  2 ENR   2 TNR 11 MNR  1 KNR1   601 KNR2   604 KNR3   605 KNR4   801 $$\n");
%!     fputs(fd, "                                  KNR5   804 KNR6   805\n");
%!     fputs(fd, "STRU BNR  2 ENR   3 TNR 12 MNR  1 KNR1   601 KNR2   606 KNR3   607 KNR4   608 $$\n");
%!     fputs(fd, "                                  KNR5   801 KNR6   806 KNR7   807 KNR8   808\n");
%!     fputs(fd, "STRU BNR  2 ENR   4 TNR 12 MNR  1 KNR1   609 KNR2   610 KNR3   603 KNR4   602 $$\n");
%!     fputs(fd, "                                  KNR5   809 KNR6   810 KNR7   803 KNR8   802\n");
%!     fputs(fd, "STRU BNR  2 ENR   5 TNR 12 MNR  1 KNR1   610 KNR2   611 KNR3   604 KNR4   603 $$\n");
%!     fputs(fd, "                                  KNR5   810 KNR6   811 KNR7   804 KNR8   803\n");
%!     fputs(fd, "STRU BNR  2 ENR   6 TNR 12 MNR  1 KNR1   604 KNR2   611 KNR3   612 KNR4   605 $$\n");
%!     fputs(fd, "                                  KNR5   804 KNR6   811 KNR7   812 KNR8   805\n");
%!     fputs(fd, "STRU BNR  2 ENR   7 TNR 11 MNR  1 KNR1   601 KNR2   605 KNR3   606 KNR4   801 $$\n");
%!     fputs(fd, "                                  KNR5   805 KNR6   806\n");
%!     fputs(fd, "STRU BNR  2 ENR   8 TNR 11 MNR  1 KNR1   605 KNR2   612 KNR3   622 KNR4   805 $$\n");
%!     fputs(fd, "                                  KNR5   812 KNR6   822\n");
%!     fputs(fd, "STRU BNR  2 ENR   9 TNR 12 MNR  1 KNR1   605 KNR2   622 KNR3   613 KNR4   606 $$\n");
%!     fputs(fd, "                                  KNR5   805 KNR6   822 KNR7   813 KNR8   806\n");
%!     fputs(fd, "STRU BNR  2 ENR  10 TNR 12 MNR  1 KNR1   608 KNR2   607 KNR3   614 KNR4   615 $$\n");
%!     fputs(fd, "                                  KNR5   808 KNR6   807 KNR7   814 KNR8   815\n");
%!     fputs(fd, "STRU BNR  2 ENR  11 TNR 12 MNR  1 KNR1   606 KNR2   613 KNR3   614 KNR4   607 $$\n");
%!     fputs(fd, "                                  KNR5   806 KNR6   813 KNR7   814 KNR8   807\n");
%!     fputs(fd, "STRU BNR  2 ENR  12 TNR 12 MNR  1 KNR1   616 KNR2   617 KNR3   610 KNR4   609 $$\n");
%!     fputs(fd, "                                  KNR5   816 KNR6   817 KNR7   810 KNR8   809\n");
%!     fputs(fd, "STRU BNR  2 ENR  13 TNR 12 MNR  1 KNR1   617 KNR2   618 KNR3   611 KNR4   610 $$\n");
%!     fputs(fd, "                                  KNR5   817 KNR6   818 KNR7   811 KNR8   810\n");
%!     fputs(fd, "STRU BNR  2 ENR  14 TNR 12 MNR  1 KNR1   611 KNR2   618 KNR3   619 KNR4   612 $$\n");
%!     fputs(fd, "                                  KNR5   811 KNR6   818 KNR7   819 KNR8   812\n");
%!     fputs(fd, "STRU BNR  2 ENR  15 TNR 12 MNR  1 KNR1   612 KNR2   619 KNR3   620 KNR4   622 $$\n");
%!     fputs(fd, "                                  KNR5   812 KNR6   819 KNR7   820 KNR8   822\n");
%!     fputs(fd, "STRU BNR  2 ENR  16 TNR 12 MNR  1 KNR1   622 KNR2   620 KNR3   621 KNR4   613 $$\n");
%!     fputs(fd, "                                  KNR5   822 KNR6   820 KNR7   821 KNR8   813\n");
%!     fputs(fd, "STRU BNR  2 ENR  17 TNR 11 MNR  1 KNR1   613 KNR2   621 KNR3   614 KNR4   813 $$\n");
%!     fputs(fd, "                                  KNR5   821 KNR6   814\n");
%!     fputs(fd, "STRU BNR  2 ENR  21 TNR 12 MNR  1 KNR1   625 KNR2   628 KNR3   627 KNR4   626 $$\n");
%!     fputs(fd, "                                  KNR5   825 KNR6   828 KNR7   827 KNR8   826\n");
%!     fputs(fd, "STRU BNR  2 ENR  22 TNR 12 MNR  1 KNR1   625 KNR2   630 KNR3   629 KNR4   628 $$\n");
%!     fputs(fd, "                                  KNR5   825 KNR6   830 KNR7   829 KNR8   828\n");
%!     fputs(fd, "STRU BNR  2 ENR  23 TNR 12 MNR  1 KNR1   632 KNR2   631 KNR3   630 KNR4   625 $$\n");
%!     fputs(fd, "                                  KNR5   832 KNR6   831 KNR7   830 KNR8   825\n");
%!     fputs(fd, "STRU BNR  2 ENR  24 TNR 12 MNR  1 KNR1   626 KNR2   627 KNR3   634 KNR4   633 $$\n");
%!     fputs(fd, "                                  KNR5   826 KNR6   827 KNR7   834 KNR8   833\n");
%!     fputs(fd, "STRU BNR  2 ENR  25 TNR 12 MNR  1 KNR1   628 KNR2   635 KNR3   634 KNR4   627 $$\n");
%!     fputs(fd, "                                  KNR5   828 KNR6   835 KNR7   834 KNR8   827\n");
%!     fputs(fd, "STRU BNR  2 ENR  26 TNR 12 MNR  1 KNR1   629 KNR2   636 KNR3   635 KNR4   628 $$\n");
%!     fputs(fd, "                                  KNR5   829 KNR6   836 KNR7   835 KNR8   828\n");
%!     fputs(fd, "STRU BNR  2 ENR  27 TNR 12 MNR  1 KNR1   630 KNR2   637 KNR3   636 KNR4   629 $$\n");
%!     fputs(fd, "                                  KNR5   830 KNR6   837 KNR7   836 KNR8   829\n");
%!     fputs(fd, "STRU BNR  2 ENR  28 TNR 11 MNR  1 KNR1   631 KNR2   637 KNR3   630 KNR4   831 $$\n");
%!     fputs(fd, "                                  KNR5   837 KNR6   830\n");
%!     fputs(fd, "STRU BNR  2 ENR  29 TNR 11 MNR  1 KNR1   638 KNR2   637 KNR3   631 KNR4   838 $$\n");
%!     fputs(fd, "                                  KNR5   837 KNR6   831\n");
%!     fputs(fd, "STRU BNR  2 ENR  30 TNR 12 MNR  1 KNR1   639 KNR2   638 KNR3   631 KNR4   632 $$\n");
%!     fputs(fd, "                                  KNR5   839 KNR6   838 KNR7   831 KNR8   832\n");
%!     fputs(fd, "STRU BNR  2 ENR  31 TNR 12 MNR  1 KNR1   633 KNR2   634 KNR3   641 KNR4   640 $$\n");
%!     fputs(fd, "                                  KNR5   833 KNR6   834 KNR7   841 KNR8   840\n");
%!     fputs(fd, "STRU BNR  2 ENR  32 TNR 12 MNR  1 KNR1   635 KNR2   642 KNR3   641 KNR4   634 $$\n");
%!     fputs(fd, "                                  KNR5   835 KNR6   842 KNR7   841 KNR8   834\n");
%!     fputs(fd, "STRU BNR  2 ENR  33 TNR 12 MNR  1 KNR1   636 KNR2   643 KNR3   642 KNR4   635 $$\n");
%!     fputs(fd, "                                  KNR5   836 KNR6   843 KNR7   842 KNR8   835\n");
%!     fputs(fd, "STRU BNR  2 ENR  34 TNR 12 MNR  1 KNR1   637 KNR2   644 KNR3   643 KNR4   636 $$\n");
%!     fputs(fd, "                                  KNR5   837 KNR6   844 KNR7   843 KNR8   836\n");
%!     fputs(fd, "STRU BNR  2 ENR  35 TNR 12 MNR  1 KNR1   638 KNR2   645 KNR3   644 KNR4   637 $$\n");
%!     fputs(fd, "                                  KNR5   838 KNR6   845 KNR7   844 KNR8   837\n");
%!     fputs(fd, "STRU BNR  2 ENR  41 TNR 11 MNR  1 KNR1   615 KNR2   650 KNR3   639 KNR4   815 $$\n");
%!     fputs(fd, "                                  KNR5   850 KNR6   839\n");
%!     fputs(fd, "STRU BNR  2 ENR  42 TNR 11 MNR  1 KNR1   615 KNR2   614 KNR3   650 KNR4   815 $$\n");
%!     fputs(fd, "                                  KNR5   814 KNR6   850\n");
%!     fputs(fd, "STRU BNR  2 ENR  43 TNR 11 MNR  1 KNR1   650 KNR2   638 KNR3   639 KNR4   850 $$\n");
%!     fputs(fd, "                                  KNR5   838 KNR6   839\n");
%!     fputs(fd, "STRU BNR  2 ENR  44 TNR 11 MNR  1 KNR1   650 KNR2   614 KNR3   638 KNR4   850 $$\n");
%!     fputs(fd, "                                  KNR5   814 KNR6   838\n");
%!     fputs(fd, "STRU BNR  2 ENR  45 TNR 12 MNR  1 KNR1   614 KNR2   621 KNR3   645 KNR4   638 $$\n");
%!     fputs(fd, "                                  KNR5   814 KNR6   821 KNR7   845 KNR8   838\n");
%!     fputs(fd, "STRU BNR  2 ENR 101 TNR 12 MNR  1 KNR1   602 KNR2   703 KNR3   704 KNR4   601 $$\n");
%!     fputs(fd, "                                  KNR5   802 KNR6   903 KNR7   904 KNR8   801\n");
%!     fputs(fd, "STRU BNR  2 ENR 102 TNR 11 MNR  1 KNR1   601 KNR2   704 KNR3   705 KNR4   801 $$\n");
%!     fputs(fd, "                                  KNR5   904 KNR6   905\n");
%!     fputs(fd, "STRU BNR  2 ENR 103 TNR 12 MNR  1 KNR1   601 KNR2   706 KNR3   707 KNR4   608 $$\n");
%!     fputs(fd, "                                  KNR5   801 KNR6   906 KNR7   907 KNR8   808\n");
%!     fputs(fd, "STRU BNR  2 ENR 104 TNR 12 MNR  1 KNR1   609 KNR2   710 KNR3   703 KNR4   602 $$\n");
%!     fputs(fd, "                                  KNR5   809 KNR6   910 KNR7   903 KNR8   802\n");
%!     fputs(fd, "STRU BNR  2 ENR 105 TNR 12 MNR  1 KNR1   710 KNR2   711 KNR3   704 KNR4   703 $$\n");
%!     fputs(fd, "                                  KNR5   910 KNR6   911 KNR7   904 KNR8   903\n");
%!     fputs(fd, "STRU BNR  2 ENR 106 TNR 12 MNR  1 KNR1   704 KNR2   711 KNR3   712 KNR4   705 $$\n");
%!     fputs(fd, "                                  KNR5   904 KNR6   911 KNR7   912 KNR8   905\n");
%!     fputs(fd, "STRU BNR  2 ENR 107 TNR 11 MNR  1 KNR1   601 KNR2   705 KNR3   706 KNR4   801 $$\n");
%!     fputs(fd, "                                  KNR5   905 KNR6   906\n");
%!     fputs(fd, "STRU BNR  2 ENR 108 TNR 11 MNR  1 KNR1   705 KNR2   712 KNR3   722 KNR4   905 $$\n");
%!     fputs(fd, "                                  KNR5   912 KNR6   922\n");
%!     fputs(fd, "STRU BNR  2 ENR 109 TNR 12 MNR  1 KNR1   705 KNR2   722 KNR3   713 KNR4   706 $$\n");
%!     fputs(fd, "                                  KNR5   905 KNR6   922 KNR7   913 KNR8   906\n");
%!     fputs(fd, "STRU BNR  2 ENR 110 TNR 12 MNR  1 KNR1   608 KNR2   707 KNR3   714 KNR4   615 $$\n");
%!     fputs(fd, "                                  KNR5   808 KNR6   907 KNR7   914 KNR8   815\n");
%!     fputs(fd, "STRU BNR  2 ENR 111 TNR 12 MNR  1 KNR1   706 KNR2   713 KNR3   714 KNR4   707 $$\n");
%!     fputs(fd, "                                  KNR5   906 KNR6   913 KNR7   914 KNR8   907\n");
%!     fputs(fd, "STRU BNR  2 ENR 112 TNR 12 MNR  1 KNR1   616 KNR2   717 KNR3   710 KNR4   609 $$\n");
%!     fputs(fd, "                                  KNR5   816 KNR6   917 KNR7   910 KNR8   809\n");
%!     fputs(fd, "STRU BNR  2 ENR 113 TNR 12 MNR  1 KNR1   717 KNR2   718 KNR3   711 KNR4   710 $$\n");
%!     fputs(fd, "                                  KNR5   917 KNR6   918 KNR7   911 KNR8   910\n");
%!     fputs(fd, "STRU BNR  2 ENR 114 TNR 12 MNR  1 KNR1   711 KNR2   718 KNR3   719 KNR4   712 $$\n");
%!     fputs(fd, "                                  KNR5   911 KNR6   918 KNR7   919 KNR8   912\n");
%!     fputs(fd, "STRU BNR  2 ENR 115 TNR 12 MNR  1 KNR1   712 KNR2   719 KNR3   720 KNR4   722 $$\n");
%!     fputs(fd, "                                  KNR5   912 KNR6   919 KNR7   920 KNR8   922\n");
%!     fputs(fd, "STRU BNR  2 ENR 116 TNR 12 MNR  1 KNR1   722 KNR2   720 KNR3   721 KNR4   713 $$\n");
%!     fputs(fd, "                                  KNR5   922 KNR6   920 KNR7   921 KNR8   913\n");
%!     fputs(fd, "STRU BNR  2 ENR 117 TNR 11 MNR  1 KNR1   713 KNR2   721 KNR3   714 KNR4   913 $$\n");
%!     fputs(fd, "                                  KNR5   921 KNR6   914\n");
%!     fputs(fd, "STRU BNR  2 ENR 121 TNR 12 MNR  1 KNR1   625 KNR2   728 KNR3   727 KNR4   626 $$\n");
%!     fputs(fd, "                                  KNR5   825 KNR6   928 KNR7   927 KNR8   826\n");
%!     fputs(fd, "STRU BNR  2 ENR 122 TNR 12 MNR  1 KNR1   625 KNR2   730 KNR3   729 KNR4   728 $$\n");
%!     fputs(fd, "                                  KNR5   825 KNR6   930 KNR7   929 KNR8   928\n");
%!     fputs(fd, "STRU BNR  2 ENR 123 TNR 12 MNR  1 KNR1   632 KNR2   731 KNR3   730 KNR4   625 $$\n");
%!     fputs(fd, "                                  KNR5   832 KNR6   931 KNR7   930 KNR8   825\n");
%!     fputs(fd, "STRU BNR  2 ENR 124 TNR 12 MNR  1 KNR1   626 KNR2   727 KNR3   734 KNR4   633 $$\n");
%!     fputs(fd, "                                  KNR5   826 KNR6   927 KNR7   934 KNR8   833\n");
%!     fputs(fd, "STRU BNR  2 ENR 125 TNR 12 MNR  1 KNR1   728 KNR2   735 KNR3   734 KNR4   727 $$\n");
%!     fputs(fd, "                                  KNR5   928 KNR6   935 KNR7   934 KNR8   927\n");
%!     fputs(fd, "STRU BNR  2 ENR 126 TNR 12 MNR  1 KNR1   729 KNR2   736 KNR3   735 KNR4   728 $$\n");
%!     fputs(fd, "                                  KNR5   929 KNR6   936 KNR7   935 KNR8   928\n");
%!     fputs(fd, "STRU BNR  2 ENR 127 TNR 12 MNR  1 KNR1   730 KNR2   737 KNR3   736 KNR4   729 $$\n");
%!     fputs(fd, "                                  KNR5   930 KNR6   937 KNR7   936 KNR8   929\n");
%!     fputs(fd, "STRU BNR  2 ENR 128 TNR 11 MNR  1 KNR1   731 KNR2   737 KNR3   730 KNR4   931 $$\n");
%!     fputs(fd, "                                  KNR5   937 KNR6   930\n");
%!     fputs(fd, "STRU BNR  2 ENR 129 TNR 11 MNR  1 KNR1   738 KNR2   737 KNR3   731 KNR4   938 $$\n");
%!     fputs(fd, "                                  KNR5   937 KNR6   931\n");
%!     fputs(fd, "STRU BNR  2 ENR 130 TNR 12 MNR  1 KNR1   639 KNR2   738 KNR3   731 KNR4   632 $$\n");
%!     fputs(fd, "                                  KNR5   839 KNR6   938 KNR7   931 KNR8   832\n");
%!     fputs(fd, "STRU BNR  2 ENR 131 TNR 12 MNR  1 KNR1   633 KNR2   734 KNR3   741 KNR4   640 $$\n");
%!     fputs(fd, "                                  KNR5   833 KNR6   934 KNR7   941 KNR8   840\n");
%!     fputs(fd, "STRU BNR  2 ENR 132 TNR 12 MNR  1 KNR1   735 KNR2   742 KNR3   741 KNR4   734 $$\n");
%!     fputs(fd, "                                  KNR5   935 KNR6   942 KNR7   941 KNR8   934\n");
%!     fputs(fd, "STRU BNR  2 ENR 133 TNR 12 MNR  1 KNR1   736 KNR2   743 KNR3   742 KNR4   735 $$\n");
%!     fputs(fd, "                                  KNR5   936 KNR6   943 KNR7   942 KNR8   935\n");
%!     fputs(fd, "STRU BNR  2 ENR 134 TNR 12 MNR  1 KNR1   737 KNR2   744 KNR3   743 KNR4   736 $$\n");
%!     fputs(fd, "                                  KNR5   937 KNR6   944 KNR7   943 KNR8   936\n");
%!     fputs(fd, "STRU BNR  2 ENR 135 TNR 12 MNR  1 KNR1   738 KNR2   745 KNR3   744 KNR4   737 $$\n");
%!     fputs(fd, "                                  KNR5   938 KNR6   945 KNR7   944 KNR8   937\n");
%!     fputs(fd, "STRU BNR  2 ENR 141 TNR 11 MNR  1 KNR1   615 KNR2   750 KNR3   639 KNR4   815 $$\n");
%!     fputs(fd, "                                  KNR5   950 KNR6   839\n");
%!     fputs(fd, "STRU BNR  2 ENR 142 TNR 11 MNR  1 KNR1   615 KNR2   714 KNR3   750 KNR4   815 $$\n");
%!     fputs(fd, "                                  KNR5   914 KNR6   950\n");
%!     fputs(fd, "STRU BNR  2 ENR 143 TNR 11 MNR  1 KNR1   750 KNR2   738 KNR3   639 KNR4   950 $$\n");
%!     fputs(fd, "                                  KNR5   938 KNR6   839\n");
%!     fputs(fd, "STRU BNR  2 ENR 144 TNR 11 MNR  1 KNR1   750 KNR2   714 KNR3   738 KNR4   950 $$\n");
%!     fputs(fd, "                                  KNR5   914 KNR6   938\n");
%!     fputs(fd, "STRU BNR  2 ENR 145 TNR 12 MNR  1 KNR1   714 KNR2   721 KNR3   745 KNR4   738 $$\n");
%!     fputs(fd, "                                  KNR5   914 KNR6   921 KNR7   945 KNR8   938\n");
%!     fputs(fd, "STRU BNR  2 ENR 201 TNR 12 MNR  1 KNR1   802 KNR2   803 KNR3   804 KNR4   801 $$\n");
%!     fputs(fd, "                                  KNR5  1002 KNR6  1003 KNR7  1004 KNR8  1001\n");
%!     fputs(fd, "STRU BNR  2 ENR 202 TNR 11 MNR  1 KNR1   801 KNR2   804 KNR3   805 KNR4  1001 $$\n");
%!     fputs(fd, "                                  KNR5  1004 KNR6  1005\n");
%!     fputs(fd, "STRU BNR  2 ENR 203 TNR 12 MNR  1 KNR1   801 KNR2   806 KNR3   807 KNR4   808 $$\n");
%!     fputs(fd, "                                  KNR5  1001 KNR6  1006 KNR7  1007 KNR8  1008\n");
%!     fputs(fd, "STRU BNR  2 ENR 204 TNR 12 MNR  1 KNR1   809 KNR2   810 KNR3   803 KNR4   802 $$\n");
%!     fputs(fd, "                                  KNR5  1009 KNR6  1010 KNR7  1003 KNR8  1002\n");
%!     fputs(fd, "STRU BNR  2 ENR 205 TNR 12 MNR  1 KNR1   810 KNR2   811 KNR3   804 KNR4   803 $$\n");
%!     fputs(fd, "                                  KNR5  1010 KNR6  1011 KNR7  1004 KNR8  1003\n");
%!     fputs(fd, "STRU BNR  2 ENR 206 TNR 12 MNR  1 KNR1   804 KNR2   811 KNR3   812 KNR4   805 $$\n");
%!     fputs(fd, "                                  KNR5  1004 KNR6  1011 KNR7  1012 KNR8  1005\n");
%!     fputs(fd, "STRU BNR  2 ENR 207 TNR 11 MNR  1 KNR1   801 KNR2   805 KNR3   806 KNR4  1001 $$\n");
%!     fputs(fd, "                                  KNR5  1005 KNR6  1006\n");
%!     fputs(fd, "STRU BNR  2 ENR 208 TNR 11 MNR  1 KNR1   805 KNR2   812 KNR3   822 KNR4  1005 $$\n");
%!     fputs(fd, "                                  KNR5  1012 KNR6  1022\n");
%!     fputs(fd, "STRU BNR  2 ENR 209 TNR 12 MNR  1 KNR1   805 KNR2   822 KNR3   813 KNR4   806 $$\n");
%!     fputs(fd, "                                  KNR5  1005 KNR6  1022 KNR7  1013 KNR8  1006\n");
%!     fputs(fd, "STRU BNR  2 ENR 210 TNR 12 MNR  1 KNR1   808 KNR2   807 KNR3   814 KNR4   815 $$\n");
%!     fputs(fd, "                                  KNR5  1008 KNR6  1007 KNR7  1014 KNR8  1015\n");
%!     fputs(fd, "STRU BNR  2 ENR 211 TNR 12 MNR  1 KNR1   806 KNR2   813 KNR3   814 KNR4   807 $$\n");
%!     fputs(fd, "                                  KNR5  1006 KNR6  1013 KNR7  1014 KNR8  1007\n");
%!     fputs(fd, "STRU BNR  2 ENR 212 TNR 12 MNR  1 KNR1   816 KNR2   817 KNR3   810 KNR4   809 $$\n");
%!     fputs(fd, "                                  KNR5  1016 KNR6  1017 KNR7  1010 KNR8  1009\n");
%!     fputs(fd, "STRU BNR  2 ENR 213 TNR 12 MNR  1 KNR1   817 KNR2   818 KNR3   811 KNR4   810 $$\n");
%!     fputs(fd, "                                  KNR5  1017 KNR6  1018 KNR7  1011 KNR8  1010\n");
%!     fputs(fd, "STRU BNR  2 ENR 214 TNR 12 MNR  1 KNR1   811 KNR2   818 KNR3   819 KNR4   812 $$\n");
%!     fputs(fd, "                                  KNR5  1011 KNR6  1018 KNR7  1019 KNR8  1012\n");
%!     fputs(fd, "STRU BNR  2 ENR 215 TNR 12 MNR  1 KNR1   812 KNR2   819 KNR3   820 KNR4   822 $$\n");
%!     fputs(fd, "                                  KNR5  1012 KNR6  1019 KNR7  1020 KNR8  1022\n");
%!     fputs(fd, "STRU BNR  2 ENR 216 TNR 12 MNR  1 KNR1   822 KNR2   820 KNR3   821 KNR4   813 $$\n");
%!     fputs(fd, "                                  KNR5  1022 KNR6  1020 KNR7  1021 KNR8  1013\n");
%!     fputs(fd, "STRU BNR  2 ENR 217 TNR 11 MNR  1 KNR1   813 KNR2   821 KNR3   814 KNR4  1013 $$\n");
%!     fputs(fd, "                                  KNR5  1021 KNR6  1014\n");
%!     fputs(fd, "STRU BNR  2 ENR 221 TNR 12 MNR  1 KNR1   825 KNR2   828 KNR3   827 KNR4   826 $$\n");
%!     fputs(fd, "                                  KNR5  1025 KNR6  1028 KNR7  1027 KNR8  1026\n");
%!     fputs(fd, "STRU BNR  2 ENR 222 TNR 12 MNR  1 KNR1   825 KNR2   830 KNR3   829 KNR4   828 $$\n");
%!     fputs(fd, "                                  KNR5  1025 KNR6  1030 KNR7  1029 KNR8  1028\n");
%!     fputs(fd, "STRU BNR  2 ENR 223 TNR 12 MNR  1 KNR1   832 KNR2   831 KNR3   830 KNR4   825 $$\n");
%!     fputs(fd, "                                  KNR5  1032 KNR6  1031 KNR7  1030 KNR8  1025\n");
%!     fputs(fd, "STRU BNR  2 ENR 224 TNR 12 MNR  1 KNR1   826 KNR2   827 KNR3   834 KNR4   833 $$\n");
%!     fputs(fd, "                                  KNR5  1026 KNR6  1027 KNR7  1034 KNR8  1033\n");
%!     fputs(fd, "STRU BNR  2 ENR 225 TNR 12 MNR  1 KNR1   828 KNR2   835 KNR3   834 KNR4   827 $$\n");
%!     fputs(fd, "                                  KNR5  1028 KNR6  1035 KNR7  1034 KNR8  1027\n");
%!     fputs(fd, "STRU BNR  2 ENR 226 TNR 12 MNR  1 KNR1   829 KNR2   836 KNR3   835 KNR4   828 $$\n");
%!     fputs(fd, "                                  KNR5  1029 KNR6  1036 KNR7  1035 KNR8  1028\n");
%!     fputs(fd, "STRU BNR  2 ENR 227 TNR 12 MNR  1 KNR1   830 KNR2   837 KNR3   836 KNR4   829 $$\n");
%!     fputs(fd, "                                  KNR5  1030 KNR6  1037 KNR7  1036 KNR8  1029\n");
%!     fputs(fd, "STRU BNR  2 ENR 228 TNR 11 MNR  1 KNR1   831 KNR2   837 KNR3   830 KNR4  1031 $$\n");
%!     fputs(fd, "                                  KNR5  1037 KNR6  1030\n");
%!     fputs(fd, "STRU BNR  2 ENR 229 TNR 11 MNR  1 KNR1   838 KNR2   837 KNR3   831 KNR4  1038 $$\n");
%!     fputs(fd, "                                  KNR5  1037 KNR6  1031\n");
%!     fputs(fd, "STRU BNR  2 ENR 230 TNR 12 MNR  1 KNR1   839 KNR2   838 KNR3   831 KNR4   832 $$\n");
%!     fputs(fd, "                                  KNR5  1039 KNR6  1038 KNR7  1031 KNR8  1032\n");
%!     fputs(fd, "STRU BNR  2 ENR 231 TNR 12 MNR  1 KNR1   833 KNR2   834 KNR3   841 KNR4   840 $$\n");
%!     fputs(fd, "                                  KNR5  1033 KNR6  1034 KNR7  1041 KNR8  1040\n");
%!     fputs(fd, "STRU BNR  2 ENR 232 TNR 12 MNR  1 KNR1   835 KNR2   842 KNR3   841 KNR4   834 $$\n");
%!     fputs(fd, "                                  KNR5  1035 KNR6  1042 KNR7  1041 KNR8  1034\n");
%!     fputs(fd, "STRU BNR  2 ENR 233 TNR 12 MNR  1 KNR1   836 KNR2   843 KNR3   842 KNR4   835 $$\n");
%!     fputs(fd, "                                  KNR5  1036 KNR6  1043 KNR7  1042 KNR8  1035\n");
%!     fputs(fd, "STRU BNR  2 ENR 234 TNR 12 MNR  1 KNR1   837 KNR2   844 KNR3   843 KNR4   836 $$\n");
%!     fputs(fd, "                                  KNR5  1037 KNR6  1044 KNR7  1043 KNR8  1036\n");
%!     fputs(fd, "STRU BNR  2 ENR 235 TNR 12 MNR  1 KNR1   838 KNR2   845 KNR3   844 KNR4   837 $$\n");
%!     fputs(fd, "                                  KNR5  1038 KNR6  1045 KNR7  1044 KNR8  1037\n");
%!     fputs(fd, "STRU BNR  2 ENR 241 TNR 11 MNR  1 KNR1   815 KNR2   850 KNR3   839 KNR4  1015 $$\n");
%!     fputs(fd, "                                  KNR5  1050 KNR6  1039\n");
%!     fputs(fd, "STRU BNR  2 ENR 242 TNR 11 MNR  1 KNR1   815 KNR2   814 KNR3   850 KNR4  1015 $$\n");
%!     fputs(fd, "                                  KNR5  1014 KNR6  1050\n");
%!     fputs(fd, "STRU BNR  2 ENR 243 TNR 11 MNR  1 KNR1   850 KNR2   838 KNR3   839 KNR4  1050 $$\n");
%!     fputs(fd, "                                  KNR5  1038 KNR6  1039\n");
%!     fputs(fd, "STRU BNR  2 ENR 244 TNR 11 MNR  1 KNR1   850 KNR2   814 KNR3   838 KNR4  1050 $$\n");
%!     fputs(fd, "                                  KNR5  1014 KNR6  1038\n");
%!     fputs(fd, "STRU BNR  2 ENR 245 TNR 12 MNR  1 KNR1   814 KNR2   821 KNR3   845 KNR4   838 $$\n");
%!     fputs(fd, "                                  KNR5  1014 KNR6  1021 KNR7  1045 KNR8  1038\n");
%!     fputs(fd, "STRU BNR  2 ENR 301 TNR 12 MNR  1 KNR1   802 KNR2   903 KNR3   904 KNR4   801 $$\n");
%!     fputs(fd, "                                  KNR5  1002 KNR6  1103 KNR7  1104 KNR8  1001\n");
%!     fputs(fd, "STRU BNR  2 ENR 302 TNR 11 MNR  1 KNR1   801 KNR2   904 KNR3   905 KNR4  1001 $$\n");
%!     fputs(fd, "                                  KNR5  1104 KNR6  1105\n");
%!     fputs(fd, "STRU BNR  2 ENR 303 TNR 12 MNR  1 KNR1   801 KNR2   906 KNR3   907 KNR4   808 $$\n");
%!     fputs(fd, "                                  KNR5  1001 KNR6  1106 KNR7  1107 KNR8  1008\n");
%!     fputs(fd, "STRU BNR  2 ENR 304 TNR 12 MNR  1 KNR1   809 KNR2   910 KNR3   903 KNR4   802 $$\n");
%!     fputs(fd, "                                  KNR5  1009 KNR6  1110 KNR7  1103 KNR8  1002\n");
%!     fputs(fd, "STRU BNR  2 ENR 305 TNR 12 MNR  1 KNR1   910 KNR2   911 KNR3   904 KNR4   903 $$\n");
%!     fputs(fd, "                                  KNR5  1110 KNR6  1111 KNR7  1104 KNR8  1103\n");
%!     fputs(fd, "STRU BNR  2 ENR 306 TNR 12 MNR  1 KNR1   904 KNR2   911 KNR3   912 KNR4   905 $$\n");
%!     fputs(fd, "                                  KNR5  1104 KNR6  1111 KNR7  1112 KNR8  1105\n");
%!     fputs(fd, "STRU BNR  2 ENR 307 TNR 11 MNR  1 KNR1   801 KNR2   905 KNR3   906 KNR4  1001 $$\n");
%!     fputs(fd, "                                  KNR5  1105 KNR6  1106\n");
%!     fputs(fd, "STRU BNR  2 ENR 308 TNR 11 MNR  1 KNR1   905 KNR2   912 KNR3   922 KNR4  1105 $$\n");
%!     fputs(fd, "                                  KNR5  1112 KNR6  1122\n");
%!     fputs(fd, "STRU BNR  2 ENR 309 TNR 12 MNR  1 KNR1   905 KNR2   922 KNR3   913 KNR4   906 $$\n");
%!     fputs(fd, "                                  KNR5  1105 KNR6  1122 KNR7  1113 KNR8  1106\n");
%!     fputs(fd, "STRU BNR  2 ENR 310 TNR 12 MNR  1 KNR1   808 KNR2   907 KNR3   914 KNR4   815 $$\n");
%!     fputs(fd, "                                  KNR5  1008 KNR6  1107 KNR7  1114 KNR8  1015\n");
%!     fputs(fd, "STRU BNR  2 ENR 311 TNR 12 MNR  1 KNR1   906 KNR2   913 KNR3   914 KNR4   907 $$\n");
%!     fputs(fd, "                                  KNR5  1106 KNR6  1113 KNR7  1114 KNR8  1107\n");
%!     fputs(fd, "STRU BNR  2 ENR 312 TNR 12 MNR  1 KNR1   816 KNR2   917 KNR3   910 KNR4   809 $$\n");
%!     fputs(fd, "                                  KNR5  1016 KNR6  1117 KNR7  1110 KNR8  1009\n");
%!     fputs(fd, "STRU BNR  2 ENR 313 TNR 12 MNR  1 KNR1   917 KNR2   918 KNR3   911 KNR4   910 $$\n");
%!     fputs(fd, "                                  KNR5  1117 KNR6  1118 KNR7  1111 KNR8  1110\n");
%!     fputs(fd, "STRU BNR  2 ENR 314 TNR 12 MNR  1 KNR1   911 KNR2   918 KNR3   919 KNR4   912 $$\n");
%!     fputs(fd, "                                  KNR5  1111 KNR6  1118 KNR7  1119 KNR8  1112\n");
%!     fputs(fd, "STRU BNR  2 ENR 315 TNR 12 MNR  1 KNR1   912 KNR2   919 KNR3   920 KNR4   922 $$\n");
%!     fputs(fd, "                                  KNR5  1112 KNR6  1119 KNR7  1120 KNR8  1122\n");
%!     fputs(fd, "STRU BNR  2 ENR 316 TNR 12 MNR  1 KNR1   922 KNR2   920 KNR3   921 KNR4   913 $$\n");
%!     fputs(fd, "                                  KNR5  1122 KNR6  1120 KNR7  1121 KNR8  1113\n");
%!     fputs(fd, "STRU BNR  2 ENR 317 TNR 11 MNR  1 KNR1   913 KNR2   921 KNR3   914 KNR4  1113 $$\n");
%!     fputs(fd, "                                  KNR5  1121 KNR6  1114\n");
%!     fputs(fd, "STRU BNR  2 ENR 321 TNR 12 MNR  1 KNR1   825 KNR2   928 KNR3   927 KNR4   826 $$\n");
%!     fputs(fd, "                                  KNR5  1025 KNR6  1128 KNR7  1127 KNR8  1026\n");
%!     fputs(fd, "STRU BNR  2 ENR 322 TNR 12 MNR  1 KNR1   825 KNR2   930 KNR3   929 KNR4   928 $$\n");
%!     fputs(fd, "                                  KNR5  1025 KNR6  1130 KNR7  1129 KNR8  1128\n");
%!     fputs(fd, "STRU BNR  2 ENR 323 TNR 12 MNR  1 KNR1   832 KNR2   931 KNR3   930 KNR4   825 $$\n");
%!     fputs(fd, "                                  KNR5  1032 KNR6  1131 KNR7  1130 KNR8  1025\n");
%!     fputs(fd, "STRU BNR  2 ENR 324 TNR 12 MNR  1 KNR1   826 KNR2   927 KNR3   934 KNR4   833 $$\n");
%!     fputs(fd, "                                  KNR5  1026 KNR6  1127 KNR7  1134 KNR8  1033\n");
%!     fputs(fd, "STRU BNR  2 ENR 325 TNR 12 MNR  1 KNR1   928 KNR2   935 KNR3   934 KNR4   927 $$\n");
%!     fputs(fd, "                                  KNR5  1128 KNR6  1135 KNR7  1134 KNR8  1127\n");
%!     fputs(fd, "STRU BNR  2 ENR 326 TNR 12 MNR  1 KNR1   929 KNR2   936 KNR3   935 KNR4   928 $$\n");
%!     fputs(fd, "                                  KNR5  1129 KNR6  1136 KNR7  1135 KNR8  1128\n");
%!     fputs(fd, "STRU BNR  2 ENR 327 TNR 12 MNR  1 KNR1   930 KNR2   937 KNR3   936 KNR4   929 $$\n");
%!     fputs(fd, "                                  KNR5  1130 KNR6  1137 KNR7  1136 KNR8  1129\n");
%!     fputs(fd, "STRU BNR  2 ENR 328 TNR 11 MNR  1 KNR1   931 KNR2   937 KNR3   930 KNR4  1131 $$\n");
%!     fputs(fd, "                                  KNR5  1137 KNR6  1130\n");
%!     fputs(fd, "STRU BNR  2 ENR 329 TNR 11 MNR  1 KNR1   938 KNR2   937 KNR3   931 KNR4  1138 $$\n");
%!     fputs(fd, "                                  KNR5  1137 KNR6  1131\n");
%!     fputs(fd, "STRU BNR  2 ENR 330 TNR 12 MNR  1 KNR1   839 KNR2   938 KNR3   931 KNR4   832 $$\n");
%!     fputs(fd, "                                  KNR5  1039 KNR6  1138 KNR7  1131 KNR8  1032\n");
%!     fputs(fd, "STRU BNR  2 ENR 331 TNR 12 MNR  1 KNR1   833 KNR2   934 KNR3   941 KNR4   840 $$\n");
%!     fputs(fd, "                                  KNR5  1033 KNR6  1134 KNR7  1141 KNR8  1040\n");
%!     fputs(fd, "STRU BNR  2 ENR 332 TNR 12 MNR  1 KNR1   935 KNR2   942 KNR3   941 KNR4   934 $$\n");
%!     fputs(fd, "                                  KNR5  1135 KNR6  1142 KNR7  1141 KNR8  1134\n");
%!     fputs(fd, "STRU BNR  2 ENR 333 TNR 12 MNR  1 KNR1   936 KNR2   943 KNR3   942 KNR4   935 $$\n");
%!     fputs(fd, "                                  KNR5  1136 KNR6  1143 KNR7  1142 KNR8  1135\n");
%!     fputs(fd, "STRU BNR  2 ENR 334 TNR 12 MNR  1 KNR1   937 KNR2   944 KNR3   943 KNR4   936 $$\n");
%!     fputs(fd, "                                  KNR5  1137 KNR6  1144 KNR7  1143 KNR8  1136\n");
%!     fputs(fd, "STRU BNR  2 ENR 335 TNR 12 MNR  1 KNR1   938 KNR2   945 KNR3   944 KNR4   937 $$\n");
%!     fputs(fd, "                                  KNR5  1138 KNR6  1145 KNR7  1144 KNR8  1137\n");
%!     fputs(fd, "STRU BNR  2 ENR 341 TNR 11 MNR  1 KNR1   815 KNR2   950 KNR3   839 KNR4  1015 $$\n");
%!     fputs(fd, "                                  KNR5  1150 KNR6  1039\n");
%!     fputs(fd, "STRU BNR  2 ENR 342 TNR 11 MNR  1 KNR1   815 KNR2   914 KNR3   950 KNR4  1015 $$\n");
%!     fputs(fd, "                                  KNR5  1114 KNR6  1150\n");
%!     fputs(fd, "STRU BNR  2 ENR 343 TNR 11 MNR  1 KNR1   950 KNR2   938 KNR3   839 KNR4  1150 $$\n");
%!     fputs(fd, "                                  KNR5  1138 KNR6  1039\n");
%!     fputs(fd, "STRU BNR  2 ENR 344 TNR 11 MNR  1 KNR1   950 KNR2   914 KNR3   938 KNR4  1150 $$\n");
%!     fputs(fd, "                                  KNR5  1114 KNR6  1138\n");
%!     fputs(fd, "STRU BNR  2 ENR 345 TNR 12 MNR  1 KNR1   914 KNR2   921 KNR3   945 KNR4   938 $$\n");
%!     fputs(fd, "                                  KNR5  1114 KNR6  1121 KNR7  1145 KNR8  1138\n");
%!     fputs(fd, "STRU BNR  2 ENR 401 TNR 12 MNR  1 KNR1  1002 KNR2  1003 KNR3  1004 KNR4  1001 $$\n");
%!     fputs(fd, "                                  KNR5  1202 KNR6  1203 KNR7  1204 KNR8  1201\n");
%!     fputs(fd, "STRU BNR  2 ENR 402 TNR 11 MNR  1 KNR1  1001 KNR2  1004 KNR3  1005 KNR4  1201 $$\n");
%!     fputs(fd, "                                  KNR5  1204 KNR6  1205\n");
%!     fputs(fd, "STRU BNR  2 ENR 403 TNR 12 MNR  1 KNR1  1001 KNR2  1006 KNR3  1007 KNR4  1008 $$\n");
%!     fputs(fd, "                                  KNR5  1201 KNR6  1206 KNR7  1207 KNR8  1208\n");
%!     fputs(fd, "STRU BNR  2 ENR 404 TNR 12 MNR  1 KNR1  1009 KNR2  1010 KNR3  1003 KNR4  1002 $$\n");
%!     fputs(fd, "                                  KNR5  1209 KNR6  1210 KNR7  1203 KNR8  1202\n");
%!     fputs(fd, "STRU BNR  2 ENR 405 TNR 12 MNR  1 KNR1  1010 KNR2  1011 KNR3  1004 KNR4  1003 $$\n");
%!     fputs(fd, "                                  KNR5  1210 KNR6  1211 KNR7  1204 KNR8  1203\n");
%!     fputs(fd, "STRU BNR  2 ENR 406 TNR 12 MNR  1 KNR1  1004 KNR2  1011 KNR3  1012 KNR4  1005 $$\n");
%!     fputs(fd, "                                  KNR5  1204 KNR6  1211 KNR7  1212 KNR8  1205\n");
%!     fputs(fd, "STRU BNR  2 ENR 407 TNR 11 MNR  1 KNR1  1001 KNR2  1005 KNR3  1006 KNR4  1201 $$\n");
%!     fputs(fd, "                                  KNR5  1205 KNR6  1206\n");
%!     fputs(fd, "STRU BNR  2 ENR 408 TNR 11 MNR  1 KNR1  1005 KNR2  1012 KNR3  1022 KNR4  1205 $$\n");
%!     fputs(fd, "                                  KNR5  1212 KNR6  1222\n");
%!     fputs(fd, "STRU BNR  2 ENR 409 TNR 12 MNR  1 KNR1  1005 KNR2  1022 KNR3  1013 KNR4  1006 $$\n");
%!     fputs(fd, "                                  KNR5  1205 KNR6  1222 KNR7  1213 KNR8  1206\n");
%!     fputs(fd, "STRU BNR  2 ENR 410 TNR 12 MNR  1 KNR1  1008 KNR2  1007 KNR3  1014 KNR4  1015 $$\n");
%!     fputs(fd, "                                  KNR5  1208 KNR6  1207 KNR7  1214 KNR8  1215\n");
%!     fputs(fd, "STRU BNR  2 ENR 411 TNR 12 MNR  1 KNR1  1006 KNR2  1013 KNR3  1014 KNR4  1007 $$\n");
%!     fputs(fd, "                                  KNR5  1206 KNR6  1213 KNR7  1214 KNR8  1207\n");
%!     fputs(fd, "STRU BNR  2 ENR 412 TNR 12 MNR  1 KNR1  1016 KNR2  1017 KNR3  1010 KNR4  1009 $$\n");
%!     fputs(fd, "                                  KNR5  1216 KNR6  1217 KNR7  1210 KNR8  1209\n");
%!     fputs(fd, "STRU BNR  2 ENR 413 TNR 12 MNR  1 KNR1  1017 KNR2  1018 KNR3  1011 KNR4  1010 $$\n");
%!     fputs(fd, "                                  KNR5  1217 KNR6  1218 KNR7  1211 KNR8  1210\n");
%!     fputs(fd, "STRU BNR  2 ENR 414 TNR 12 MNR  1 KNR1  1011 KNR2  1018 KNR3  1019 KNR4  1012 $$\n");
%!     fputs(fd, "                                  KNR5  1211 KNR6  1218 KNR7  1219 KNR8  1212\n");
%!     fputs(fd, "STRU BNR  2 ENR 415 TNR 12 MNR  1 KNR1  1012 KNR2  1019 KNR3  1020 KNR4  1022 $$\n");
%!     fputs(fd, "                                  KNR5  1212 KNR6  1219 KNR7  1220 KNR8  1222\n");
%!     fputs(fd, "STRU BNR  2 ENR 416 TNR 12 MNR  1 KNR1  1022 KNR2  1020 KNR3  1021 KNR4  1013 $$\n");
%!     fputs(fd, "                                  KNR5  1222 KNR6  1220 KNR7  1221 KNR8  1213\n");
%!     fputs(fd, "STRU BNR  2 ENR 417 TNR 11 MNR  1 KNR1  1013 KNR2  1021 KNR3  1014 KNR4  1213 $$\n");
%!     fputs(fd, "                                  KNR5  1221 KNR6  1214\n");
%!     fputs(fd, "STRU BNR  2 ENR 421 TNR 12 MNR  1 KNR1  1025 KNR2  1028 KNR3  1027 KNR4  1026 $$\n");
%!     fputs(fd, "                                  KNR5  1225 KNR6  1228 KNR7  1227 KNR8  1226\n");
%!     fputs(fd, "STRU BNR  2 ENR 422 TNR 12 MNR  1 KNR1  1025 KNR2  1030 KNR3  1029 KNR4  1028 $$\n");
%!     fputs(fd, "                                  KNR5  1225 KNR6  1230 KNR7  1229 KNR8  1228\n");
%!     fputs(fd, "STRU BNR  2 ENR 423 TNR 12 MNR  1 KNR1  1032 KNR2  1031 KNR3  1030 KNR4  1025 $$\n");
%!     fputs(fd, "                                  KNR5  1232 KNR6  1231 KNR7  1230 KNR8  1225\n");
%!     fputs(fd, "STRU BNR  2 ENR 424 TNR 12 MNR  1 KNR1  1026 KNR2  1027 KNR3  1034 KNR4  1033 $$\n");
%!     fputs(fd, "                                  KNR5  1226 KNR6  1227 KNR7  1234 KNR8  1233\n");
%!     fputs(fd, "STRU BNR  2 ENR 425 TNR 12 MNR  1 KNR1  1028 KNR2  1035 KNR3  1034 KNR4  1027 $$\n");
%!     fputs(fd, "                                  KNR5  1228 KNR6  1235 KNR7  1234 KNR8  1227\n");
%!     fputs(fd, "STRU BNR  2 ENR 426 TNR 12 MNR  1 KNR1  1029 KNR2  1036 KNR3  1035 KNR4  1028 $$\n");
%!     fputs(fd, "                                  KNR5  1229 KNR6  1236 KNR7  1235 KNR8  1228\n");
%!     fputs(fd, "STRU BNR  2 ENR 427 TNR 12 MNR  1 KNR1  1030 KNR2  1037 KNR3  1036 KNR4  1029 $$\n");
%!     fputs(fd, "                                  KNR5  1230 KNR6  1237 KNR7  1236 KNR8  1229\n");
%!     fputs(fd, "STRU BNR  2 ENR 428 TNR 11 MNR  1 KNR1  1031 KNR2  1037 KNR3  1030 KNR4  1231 $$\n");
%!     fputs(fd, "                                  KNR5  1237 KNR6  1230\n");
%!     fputs(fd, "STRU BNR  2 ENR 429 TNR 11 MNR  1 KNR1  1038 KNR2  1037 KNR3  1031 KNR4  1238 $$\n");
%!     fputs(fd, "                                  KNR5  1237 KNR6  1231\n");
%!     fputs(fd, "STRU BNR  2 ENR 430 TNR 12 MNR  1 KNR1  1039 KNR2  1038 KNR3  1031 KNR4  1032 $$\n");
%!     fputs(fd, "                                  KNR5  1239 KNR6  1238 KNR7  1231 KNR8  1232\n");
%!     fputs(fd, "STRU BNR  2 ENR 431 TNR 12 MNR  1 KNR1  1033 KNR2  1034 KNR3  1041 KNR4  1040 $$\n");
%!     fputs(fd, "                                  KNR5  1233 KNR6  1234 KNR7  1241 KNR8  1240\n");
%!     fputs(fd, "STRU BNR  2 ENR 432 TNR 12 MNR  1 KNR1  1035 KNR2  1042 KNR3  1041 KNR4  1034 $$\n");
%!     fputs(fd, "                                  KNR5  1235 KNR6  1242 KNR7  1241 KNR8  1234\n");
%!     fputs(fd, "STRU BNR  2 ENR 433 TNR 12 MNR  1 KNR1  1036 KNR2  1043 KNR3  1042 KNR4  1035 $$\n");
%!     fputs(fd, "                                  KNR5  1236 KNR6  1243 KNR7  1242 KNR8  1235\n");
%!     fputs(fd, "STRU BNR  2 ENR 434 TNR 12 MNR  1 KNR1  1037 KNR2  1044 KNR3  1043 KNR4  1036 $$\n");
%!     fputs(fd, "                                  KNR5  1237 KNR6  1244 KNR7  1243 KNR8  1236\n");
%!     fputs(fd, "STRU BNR  2 ENR 435 TNR 12 MNR  1 KNR1  1038 KNR2  1045 KNR3  1044 KNR4  1037 $$\n");
%!     fputs(fd, "                                  KNR5  1238 KNR6  1245 KNR7  1244 KNR8  1237\n");
%!     fputs(fd, "STRU BNR  2 ENR 441 TNR 11 MNR  1 KNR1  1015 KNR2  1050 KNR3  1039 KNR4  1215 $$\n");
%!     fputs(fd, "                                  KNR5  1250 KNR6  1239\n");
%!     fputs(fd, "STRU BNR  2 ENR 442 TNR 11 MNR  1 KNR1  1015 KNR2  1014 KNR3  1050 KNR4  1215 $$\n");
%!     fputs(fd, "                                  KNR5  1214 KNR6  1250\n");
%!     fputs(fd, "STRU BNR  2 ENR 443 TNR 11 MNR  1 KNR1  1050 KNR2  1038 KNR3  1039 KNR4  1250 $$\n");
%!     fputs(fd, "                                  KNR5  1238 KNR6  1239\n");
%!     fputs(fd, "STRU BNR  2 ENR 444 TNR 11 MNR  1 KNR1  1050 KNR2  1014 KNR3  1038 KNR4  1250 $$\n");
%!     fputs(fd, "                                  KNR5  1214 KNR6  1238\n");
%!     fputs(fd, "STRU BNR  2 ENR 445 TNR 12 MNR  1 KNR1  1014 KNR2  1021 KNR3  1045 KNR4  1038 $$\n");
%!     fputs(fd, "                                  KNR5  1214 KNR6  1221 KNR7  1245 KNR8  1238\n");
%!     fputs(fd, "STRU BNR  2 ENR 501 TNR 12 MNR  1 KNR1  1002 KNR2  1103 KNR3  1104 KNR4  1001 $$\n");
%!     fputs(fd, "                                  KNR5  1202 KNR6  1303 KNR7  1304 KNR8  1201\n");
%!     fputs(fd, "STRU BNR  2 ENR 502 TNR 11 MNR  1 KNR1  1001 KNR2  1104 KNR3  1105 KNR4  1201 $$\n");
%!     fputs(fd, "                                  KNR5  1304 KNR6  1305\n");
%!     fputs(fd, "STRU BNR  2 ENR 503 TNR 12 MNR  1 KNR1  1001 KNR2  1106 KNR3  1107 KNR4  1008 $$\n");
%!     fputs(fd, "                                  KNR5  1201 KNR6  1306 KNR7  1307 KNR8  1208\n");
%!     fputs(fd, "STRU BNR  2 ENR 504 TNR 12 MNR  1 KNR1  1009 KNR2  1110 KNR3  1103 KNR4  1002 $$\n");
%!     fputs(fd, "                                  KNR5  1209 KNR6  1310 KNR7  1303 KNR8  1202\n");
%!     fputs(fd, "STRU BNR  2 ENR 505 TNR 12 MNR  1 KNR1  1110 KNR2  1111 KNR3  1104 KNR4  1103 $$\n");
%!     fputs(fd, "                                  KNR5  1310 KNR6  1311 KNR7  1304 KNR8  1303\n");
%!     fputs(fd, "STRU BNR  2 ENR 506 TNR 12 MNR  1 KNR1  1104 KNR2  1111 KNR3  1112 KNR4  1105 $$\n");
%!     fputs(fd, "                                  KNR5  1304 KNR6  1311 KNR7  1312 KNR8  1305\n");
%!     fputs(fd, "STRU BNR  2 ENR 507 TNR 11 MNR  1 KNR1  1001 KNR2  1105 KNR3  1106 KNR4  1201 $$\n");
%!     fputs(fd, "                                  KNR5  1305 KNR6  1306\n");
%!     fputs(fd, "STRU BNR  2 ENR 508 TNR 11 MNR  1 KNR1  1105 KNR2  1112 KNR3  1122 KNR4  1305 $$\n");
%!     fputs(fd, "                                  KNR5  1312 KNR6  1322\n");
%!     fputs(fd, "STRU BNR  2 ENR 509 TNR 12 MNR  1 KNR1  1105 KNR2  1122 KNR3  1113 KNR4  1106 $$\n");
%!     fputs(fd, "                                  KNR5  1305 KNR6  1322 KNR7  1313 KNR8  1306\n");
%!     fputs(fd, "STRU BNR  2 ENR 510 TNR 12 MNR  1 KNR1  1008 KNR2  1107 KNR3  1114 KNR4  1015 $$\n");
%!     fputs(fd, "                                  KNR5  1208 KNR6  1307 KNR7  1314 KNR8  1215\n");
%!     fputs(fd, "STRU BNR  2 ENR 511 TNR 12 MNR  1 KNR1  1106 KNR2  1113 KNR3  1114 KNR4  1107 $$\n");
%!     fputs(fd, "                                  KNR5  1306 KNR6  1313 KNR7  1314 KNR8  1307\n");
%!     fputs(fd, "STRU BNR  2 ENR 512 TNR 12 MNR  1 KNR1  1016 KNR2  1117 KNR3  1110 KNR4  1009 $$\n");
%!     fputs(fd, "                                  KNR5  1216 KNR6  1317 KNR7  1310 KNR8  1209\n");
%!     fputs(fd, "STRU BNR  2 ENR 513 TNR 12 MNR  1 KNR1  1117 KNR2  1118 KNR3  1111 KNR4  1110 $$\n");
%!     fputs(fd, "                                  KNR5  1317 KNR6  1318 KNR7  1311 KNR8  1310\n");
%!     fputs(fd, "STRU BNR  2 ENR 514 TNR 12 MNR  1 KNR1  1111 KNR2  1118 KNR3  1119 KNR4  1112 $$\n");
%!     fputs(fd, "                                  KNR5  1311 KNR6  1318 KNR7  1319 KNR8  1312\n");
%!     fputs(fd, "STRU BNR  2 ENR 515 TNR 12 MNR  1 KNR1  1112 KNR2  1119 KNR3  1120 KNR4  1122 $$\n");
%!     fputs(fd, "                                  KNR5  1312 KNR6  1319 KNR7  1320 KNR8  1322\n");
%!     fputs(fd, "STRU BNR  2 ENR 516 TNR 12 MNR  1 KNR1  1122 KNR2  1120 KNR3  1121 KNR4  1113 $$\n");
%!     fputs(fd, "                                  KNR5  1322 KNR6  1320 KNR7  1321 KNR8  1313\n");
%!     fputs(fd, "STRU BNR  2 ENR 517 TNR 11 MNR  1 KNR1  1113 KNR2  1121 KNR3  1114 KNR4  1313 $$\n");
%!     fputs(fd, "                                  KNR5  1321 KNR6  1314\n");
%!     fputs(fd, "STRU BNR  2 ENR 521 TNR 12 MNR  1 KNR1  1025 KNR2  1128 KNR3  1127 KNR4  1026 $$\n");
%!     fputs(fd, "                                  KNR5  1225 KNR6  1328 KNR7  1327 KNR8  1226\n");
%!     fputs(fd, "STRU BNR  2 ENR 522 TNR 12 MNR  1 KNR1  1025 KNR2  1130 KNR3  1129 KNR4  1128 $$\n");
%!     fputs(fd, "                                  KNR5  1225 KNR6  1330 KNR7  1329 KNR8  1328\n");
%!     fputs(fd, "STRU BNR  2 ENR 523 TNR 12 MNR  1 KNR1  1032 KNR2  1131 KNR3  1130 KNR4  1025 $$\n");
%!     fputs(fd, "                                  KNR5  1232 KNR6  1331 KNR7  1330 KNR8  1225\n");
%!     fputs(fd, "STRU BNR  2 ENR 524 TNR 12 MNR  1 KNR1  1026 KNR2  1127 KNR3  1134 KNR4  1033 $$\n");
%!     fputs(fd, "                                  KNR5  1226 KNR6  1327 KNR7  1334 KNR8  1233\n");
%!     fputs(fd, "STRU BNR  2 ENR 525 TNR 12 MNR  1 KNR1  1128 KNR2  1135 KNR3  1134 KNR4  1127 $$\n");
%!     fputs(fd, "                                  KNR5  1328 KNR6  1335 KNR7  1334 KNR8  1327\n");
%!     fputs(fd, "STRU BNR  2 ENR 526 TNR 12 MNR  1 KNR1  1129 KNR2  1136 KNR3  1135 KNR4  1128 $$\n");
%!     fputs(fd, "                                  KNR5  1329 KNR6  1336 KNR7  1335 KNR8  1328\n");
%!     fputs(fd, "STRU BNR  2 ENR 527 TNR 12 MNR  1 KNR1  1130 KNR2  1137 KNR3  1136 KNR4  1129 $$\n");
%!     fputs(fd, "                                  KNR5  1330 KNR6  1337 KNR7  1336 KNR8  1329\n");
%!     fputs(fd, "STRU BNR  2 ENR 528 TNR 11 MNR  1 KNR1  1131 KNR2  1137 KNR3  1130 KNR4  1331 $$\n");
%!     fputs(fd, "                                  KNR5  1337 KNR6  1330\n");
%!     fputs(fd, "STRU BNR  2 ENR 529 TNR 11 MNR  1 KNR1  1138 KNR2  1137 KNR3  1131 KNR4  1338 $$\n");
%!     fputs(fd, "                                  KNR5  1337 KNR6  1331\n");
%!     fputs(fd, "STRU BNR  2 ENR 530 TNR 12 MNR  1 KNR1  1039 KNR2  1138 KNR3  1131 KNR4  1032 $$\n");
%!     fputs(fd, "                                  KNR5  1239 KNR6  1338 KNR7  1331 KNR8  1232\n");
%!     fputs(fd, "STRU BNR  2 ENR 531 TNR 12 MNR  1 KNR1  1033 KNR2  1134 KNR3  1141 KNR4  1040 $$\n");
%!     fputs(fd, "                                  KNR5  1233 KNR6  1334 KNR7  1341 KNR8  1240\n");
%!     fputs(fd, "STRU BNR  2 ENR 532 TNR 12 MNR  1 KNR1  1135 KNR2  1142 KNR3  1141 KNR4  1134 $$\n");
%!     fputs(fd, "                                  KNR5  1335 KNR6  1342 KNR7  1341 KNR8  1334\n");
%!     fputs(fd, "STRU BNR  2 ENR 533 TNR 12 MNR  1 KNR1  1136 KNR2  1143 KNR3  1142 KNR4  1135 $$\n");
%!     fputs(fd, "                                  KNR5  1336 KNR6  1343 KNR7  1342 KNR8  1335\n");
%!     fputs(fd, "STRU BNR  2 ENR 534 TNR 12 MNR  1 KNR1  1137 KNR2  1144 KNR3  1143 KNR4  1136 $$\n");
%!     fputs(fd, "                                  KNR5  1337 KNR6  1344 KNR7  1343 KNR8  1336\n");
%!     fputs(fd, "STRU BNR  2 ENR 535 TNR 12 MNR  1 KNR1  1138 KNR2  1145 KNR3  1144 KNR4  1137 $$\n");
%!     fputs(fd, "                                  KNR5  1338 KNR6  1345 KNR7  1344 KNR8  1337\n");
%!     fputs(fd, "STRU BNR  2 ENR 541 TNR 11 MNR  1 KNR1  1015 KNR2  1150 KNR3  1039 KNR4  1215 $$\n");
%!     fputs(fd, "                                  KNR5  1350 KNR6  1239\n");
%!     fputs(fd, "STRU BNR  2 ENR 542 TNR 11 MNR  1 KNR1  1015 KNR2  1114 KNR3  1150 KNR4  1215 $$\n");
%!     fputs(fd, "                                  KNR5  1314 KNR6  1350\n");
%!     fputs(fd, "STRU BNR  2 ENR 543 TNR 11 MNR  1 KNR1  1150 KNR2  1138 KNR3  1039 KNR4  1350 $$\n");
%!     fputs(fd, "                                  KNR5  1338 KNR6  1239\n");
%!     fputs(fd, "STRU BNR  2 ENR 544 TNR 11 MNR  1 KNR1  1150 KNR2  1114 KNR3  1138 KNR4  1350 $$\n");
%!     fputs(fd, "                                  KNR5  1314 KNR6  1338\n");
%!     fputs(fd, "STRU BNR  2 ENR 545 TNR 12 MNR  1 KNR1  1114 KNR2  1121 KNR3  1145 KNR4  1138 $$\n");
%!     fputs(fd, "                                  KNR5  1314 KNR6  1321 KNR7  1345 KNR8  1338\n");
%!     fputs(fd, "STRU BNR  3 ENR  21 TNR 12 MNR  1 KNR1  1225 KNR2  1228 KNR3  1227 KNR4  1226 $$\n");
%!     fputs(fd, "                                  KNR5  1425 KNR6  1428 KNR7  1427 KNR8  1426\n");
%!     fputs(fd, "STRU BNR  3 ENR  22 TNR 12 MNR  1 KNR1  1225 KNR2  1230 KNR3  1229 KNR4  1228 $$\n");
%!     fputs(fd, "                                  KNR5  1425 KNR6  1430 KNR7  1429 KNR8  1428\n");
%!     fputs(fd, "STRU BNR  3 ENR  23 TNR 12 MNR  1 KNR1  1232 KNR2  1231 KNR3  1230 KNR4  1225 $$\n");
%!     fputs(fd, "                                  KNR5  1432 KNR6  1431 KNR7  1430 KNR8  1425\n");
%!     fputs(fd, "STRU BNR  3 ENR  24 TNR 12 MNR  1 KNR1  1226 KNR2  1227 KNR3  1234 KNR4  1233 $$\n");
%!     fputs(fd, "                                  KNR5  1426 KNR6  1427 KNR7  1434 KNR8  1433\n");
%!     fputs(fd, "STRU BNR  3 ENR  25 TNR 12 MNR  1 KNR1  1228 KNR2  1235 KNR3  1234 KNR4  1227 $$\n");
%!     fputs(fd, "                                  KNR5  1428 KNR6  1435 KNR7  1434 KNR8  1427\n");
%!     fputs(fd, "STRU BNR  3 ENR  26 TNR 12 MNR  1 KNR1  1229 KNR2  1236 KNR3  1235 KNR4  1228 $$\n");
%!     fputs(fd, "                                  KNR5  1429 KNR6  1436 KNR7  1435 KNR8  1428\n");
%!     fputs(fd, "STRU BNR  3 ENR  27 TNR 12 MNR  1 KNR1  1230 KNR2  1237 KNR3  1236 KNR4  1229 $$\n");
%!     fputs(fd, "                                  KNR5  1430 KNR6  1437 KNR7  1436 KNR8  1429\n");
%!     fputs(fd, "STRU BNR  3 ENR  28 TNR 11 MNR  1 KNR1  1231 KNR2  1237 KNR3  1230 KNR4  1431 $$\n");
%!     fputs(fd, "                                  KNR5  1437 KNR6  1430\n");
%!     fputs(fd, "STRU BNR  3 ENR  29 TNR 11 MNR  1 KNR1  1238 KNR2  1237 KNR3  1231 KNR4  1438 $$\n");
%!     fputs(fd, "                                  KNR5  1437 KNR6  1431\n");
%!     fputs(fd, "STRU BNR  3 ENR  30 TNR 12 MNR  1 KNR1  1239 KNR2  1238 KNR3  1231 KNR4  1232 $$\n");
%!     fputs(fd, "                                  KNR5  1439 KNR6  1438 KNR7  1431 KNR8  1432\n");
%!     fputs(fd, "STRU BNR  3 ENR  41 TNR 11 MNR  1 KNR1  1215 KNR2  1250 KNR3  1239 KNR4  1415 $$\n");
%!     fputs(fd, "                                  KNR5  1450 KNR6  1439\n");
%!     fputs(fd, "STRU BNR  3 ENR  43 TNR 11 MNR  1 KNR1  1250 KNR2  1238 KNR3  1239 KNR4  1450 $$\n");
%!     fputs(fd, "                                  KNR5  1438 KNR6  1439\n");
%!     fputs(fd, "STRU BNR  3 ENR 121 TNR 12 MNR  1 KNR1  1225 KNR2  1328 KNR3  1327 KNR4  1226 $$\n");
%!     fputs(fd, "                                  KNR5  1425 KNR6  1528 KNR7  1527 KNR8  1426\n");
%!     fputs(fd, "STRU BNR  3 ENR 122 TNR 12 MNR  1 KNR1  1225 KNR2  1330 KNR3  1329 KNR4  1328 $$\n");
%!     fputs(fd, "                                  KNR5  1425 KNR6  1530 KNR7  1529 KNR8  1528\n");
%!     fputs(fd, "STRU BNR  3 ENR 123 TNR 12 MNR  1 KNR1  1232 KNR2  1331 KNR3  1330 KNR4  1225 $$\n");
%!     fputs(fd, "                                  KNR5  1432 KNR6  1531 KNR7  1530 KNR8  1425\n");
%!     fputs(fd, "STRU BNR  3 ENR 124 TNR 12 MNR  1 KNR1  1226 KNR2  1327 KNR3  1334 KNR4  1233 $$\n");
%!     fputs(fd, "                                  KNR5  1426 KNR6  1527 KNR7  1534 KNR8  1433\n");
%!     fputs(fd, "STRU BNR  3 ENR 125 TNR 12 MNR  1 KNR1  1328 KNR2  1335 KNR3  1334 KNR4  1327 $$\n");
%!     fputs(fd, "                                  KNR5  1528 KNR6  1535 KNR7  1534 KNR8  1527\n");
%!     fputs(fd, "STRU BNR  3 ENR 126 TNR 12 MNR  1 KNR1  1329 KNR2  1336 KNR3  1335 KNR4  1328 $$\n");
%!     fputs(fd, "                                  KNR5  1529 KNR6  1536 KNR7  1535 KNR8  1528\n");
%!     fputs(fd, "STRU BNR  3 ENR 127 TNR 12 MNR  1 KNR1  1330 KNR2  1337 KNR3  1336 KNR4  1329 $$\n");
%!     fputs(fd, "                                  KNR5  1530 KNR6  1537 KNR7  1536 KNR8  1529\n");
%!     fputs(fd, "STRU BNR  3 ENR 128 TNR 11 MNR  1 KNR1  1331 KNR2  1337 KNR3  1330 KNR4  1531 $$\n");
%!     fputs(fd, "                                  KNR5  1537 KNR6  1530\n");
%!     fputs(fd, "STRU BNR  3 ENR 129 TNR 11 MNR  1 KNR1  1338 KNR2  1337 KNR3  1331 KNR4  1538 $$\n");
%!     fputs(fd, "                                  KNR5  1537 KNR6  1531\n");
%!     fputs(fd, "STRU BNR  3 ENR 130 TNR 12 MNR  1 KNR1  1239 KNR2  1338 KNR3  1331 KNR4  1232 $$\n");
%!     fputs(fd, "                                  KNR5  1439 KNR6  1538 KNR7  1531 KNR8  1432\n");
%!     fputs(fd, "STRU BNR  3 ENR 141 TNR 11 MNR  1 KNR1  1215 KNR2  1350 KNR3  1239 KNR4  1415 $$\n");
%!     fputs(fd, "                                  KNR5  1550 KNR6  1439\n");
%!     fputs(fd, "STRU BNR  3 ENR 143 TNR 11 MNR  1 KNR1  1350 KNR2  1338 KNR3  1239 KNR4  1550 $$\n");
%!     fputs(fd, "                                  KNR5  1538 KNR6  1439\n");
%!     fputs(fd, "STRU BNR  3 ENR 221 TNR 12 MNR  1 KNR1  1425 KNR2  1428 KNR3  1427 KNR4  1426 $$\n");
%!     fputs(fd, "                                  KNR5  1625 KNR6  1628 KNR7  1627 KNR8  1626\n");
%!     fputs(fd, "STRU BNR  3 ENR 222 TNR 12 MNR  1 KNR1  1425 KNR2  1430 KNR3  1429 KNR4  1428 $$\n");
%!     fputs(fd, "                                  KNR5  1625 KNR6  1630 KNR7  1629 KNR8  1628\n");
%!     fputs(fd, "STRU BNR  3 ENR 223 TNR 12 MNR  1 KNR1  1432 KNR2  1431 KNR3  1430 KNR4  1425 $$\n");
%!     fputs(fd, "                                  KNR5  1632 KNR6  1631 KNR7  1630 KNR8  1625\n");
%!     fputs(fd, "STRU BNR  3 ENR 224 TNR 12 MNR  1 KNR1  1426 KNR2  1427 KNR3  1434 KNR4  1433 $$\n");
%!     fputs(fd, "                                  KNR5  1626 KNR6  1627 KNR7  1634 KNR8  1633\n");
%!     fputs(fd, "STRU BNR  3 ENR 225 TNR 12 MNR  1 KNR1  1428 KNR2  1435 KNR3  1434 KNR4  1427 $$\n");
%!     fputs(fd, "                                  KNR5  1628 KNR6  1635 KNR7  1634 KNR8  1627\n");
%!     fputs(fd, "STRU BNR  3 ENR 226 TNR 12 MNR  1 KNR1  1429 KNR2  1436 KNR3  1435 KNR4  1428 $$\n");
%!     fputs(fd, "                                  KNR5  1629 KNR6  1636 KNR7  1635 KNR8  1628\n");
%!     fputs(fd, "STRU BNR  3 ENR 227 TNR 12 MNR  1 KNR1  1430 KNR2  1437 KNR3  1436 KNR4  1429 $$\n");
%!     fputs(fd, "                                  KNR5  1630 KNR6  1637 KNR7  1636 KNR8  1629\n");
%!     fputs(fd, "STRU BNR  3 ENR 228 TNR 11 MNR  1 KNR1  1431 KNR2  1437 KNR3  1430 KNR4  1631 $$\n");
%!     fputs(fd, "                                  KNR5  1637 KNR6  1630\n");
%!     fputs(fd, "STRU BNR  3 ENR 229 TNR 11 MNR  1 KNR1  1438 KNR2  1437 KNR3  1431 KNR4  1638 $$\n");
%!     fputs(fd, "                                  KNR5  1637 KNR6  1631\n");
%!     fputs(fd, "STRU BNR  3 ENR 230 TNR 12 MNR  1 KNR1  1439 KNR2  1438 KNR3  1431 KNR4  1432 $$\n");
%!     fputs(fd, "                                  KNR5  1639 KNR6  1638 KNR7  1631 KNR8  1632\n");
%!     fputs(fd, "STRU BNR  3 ENR 241 TNR 11 MNR  1 KNR1  1415 KNR2  1450 KNR3  1439 KNR4  1615 $$\n");
%!     fputs(fd, "                                  KNR5  1650 KNR6  1639\n");
%!     fputs(fd, "STRU BNR  3 ENR 243 TNR 11 MNR  1 KNR1  1450 KNR2  1438 KNR3  1439 KNR4  1650 $$\n");
%!     fputs(fd, "                                  KNR5  1638 KNR6  1639\n");
%!     fputs(fd, "STRU BNR  3 ENR 321 TNR 12 MNR  1 KNR1  1425 KNR2  1528 KNR3  1527 KNR4  1426 $$\n");
%!     fputs(fd, "                                  KNR5  1625 KNR6  1728 KNR7  1727 KNR8  1626\n");
%!     fputs(fd, "STRU BNR  3 ENR 322 TNR 12 MNR  1 KNR1  1425 KNR2  1530 KNR3  1529 KNR4  1528 $$\n");
%!     fputs(fd, "                                  KNR5  1625 KNR6  1730 KNR7  1729 KNR8  1728\n");
%!     fputs(fd, "STRU BNR  3 ENR 323 TNR 12 MNR  1 KNR1  1432 KNR2  1531 KNR3  1530 KNR4  1425 $$\n");
%!     fputs(fd, "                                  KNR5  1632 KNR6  1731 KNR7  1730 KNR8  1625\n");
%!     fputs(fd, "STRU BNR  3 ENR 324 TNR 12 MNR  1 KNR1  1426 KNR2  1527 KNR3  1534 KNR4  1433 $$\n");
%!     fputs(fd, "                                  KNR5  1626 KNR6  1727 KNR7  1734 KNR8  1633\n");
%!     fputs(fd, "STRU BNR  3 ENR 325 TNR 12 MNR  1 KNR1  1528 KNR2  1535 KNR3  1534 KNR4  1527 $$\n");
%!     fputs(fd, "                                  KNR5  1728 KNR6  1735 KNR7  1734 KNR8  1727\n");
%!     fputs(fd, "STRU BNR  3 ENR 326 TNR 12 MNR  1 KNR1  1529 KNR2  1536 KNR3  1535 KNR4  1528 $$\n");
%!     fputs(fd, "                                  KNR5  1729 KNR6  1736 KNR7  1735 KNR8  1728\n");
%!     fputs(fd, "STRU BNR  3 ENR 327 TNR 12 MNR  1 KNR1  1530 KNR2  1537 KNR3  1536 KNR4  1529 $$\n");
%!     fputs(fd, "                                  KNR5  1730 KNR6  1737 KNR7  1736 KNR8  1729\n");
%!     fputs(fd, "STRU BNR  3 ENR 328 TNR 11 MNR  1 KNR1  1531 KNR2  1537 KNR3  1530 KNR4  1731 $$\n");
%!     fputs(fd, "                                  KNR5  1737 KNR6  1730\n");
%!     fputs(fd, "STRU BNR  3 ENR 329 TNR 11 MNR  1 KNR1  1538 KNR2  1537 KNR3  1531 KNR4  1738 $$\n");
%!     fputs(fd, "                                  KNR5  1737 KNR6  1731\n");
%!     fputs(fd, "STRU BNR  3 ENR 330 TNR 12 MNR  1 KNR1  1439 KNR2  1538 KNR3  1531 KNR4  1432 $$\n");
%!     fputs(fd, "                                  KNR5  1639 KNR6  1738 KNR7  1731 KNR8  1632\n");
%!     fputs(fd, "STRU BNR  3 ENR 341 TNR 11 MNR  1 KNR1  1415 KNR2  1550 KNR3  1439 KNR4  1615 $$\n");
%!     fputs(fd, "                                  KNR5  1750 KNR6  1639\n");
%!     fputs(fd, "STRU BNR  3 ENR 343 TNR 11 MNR  1 KNR1  1550 KNR2  1538 KNR3  1439 KNR4  1750 $$\n");
%!     fputs(fd, "                                  KNR5  1738 KNR6  1639\n");
%!     fputs(fd, "STRU BNR  3 ENR 421 TNR 12 MNR  1 KNR1  1625 KNR2  1628 KNR3  1627 KNR4  1626 $$\n");
%!     fputs(fd, "                                  KNR5  1825 KNR6  1828 KNR7  1827 KNR8  1826\n");
%!     fputs(fd, "STRU BNR  3 ENR 422 TNR 12 MNR  1 KNR1  1625 KNR2  1630 KNR3  1629 KNR4  1628 $$\n");
%!     fputs(fd, "                                  KNR5  1825 KNR6  1830 KNR7  1829 KNR8  1828\n");
%!     fputs(fd, "STRU BNR  3 ENR 423 TNR 12 MNR  1 KNR1  1632 KNR2  1631 KNR3  1630 KNR4  1625 $$\n");
%!     fputs(fd, "                                  KNR5  1832 KNR6  1831 KNR7  1830 KNR8  1825\n");
%!     fputs(fd, "STRU BNR  3 ENR 424 TNR 12 MNR  1 KNR1  1626 KNR2  1627 KNR3  1634 KNR4  1633 $$\n");
%!     fputs(fd, "                                  KNR5  1826 KNR6  1827 KNR7  1834 KNR8  1833\n");
%!     fputs(fd, "STRU BNR  3 ENR 425 TNR 12 MNR  1 KNR1  1628 KNR2  1635 KNR3  1634 KNR4  1627 $$\n");
%!     fputs(fd, "                                  KNR5  1828 KNR6  1835 KNR7  1834 KNR8  1827\n");
%!     fputs(fd, "STRU BNR  3 ENR 426 TNR 12 MNR  1 KNR1  1629 KNR2  1636 KNR3  1635 KNR4  1628 $$\n");
%!     fputs(fd, "                                  KNR5  1829 KNR6  1836 KNR7  1835 KNR8  1828\n");
%!     fputs(fd, "STRU BNR  3 ENR 427 TNR 12 MNR  1 KNR1  1630 KNR2  1637 KNR3  1636 KNR4  1629 $$\n");
%!     fputs(fd, "                                  KNR5  1830 KNR6  1837 KNR7  1836 KNR8  1829\n");
%!     fputs(fd, "STRU BNR  3 ENR 428 TNR 11 MNR  1 KNR1  1631 KNR2  1637 KNR3  1630 KNR4  1831 $$\n");
%!     fputs(fd, "                                  KNR5  1837 KNR6  1830\n");
%!     fputs(fd, "STRU BNR  3 ENR 429 TNR 11 MNR  1 KNR1  1638 KNR2  1637 KNR3  1631 KNR4  1838 $$\n");
%!     fputs(fd, "                                  KNR5  1837 KNR6  1831\n");
%!     fputs(fd, "STRU BNR  3 ENR 430 TNR 12 MNR  1 KNR1  1639 KNR2  1638 KNR3  1631 KNR4  1632 $$\n");
%!     fputs(fd, "                                  KNR5  1839 KNR6  1838 KNR7  1831 KNR8  1832\n");
%!     fputs(fd, "STRU BNR  3 ENR 441 TNR 11 MNR  1 KNR1  1615 KNR2  1650 KNR3  1639 KNR4  1815 $$\n");
%!     fputs(fd, "                                  KNR5  1850 KNR6  1839\n");
%!     fputs(fd, "STRU BNR  3 ENR 443 TNR 11 MNR  1 KNR1  1650 KNR2  1638 KNR3  1639 KNR4  1850 $$\n");
%!     fputs(fd, "                                  KNR5  1838 KNR6  1839\n");
%!     fputs(fd, "STRU BNR  3 ENR 521 TNR 12 MNR  1 KNR1  1625 KNR2  1728 KNR3  1727 KNR4  1626 $$\n");
%!     fputs(fd, "                                  KNR5  1825 KNR6  1928 KNR7  1927 KNR8  1826\n");
%!     fputs(fd, "STRU BNR  3 ENR 522 TNR 12 MNR  1 KNR1  1625 KNR2  1730 KNR3  1729 KNR4  1728 $$\n");
%!     fputs(fd, "                                  KNR5  1825 KNR6  1930 KNR7  1929 KNR8  1928\n");
%!     fputs(fd, "STRU BNR  3 ENR 523 TNR 12 MNR  1 KNR1  1632 KNR2  1731 KNR3  1730 KNR4  1625 $$\n");
%!     fputs(fd, "                                  KNR5  1832 KNR6  1931 KNR7  1930 KNR8  1825\n");
%!     fputs(fd, "STRU BNR  3 ENR 524 TNR 12 MNR  1 KNR1  1626 KNR2  1727 KNR3  1734 KNR4  1633 $$\n");
%!     fputs(fd, "                                  KNR5  1826 KNR6  1927 KNR7  1934 KNR8  1833\n");
%!     fputs(fd, "STRU BNR  3 ENR 525 TNR 12 MNR  1 KNR1  1728 KNR2  1735 KNR3  1734 KNR4  1727 $$\n");
%!     fputs(fd, "                                  KNR5  1928 KNR6  1935 KNR7  1934 KNR8  1927\n");
%!     fputs(fd, "STRU BNR  3 ENR 526 TNR 12 MNR  1 KNR1  1729 KNR2  1736 KNR3  1735 KNR4  1728 $$\n");
%!     fputs(fd, "                                  KNR5  1929 KNR6  1936 KNR7  1935 KNR8  1928\n");
%!     fputs(fd, "STRU BNR  3 ENR 527 TNR 12 MNR  1 KNR1  1730 KNR2  1737 KNR3  1736 KNR4  1729 $$\n");
%!     fputs(fd, "                                  KNR5  1930 KNR6  1937 KNR7  1936 KNR8  1929\n");
%!     fputs(fd, "STRU BNR  3 ENR 528 TNR 11 MNR  1 KNR1  1731 KNR2  1737 KNR3  1730 KNR4  1931 $$\n");
%!     fputs(fd, "                                  KNR5  1937 KNR6  1930\n");
%!     fputs(fd, "STRU BNR  3 ENR 529 TNR 11 MNR  1 KNR1  1738 KNR2  1737 KNR3  1731 KNR4  1938 $$\n");
%!     fputs(fd, "                                  KNR5  1937 KNR6  1931\n");
%!     fputs(fd, "STRU BNR  3 ENR 530 TNR 12 MNR  1 KNR1  1639 KNR2  1738 KNR3  1731 KNR4  1632 $$\n");
%!     fputs(fd, "                                  KNR5  1839 KNR6  1938 KNR7  1931 KNR8  1832\n");
%!     fputs(fd, "STRU BNR  3 ENR 541 TNR 11 MNR  1 KNR1  1615 KNR2  1750 KNR3  1639 KNR4  1815 $$\n");
%!     fputs(fd, "                                  KNR5  1950 KNR6  1839\n");
%!     fputs(fd, "STRU BNR  3 ENR 543 TNR 11 MNR  1 KNR1  1750 KNR2  1738 KNR3  1639 KNR4  1950 $$\n");
%!     fputs(fd, "                                  KNR5  1938 KNR6  1839\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ BANDBREITENOPTIMIERUNG EINSCHALTEN\n");
%!     fputs(fd, "BAND OPTI 1\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ RANDBEDINGUNGEN\n");
%!     fputs(fd, "RBED KNR  1815 VX 1 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1825 VX 1 VY 2 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1826 VX 1 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1827 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1828 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1829 VX 0 VY 2 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1830 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1831 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1832 VX 1 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1833 VX 1 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1834 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1835 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1836 VX 0 VY 2 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1837 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1838 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1839 VX 1 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1850 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1927 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1928 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1929 VX 0 VY 2 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1930 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1931 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1934 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1935 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1936 VX 0 VY 2 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1937 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1938 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "RBED KNR  1950 VX 0 VY 0 VZ 3 DX 0 DY 0 DZ 0\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ KNOTENKRAEFTE FUER DIE LASTFAELLE 2 - 4\n");
%!     fputs(fd, "LAST LFNR  2 KNR     1 FX        648.386 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     2 FX        583.040 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     3 FX        461.814 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     4 FX        482.666 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     5 FX        409.096 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     6 FX        295.039 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     7 FX        281.992 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     8 FX        442.137 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR     9 FX        340.588 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    10 FX        340.588 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    11 FX        300.295 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    12 FX        193.124 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    13 FX        152.488 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    14 FX        292.916 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    15 FX        666.018 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    22 FX        160.317 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    39 FX        202.032 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR    50 FX        188.580 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   103 FX        461.814 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   104 FX        482.666 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   105 FX        409.096 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   106 FX        295.039 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   107 FX        281.992 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   110 FX        340.588 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   111 FX        300.295 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   112 FX        193.124 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   113 FX        152.488 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   114 FX        292.916 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   122 FX        160.317 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  2 KNR   150 FX        188.580 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     1 FX          0.000 FY        648.386 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     2 FX          0.000 FY        583.040 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     3 FX          0.000 FY        461.814 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     4 FX          0.000 FY        482.666 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     5 FX          0.000 FY        409.096 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     6 FX          0.000 FY        295.039 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     7 FX          0.000 FY        281.992 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     8 FX          0.000 FY        442.137 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR     9 FX          0.000 FY        340.588 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    10 FX          0.000 FY        340.588 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    11 FX          0.000 FY        300.295 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    12 FX          0.000 FY        193.124 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    13 FX          0.000 FY        152.488 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    14 FX          0.000 FY        292.916 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    15 FX          0.000 FY        666.018 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    22 FX          0.000 FY        160.317 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    39 FX          0.000 FY        202.032 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR    50 FX          0.000 FY        188.580 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   103 FX          0.000 FY        461.814 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   104 FX          0.000 FY        482.666 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   105 FX          0.000 FY        409.096 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   106 FX          0.000 FY        295.039 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   107 FX          0.000 FY        281.992 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   110 FX          0.000 FY        340.588 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   111 FX          0.000 FY        300.295 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   112 FX          0.000 FY        193.124 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   113 FX          0.000 FY        152.488 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   114 FX          0.000 FY        292.916 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   122 FX          0.000 FY        160.317 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  3 KNR   150 FX          0.000 FY        188.580 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  4 KNR     1 FX          0.000 FY          0.000 FZ        648.386\n");
%!     fputs(fd, "LAST LFNR  4 KNR     2 FX          0.000 FY          0.000 FZ        583.040\n");
%!     fputs(fd, "LAST LFNR  4 KNR     3 FX          0.000 FY          0.000 FZ        461.814\n");
%!     fputs(fd, "LAST LFNR  4 KNR     4 FX          0.000 FY          0.000 FZ        482.666\n");
%!     fputs(fd, "LAST LFNR  4 KNR     5 FX          0.000 FY          0.000 FZ        409.096\n");
%!     fputs(fd, "LAST LFNR  4 KNR     6 FX          0.000 FY          0.000 FZ        295.039\n");
%!     fputs(fd, "LAST LFNR  4 KNR     7 FX          0.000 FY          0.000 FZ        281.992\n");
%!     fputs(fd, "LAST LFNR  4 KNR     8 FX          0.000 FY          0.000 FZ        442.137\n");
%!     fputs(fd, "LAST LFNR  4 KNR     9 FX          0.000 FY          0.000 FZ        340.588\n");
%!     fputs(fd, "LAST LFNR  4 KNR    10 FX          0.000 FY          0.000 FZ        340.588\n");
%!     fputs(fd, "LAST LFNR  4 KNR    11 FX          0.000 FY          0.000 FZ        300.295\n");
%!     fputs(fd, "LAST LFNR  4 KNR    12 FX          0.000 FY          0.000 FZ        193.124\n");
%!     fputs(fd, "LAST LFNR  4 KNR    13 FX          0.000 FY          0.000 FZ        152.488\n");
%!     fputs(fd, "LAST LFNR  4 KNR    14 FX          0.000 FY          0.000 FZ        292.916\n");
%!     fputs(fd, "LAST LFNR  4 KNR    15 FX          0.000 FY          0.000 FZ        666.018\n");
%!     fputs(fd, "LAST LFNR  4 KNR    22 FX          0.000 FY          0.000 FZ        160.317\n");
%!     fputs(fd, "LAST LFNR  4 KNR    39 FX          0.000 FY          0.000 FZ        202.032\n");
%!     fputs(fd, "LAST LFNR  4 KNR    50 FX          0.000 FY          0.000 FZ        188.580\n");
%!     fputs(fd, "LAST LFNR  4 KNR   103 FX          0.000 FY          0.000 FZ        461.814\n");
%!     fputs(fd, "LAST LFNR  4 KNR   104 FX          0.000 FY          0.000 FZ        482.666\n");
%!     fputs(fd, "LAST LFNR  4 KNR   105 FX          0.000 FY          0.000 FZ        409.096\n");
%!     fputs(fd, "LAST LFNR  4 KNR   106 FX          0.000 FY          0.000 FZ        295.039\n");
%!     fputs(fd, "LAST LFNR  4 KNR   107 FX          0.000 FY          0.000 FZ        281.992\n");
%!     fputs(fd, "LAST LFNR  4 KNR   110 FX          0.000 FY          0.000 FZ        340.588\n");
%!     fputs(fd, "LAST LFNR  4 KNR   111 FX          0.000 FY          0.000 FZ        300.295\n");
%!     fputs(fd, "LAST LFNR  4 KNR   112 FX          0.000 FY          0.000 FZ        193.124\n");
%!     fputs(fd, "LAST LFNR  4 KNR   113 FX          0.000 FY          0.000 FZ        152.488\n");
%!     fputs(fd, "LAST LFNR  4 KNR   114 FX          0.000 FY          0.000 FZ        292.916\n");
%!     fputs(fd, "LAST LFNR  4 KNR   122 FX          0.000 FY          0.000 FZ        160.317\n");
%!     fputs(fd, "LAST LFNR  4 KNR   150 FX          0.000 FY          0.000 FZ        188.580\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "$ KNOTENKRAEFTE FUER DIE LASTFAELLE 5 - 7\n");
%!     fputs(fd, "LAST LFNR  5 KNR     2 FX          0.000 FY          0.000 FZ      -2666.533\n");
%!     fputs(fd, "LAST LFNR  5 KNR     3 FX          0.000 FY          0.000 FZ      -2112.106\n");
%!     fputs(fd, "LAST LFNR  5 KNR     4 FX          0.000 FY          0.000 FZ      -2207.471\n");
%!     fputs(fd, "LAST LFNR  5 KNR     5 FX          0.000 FY          0.000 FZ       1457.126\n");
%!     fputs(fd, "LAST LFNR  5 KNR     6 FX          0.000 FY          0.000 FZ       1050.874\n");
%!     fputs(fd, "LAST LFNR  5 KNR     7 FX          0.000 FY          0.000 FZ       1004.404\n");
%!     fputs(fd, "LAST LFNR  5 KNR     8 FX          0.000 FY          0.000 FZ       1574.811\n");
%!     fputs(fd, "LAST LFNR  5 KNR     9 FX          0.000 FY          0.000 FZ      -1557.678\n");
%!     fputs(fd, "LAST LFNR  5 KNR    10 FX          0.000 FY          0.000 FZ      -1557.678\n");
%!     fputs(fd, "LAST LFNR  5 KNR    11 FX          0.000 FY          0.000 FZ      -1373.400\n");
%!     fputs(fd, "LAST LFNR  5 KNR    12 FX          0.000 FY          0.000 FZ        687.873\n");
%!     fputs(fd, "LAST LFNR  5 KNR    13 FX          0.000 FY          0.000 FZ        543.135\n");
%!     fputs(fd, "LAST LFNR  5 KNR    14 FX          0.000 FY          0.000 FZ       1043.314\n");
%!     fputs(fd, "LAST LFNR  5 KNR    15 FX          0.000 FY          0.000 FZ       2372.234\n");
%!     fputs(fd, "LAST LFNR  5 KNR    22 FX          0.000 FY          0.000 FZ        571.021\n");
%!     fputs(fd, "LAST LFNR  5 KNR    39 FX          0.000 FY          0.000 FZ        719.602\n");
%!     fputs(fd, "LAST LFNR  5 KNR    50 FX          0.000 FY          0.000 FZ        671.689\n");
%!     fputs(fd, "LAST LFNR  5 KNR   103 FX          0.000 FY          0.000 FZ      -2112.106\n");
%!     fputs(fd, "LAST LFNR  5 KNR   104 FX          0.000 FY          0.000 FZ      -2207.471\n");
%!     fputs(fd, "LAST LFNR  5 KNR   105 FX          0.000 FY          0.000 FZ       1457.126\n");
%!     fputs(fd, "LAST LFNR  5 KNR   106 FX          0.000 FY          0.000 FZ       1050.874\n");
%!     fputs(fd, "LAST LFNR  5 KNR   107 FX          0.000 FY          0.000 FZ       1004.404\n");
%!     fputs(fd, "LAST LFNR  5 KNR   110 FX          0.000 FY          0.000 FZ      -1557.678\n");
%!     fputs(fd, "LAST LFNR  5 KNR   111 FX          0.000 FY          0.000 FZ      -1373.400\n");
%!     fputs(fd, "LAST LFNR  5 KNR   112 FX          0.000 FY          0.000 FZ        687.873\n");
%!     fputs(fd, "LAST LFNR  5 KNR   113 FX          0.000 FY          0.000 FZ        543.135\n");
%!     fputs(fd, "LAST LFNR  5 KNR   114 FX          0.000 FY          0.000 FZ       1043.314\n");
%!     fputs(fd, "LAST LFNR  5 KNR   122 FX          0.000 FY          0.000 FZ        571.021\n");
%!     fputs(fd, "LAST LFNR  5 KNR   150 FX          0.000 FY          0.000 FZ        671.689\n");
%!     fputs(fd, "LAST LFNR  6 KNR     3 FX          0.000 FY          0.000 FZ      -1804.210\n");
%!     fputs(fd, "LAST LFNR  6 KNR     4 FX          0.000 FY          0.000 FZ      -1885.674\n");
%!     fputs(fd, "LAST LFNR  6 KNR     5 FX          0.000 FY          0.000 FZ      -1670.597\n");
%!     fputs(fd, "LAST LFNR  6 KNR     6 FX          0.000 FY          0.000 FZ      -1204.829\n");
%!     fputs(fd, "LAST LFNR  6 KNR     7 FX          0.000 FY          0.000 FZ      -1151.550\n");
%!     fputs(fd, "LAST LFNR  6 KNR    10 FX          0.000 FY          0.000 FZ      -1330.605\n");
%!     fputs(fd, "LAST LFNR  6 KNR    11 FX          0.000 FY          0.000 FZ      -1173.190\n");
%!     fputs(fd, "LAST LFNR  6 KNR    12 FX          0.000 FY          0.000 FZ       -788.647\n");
%!     fputs(fd, "LAST LFNR  6 KNR    13 FX          0.000 FY          0.000 FZ       -622.705\n");
%!     fputs(fd, "LAST LFNR  6 KNR    14 FX          0.000 FY          0.000 FZ      -1196.161\n");
%!     fputs(fd, "LAST LFNR  6 KNR    22 FX          0.000 FY          0.000 FZ       -654.676\n");
%!     fputs(fd, "LAST LFNR  6 KNR    50 FX          0.000 FY          0.000 FZ       -770.092\n");
%!     fputs(fd, "LAST LFNR  6 KNR   103 FX          0.000 FY          0.000 FZ       1804.210\n");
%!     fputs(fd, "LAST LFNR  6 KNR   104 FX          0.000 FY          0.000 FZ       1885.674\n");
%!     fputs(fd, "LAST LFNR  6 KNR   105 FX          0.000 FY          0.000 FZ       1670.597\n");
%!     fputs(fd, "LAST LFNR  6 KNR   106 FX          0.000 FY          0.000 FZ       1204.829\n");
%!     fputs(fd, "LAST LFNR  6 KNR   107 FX          0.000 FY          0.000 FZ       1151.550\n");
%!     fputs(fd, "LAST LFNR  6 KNR   110 FX          0.000 FY          0.000 FZ       1330.605\n");
%!     fputs(fd, "LAST LFNR  6 KNR   111 FX          0.000 FY          0.000 FZ       1173.190\n");
%!     fputs(fd, "LAST LFNR  6 KNR   112 FX          0.000 FY          0.000 FZ        788.647\n");
%!     fputs(fd, "LAST LFNR  6 KNR   113 FX          0.000 FY          0.000 FZ        622.705\n");
%!     fputs(fd, "LAST LFNR  6 KNR   114 FX          0.000 FY          0.000 FZ       1196.161\n");
%!     fputs(fd, "LAST LFNR  6 KNR   122 FX          0.000 FY          0.000 FZ        654.676\n");
%!     fputs(fd, "LAST LFNR  6 KNR   150 FX          0.000 FY          0.000 FZ        770.092\n");
%!     fputs(fd, "LAST LFNR  7 KNR     2 FX       1453.259 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     3 FX        964.914 FY        614.739 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     4 FX        487.663 FY       1077.152 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     5 FX         -8.063 FY       1043.868 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     6 FX       -348.026 FY        663.801 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     7 FX       -452.404 FY        551.895 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     8 FX      -1105.445 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR     9 FX        848.934 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    10 FX        711.624 FY        453.370 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    11 FX        303.404 FY        670.161 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    12 FX         -7.250 FY        492.744 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    13 FX       -179.874 FY        343.080 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    14 FX       -469.930 FY        573.275 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    15 FX      -1665.200 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    22 FX       -113.073 FY        392.471 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    39 FX       -505.128 FY          0.000 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR    50 FX       -402.194 FY        251.131 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   103 FX        964.914 FY       -614.739 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   104 FX        487.663 FY      -1077.152 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   105 FX         -8.063 FY      -1043.868 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   106 FX       -348.026 FY       -663.801 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   107 FX       -452.404 FY       -551.895 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   110 FX        711.624 FY       -453.370 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   111 FX        303.404 FY       -670.161 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   112 FX         -7.250 FY       -492.744 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   113 FX       -179.874 FY       -343.080 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   114 FX       -469.930 FY       -573.275 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   122 FX       -113.073 FY       -392.471 FZ          0.000\n");
%!     fputs(fd, "LAST LFNR  7 KNR   150 FX       -402.194 FY       -251.131 FZ          0.000\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "LVAR LFNR 2 'load case 2'\n");
%!     fputs(fd, "LVAR LFNR 3 'load case 3'\n");
%!     fputs(fd, "LVAR LFNR 4 'load case 4'\n");
%!     fputs(fd, "LVAR LFNR 5 'load case 5'\n");
%!     fputs(fd, "LVAR LFNR 6 'load case 6'\n");
%!     fputs(fd, "LVAR LFNR 7 'load case 7'\n");
%!     fputs(fd, "$\n");
%!     fputs(fd, "ENDE\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [mesh, load_case] = fem_pre_mesh_import(filename, "eossp");
%!   for i=1:rows(mesh.elements.iso8)
%!     X = mesh.nodes(mesh.elements.iso8(i, :), 1:3);
%!     ## Half of the elements our mesh are inverted, which is not allowed for this package.
%!     if (all(X(:, 1) >= 0))
%!       mesh.elements.iso8(i, :) = mesh.elements.iso8(i, [4,3,2,1,8,7,6,5]);
%!     endif
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, load_case(1));
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_MASS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   num_modes = 10;
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, num_modes);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     unlink(filename);
%!   endif
%! end_unwind_protect

%!test
%! ### TEST18
%! close all;
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 30e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 2e-3;
%!   p = 25e6;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6}; Layers{Ceil(c/h)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
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
%!     group_defs(1).elem_type = "quad8";
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
%!     group_defs(2).elem_type = "quad8";
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
%!     group_defs(3).elem_type = "quad8";
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
%!     group_defs(4).elem_type = "quad8";
%!     groups = fem_pre_mesh_groups_create(mesh, group_defs, sqrt(eps));
%!     assert(numel(groups.quad8), 4);
%!     assert([groups.quad8.id], [group_defs.id]);
%!     assert(groups.quad8(1).nodes, mesh.groups.quad8(1).nodes);
%!     assert(groups.quad8(2).nodes, mesh.groups.quad8(1).nodes);
%!     assert(groups.quad8(3).nodes, mesh.groups.quad8(2).nodes);
%!     assert(groups.quad8(4).nodes, mesh.groups.quad8(2).nodes);
%!     assert(groups.quad8(1).elements, mesh.groups.quad8(1).elements);
%!     assert(groups.quad8(2).elements, mesh.groups.quad8(1).elements);
%!     assert(groups.quad8(3).elements, mesh.groups.quad8(2).elements);
%!     assert(groups.quad8(4).elements, mesh.groups.quad8(2).elements);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.quad8(find([[mesh.groups.quad8].id] == 1)).nodes, :) = true;
%!   load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(find([mesh.groups.quad8.id] == 2)).elements, :);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_SCA_TOT_MASS], load_case);
%!   assert(mtot, a * b * c * mesh.material_data.rho, sqrt(eps) * a * b * c * mesh.material_data.rho);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   sol_eig.stress = fem_ass_matrix(mesh, ...
%!                                   dof_map, ...
%!                                   [FEM_VEC_STRESS_CAUCH], ...
%!                                   load_case, ...
%!                                   sol_eig);
%!   X = mesh.nodes(unique(load_case.pressure.quad8.elements), 1:3).';
%!   dof_idx = dof_map.ndof(unique(load_case.pressure.quad8.elements), 1:3);
%!   f = sol_eig.f(:);
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
%!     fprintf(stderr, "mode %d f=%.0f\n", i, f(i));
%!   endfor
%!   assert(f, f_ref, tol * max(f_ref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ### TEST19
%! close all;
%! number_of_modes = 10;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! do_rotate = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 30e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   h = 3e-3;
%!   p = 25e6;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6}; Layers{Ceil(c/h)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   if (do_rotate)
%!     fputs(fd, "Rotate{{1, 0, 0},{0, 0, 0}, 30 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 1, 0},{0, 0, 0}, 7 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!     fputs(fd, "Rotate{{0, 0, 1},{0, 0, 0}, 15 * Pi / 180}{ Volume{tmp[1]}; }\n");
%!   endif
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
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
%!     group_defs(1).elem_type = "quad8";
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
%!     group_defs(2).elem_type = "quad8";
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
%!     group_defs(3).elem_type = "quad8";
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
%!     group_defs(4).elem_type = "quad8";
%!     groups = fem_pre_mesh_groups_create(mesh, group_defs, sqrt(eps));
%!     assert(numel(groups.quad8), 4);
%!     assert([groups.quad8.id], [group_defs.id]);
%!     assert(groups.quad8(1).nodes, mesh.groups.quad8(1).nodes);
%!     assert(groups.quad8(2).nodes, mesh.groups.quad8(1).nodes);
%!     assert(groups.quad8(3).nodes, mesh.groups.quad8(2).nodes);
%!     assert(groups.quad8(4).nodes, mesh.groups.quad8(2).nodes);
%!     assert(groups.quad8(1).elements, mesh.groups.quad8(1).elements);
%!     assert(groups.quad8(2).elements, mesh.groups.quad8(1).elements);
%!     assert(groups.quad8(3).elements, mesh.groups.quad8(2).elements);
%!     assert(groups.quad8(4).elements, mesh.groups.quad8(2).elements);
%!   endif
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.locked_dof(mesh.groups.quad8(find([[mesh.groups.quad8].id] == 1)).nodes, :) = true;
%!   load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(find([mesh.groups.quad8.id] == 2)).elements, :);
%!   Xp = mesh.nodes(load_case.pressure.quad8.elements, 1:3);
%!   xp = reshape(Xp(:, 1), rows(load_case.pressure.quad8.elements), columns(load_case.pressure.quad8.elements));
%!   yp = reshape(Xp(:, 2), rows(load_case.pressure.quad8.elements), columns(load_case.pressure.quad8.elements));
%!   zp = reshape(Xp(:, 3), rows(load_case.pressure.quad8.elements), columns(load_case.pressure.quad8.elements));
%!   load_case.pressure.quad8.p = p / 2 * (yp / b + zp / c); #repmat(p, rows(load_case.pressure.quad8.elements), columns(load_case.pressure.quad8.elements));
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M, ...
%!    mat_ass.Mlumped, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped, ...
%!    mtot] = fem_ass_matrix(mesh, ...
%!                           dof_map, ...
%!                           [FEM_MAT_STIFFNESS, ...
%!                            FEM_MAT_MASS, ...
%!                            FEM_MAT_MASS_LUMPED, ...
%!                            FEM_VEC_LOAD_CONSISTENT, ...
%!                            FEM_VEC_LOAD_LUMPED, ...
%!                            FEM_SCA_TOT_MASS], load_case);
%!   assert(mtot, a * b * c * mesh.material_data.rho, sqrt(eps) * a * b * c * mesh.material_data.rho);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, number_of_modes);
%!   [sol_eig_lumped] = fem_sol_modal(mesh, dof_map, setfield(mat_ass, "M", mat_ass.Mlumped), number_of_modes);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%!   X = mesh.nodes(unique(load_case.pressure.quad8.elements), 1:3).';
%!   dof_idx = dof_map.ndof(unique(load_case.pressure.quad8.elements), 1:3);
%!   F_con = full(mat_ass.R(dof_idx)).';
%!   F_lumped = full(mat_ass.Rlumped(dof_idx)).';
%!   M_con = cross(X, F_con);
%!   M_lumped = cross(X, F_lumped);
%!   Ftot_con = sum(F_con, 2);
%!   Mtot_con = sum(M_con, 2);
%!   Ftot_lumped = sum(F_lumped, 2);
%!   Mtot_lumped = sum(M_lumped, 2);
%!   Fx_an = -(b * c * p) / 2; ## grind(integrate(integrate(-1/2*p*(y/b+z/c),z,0,c),y,0,b));
%!   Mz_an = (7 * b^2 * c * p) / 24; ## grind(integrate(integrate(1/2*p*(y/b+z/c)*y,z,0,c),y,0,b));
%!   My_an = -(7 * b * c^2 * p) / 24; ## grind(integrate(integrate(-1/2*p*(y/b+z/c)*z,z,0,c),y,0,b));
%!   F_an = [Fx_an; 0; 0];
%!   M_an = [0; My_an; Mz_an];
%!   assert(Ftot_con, F_an, eps^0.9 * norm(F_an));
%!   assert(Ftot_lumped, F_an, eps^0.9 * norm(F_an));
%!   assert(Mtot_con, M_an, eps^0.9 * norm(M_an));
%!   assert(Mtot_lumped, M_an, 5e-3 * norm(M_an));
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
%!   assert(all(f_lumped <= f));
%!   assert(f, f_ref, tol * max(f_ref));

%!   figure("visible","off");
%!   hold on;
%!   fem_post_sol_plot(mesh);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title('undeformed mesh');
%!   opts_plot.elem_types = {"quad8", "iso20"};
%!   opts_plot.elem_groups.quad8 = [mesh.groups.quad8.id];
%!   opts_plot.elem_groups.iso20 = [mesh.groups.iso20.id];
%!   for i=1:min(number_of_modes, length(sol_eig.f))
%!     figure("visible", "off");
%!     hold on;
%!     fem_post_sol_plot(mesh, sol_eig, scale_eig/max(norm(sol_eig.def(:, :, i), "rows")), i, opts_plot);
%!     view(30,30);
%!     xlabel('x [m]');
%!     ylabel('y [m]');
%!     zlabel('z [m]');
%!     grid on;
%!     grid minor on;
%!     title(sprintf("%d. eigenmode: %gHz",i,sol_eig.f(i)));
%!   endfor

%!   figure("visible", "off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat, scale_eig/max(norm(sol_stat.def, "rows")), 1, opts_plot);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title("deformed mesh consistent load vector");

%!   figure("visible", "off");
%!   hold on;
%!   fem_post_sol_plot(mesh, sol_stat_lumped, scale_eig/max(norm(sol_stat_lumped.def, "rows")), 1, opts_plot);
%!   view(30,30);
%!   xlabel('x [m]');
%!   ylabel('y [m]');
%!   zlabel('z [m]');
%!   grid on;
%!   grid minor on;
%!   title("deformed mesh lumped load vector");

%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO 3
%! ## K.J.Bathe 2002, page 328 4.20a
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   animate = true;
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   mesh_size = 0.3;
%!   p1 = 0.006;
%!   E = 55;
%!   nu = 0.3;
%!   rho = 1000e-12;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%!   fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Circle(3) = {3,11,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Circle(7) = {7,12,8};\n");
%!   fputs(fd, "Line(8) = {8,9};\n");
%!   fputs(fd, "Line(9) = {9,10};\n");
%!   fputs(fd, "Line(10) = {10,1};\n");
%!   fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%!   fputs(fd, "Plane Surface(14) = {11};\n");
%!   fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; Layers{1}; Recombine; };\n", mesh_size);
%!   fputs(fd, "Recombine Surface{14, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%!   fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%!   fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%!   fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.quad8].id] == 2);
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp).nodes, 1:3) = true;
%!   grp_id_p1 = find([[mesh.groups.quad8].id] == 3);
%!   elem_id_p1 = mesh.groups.quad8(grp_id_p1).elements;
%!   elno_p1 = mesh.elements.quad8(elem_id_p1, :);

%!   load_case.pressure.quad8.elements = elno_p1;
%!   load_case.pressure.quad8.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_VEC_LOAD_LUMPED], ...
%!                                load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   grp_id_displacement = find([[mesh.groups.quad8].id] == 4);
%!   elem_id_displacement = mesh.groups.quad8(grp_id_displacement).elements;
%!   elno_id_displacement = mesh.elements.quad8(elem_id_displacement, :);
%!   delta = mean(sol_stat.def(elno_id_displacement, 3));
%!   grp_id_stress = find([[mesh.groups.quad8].id] == 5);
%!   elem_id_stress = mesh.groups.quad8(grp_id_stress).elements;
%!   elno_id_stress = mesh.elements.quad8(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(mesh.elements.iso20 == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(sol_stat.stress.taum.iso20(ridx(j), cidx(j), :), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   sigma1_max = 0;
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max = max(sigma1_max, max(eig(TAU)));
%!   endfor
%!   fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   ## K.J.Bathe page 329 4.20b
%!   sigma1_max_ref = 0.6056;
%!   delta_ref = -1.669;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   if (animate)
%!     opt_anim.scale_def = 10;
%!     opt_anim.animation_delay = 1;
%!     opt_anim.print_and_exit = true;
%!     opt_anim.print_to_file = filename;
%!     opt_anim.rotation_angle = [90, 0, 0] * pi / 180;
%!     opt_anim.skin_only = true;
%!     opt_anim.show_element = false;
%!     unwind_protect
%!       fem_post_sol_external(mesh, sol_stat, opt_anim);
%!       [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title("Gmsh - deformed mesh / continuous stress tensor");
%!     unwind_protect_cleanup
%!       unlink([opt_anim.print_to_file, "_001.jpg"]);
%!     end_unwind_protect
%!   endif
%!   assert(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%!   assert(delta, delta_ref, 0.04 * abs(delta_ref));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO 4
%! ## K.J.Bathe 2002, page 328 4.20a
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   animate = true;
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   mesh_size = 0.15;
%!   p1 = 0.006;
%!   E = 55;
%!   nu = 0.3;
%!   rho = 1000e-12;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fputs(fd, "Point(1) = { 0.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(2) = {10.0, 0.0,-20.0};\n");
%!   fputs(fd, "Point(3) = {10.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(4) = {15.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(5) = {65.0, 0.0, -5.0};\n");
%!   fputs(fd, "Point(6) = {65.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(7) = {15.0, 0.0,  5.0};\n");
%!   fputs(fd, "Point(8) = {10.0, 0.0, 10.0};\n");
%!   fputs(fd, "Point(9) = {10.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(10)= { 0.0, 0.0, 20.0};\n");
%!   fputs(fd, "Point(11)= {15.0, 0.0,-10.0};\n");
%!   fputs(fd, "Point(12)= {15.0, 0.0, 10.0};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Circle(3) = {3,11,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Circle(7) = {7,12,8};\n");
%!   fputs(fd, "Line(8) = {8,9};\n");
%!   fputs(fd, "Line(9) = {9,10};\n");
%!   fputs(fd, "Line(10) = {10,1};\n");
%!   fputs(fd, "Line Loop(11) = {1,2,3,4,5,6,7,8,9,10};\n");
%!   fputs(fd, "Plane Surface(14) = {11};\n");
%!   fprintf(fd, "tmp[] = Extrude {0, %g, 0}{ Surface{14}; Layers{1}; Recombine; };\n", mesh_size);
%!   fputs(fd, "Recombine Surface{14, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[11]};\n");
%!   fputs(fd, "Physical Surface(\"pressure\",3) = {tmp[7],tmp[8],tmp[9]};\n");
%!   fputs(fd, "Physical Surface(\"displacement\",4) = {tmp[6]};\n");
%!   fputs(fd, "Physical Surface(\"stress\",5) = {tmp[8]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 * mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.iso4].id] == 2);
%!   load_case.locked_dof(mesh.groups.iso4(grp_id_clamp).nodes, 1:3) = true;
%!   grp_id_p1 = find([[mesh.groups.iso4].id] == 3);
%!   elem_id_p1 = mesh.groups.iso4(grp_id_p1).elements;
%!   elno_p1 = mesh.elements.iso4(elem_id_p1, :);

%!   load_case.pressure.iso4.elements = elno_p1;
%!   load_case.pressure.iso4.p = [repmat(p1, rows(elno_p1), columns(elno_p1))];

%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.material_data.E = E;
%!   mesh.material_data.nu = nu;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_VEC_LOAD_LUMPED], ...
%!                                load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_stat.stress] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_VEC_STRESS_CAUCH], ...
%!                                      load_case, ...
%!                                      sol_stat);
%!   grp_id_displacement = find([[mesh.groups.iso4].id] == 4);
%!   elem_id_displacement = mesh.groups.iso4(grp_id_displacement).elements;
%!   elno_id_displacement = mesh.elements.iso4(elem_id_displacement, :);
%!   delta = mean(sol_stat.def(elno_id_displacement, 3));
%!   grp_id_stress = find([[mesh.groups.iso4].id] == 5);
%!   elem_id_stress = mesh.groups.iso4(grp_id_stress).elements;
%!   elno_id_stress = mesh.elements.iso4(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(mesh.elements.iso8 == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(sol_stat.stress.taum.iso8(ridx(j), cidx(j), :), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   sigma1_max = 0;
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max = max(sigma1_max, max(eig(TAU)));
%!   endfor
%!   fprintf(stderr, "mesh size=%.1f\n", mesh_size);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   ## K.J.Bathe page 329 4.20b
%!   sigma1_max_ref = 0.6056;
%!   delta_ref = -1.669;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   if (animate)
%!     opt_anim.scale_def = 10;
%!     opt_anim.animation_delay = 1;
%!     opt_anim.print_and_exit = true;
%!     opt_anim.print_to_file = filename;
%!     opt_anim.rotation_angle = [90, 0, 0] * pi / 180;
%!     opt_anim.skin_only = true;
%!     opt_anim.show_element = false;
%!     unwind_protect
%!       fem_post_sol_external(mesh, sol_stat, opt_anim);
%!       [img, map, alpha] = imread([opt_anim.print_to_file, "_001.jpg"]);
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title("Gmsh - deformed mesh / continuous stress tensor");
%!     unwind_protect_cleanup
%!       unlink([opt_anim.print_to_file, "_001.jpg"]);
%!     end_unwind_protect
%!   endif
%!   assert(sigma1_max, sigma1_max_ref, 0.02 * abs(sigma1_max_ref));
%!   assert(delta, delta_ref, 0.04 * abs(delta_ref));
%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 20
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   ri = 8e-3;
%!   ro = 10e-3;
%!   h = 12e-3;
%!   c = 2e-3;
%!   b = h - 2 * c;
%!   p1 = 25.79e6;
%!   p2 = 7.83e6;
%!   p3 = 1.3758e6;
%!   scale_def = 5e-3;
%!   mesh_size = 1e-3;
%!   enable_linear_dist = false;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "ri = %g;\n", ri);
%!   fprintf(fd, "ro = %g;\n", ro);
%!   fprintf(fd, "h = %g;\n", h);
%!   fprintf(fd, "c = %g;\n", c);
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!   fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!   fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!   fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!   fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!   fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!   fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!   fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,5};\n");
%!   fputs(fd, "Line(5) = {5,6};\n");
%!   fputs(fd, "Line(6) = {6,7};\n");
%!   fputs(fd, "Line(7) = {7,8};\n");
%!   fputs(fd, "Line(8) = {8,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fprintf(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{6}; Layers{Ceil(ro * Pi / 2 / %g)}; Recombine; };\n", mesh_size);
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!   fputs(fd, "Physical Surface(\"load3\",4) = {tmp[6]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp = find([[mesh.groups.quad8].id] == 1);
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp).nodes, :) = true;
%!   grp_id_p1 = find([[mesh.groups.quad8].id] == 3);
%!   grp_id_p2 = find([[mesh.groups.quad8].id] == 2);
%!   grp_id_p3 = find([[mesh.groups.quad8].id] == 4);
%!   elem_id_p1 = mesh.groups.quad8(grp_id_p1).elements;
%!   elem_id_p2 = mesh.groups.quad8(grp_id_p2).elements;
%!   elem_id_p3 = mesh.groups.quad8(grp_id_p3).elements;
%!   elno_p1 = mesh.elements.quad8(elem_id_p1, :);
%!   elno_p2 = mesh.elements.quad8(elem_id_p2, :);
%!   elno_p3 = mesh.elements.quad8(elem_id_p3, :);
%!   x1 = mesh.nodes(:, 1)(elno_p1);
%!   y1 = mesh.nodes(:, 2)(elno_p1);
%!   z1 = mesh.nodes(:, 3)(elno_p1);
%!   Phi1 = atan2(y1, x1);

%!   x2 = mesh.nodes(:, 1)(elno_p2);
%!   y2 = mesh.nodes(:, 2)(elno_p2);
%!   z2 = mesh.nodes(:, 3)(elno_p2);
%!   Phi2 = atan2(y2, x2);

%!   load_case.pressure.quad8.elements = [elno_p1; elno_p2; elno_p3];
%!   load_case.pressure.quad8.p = [p1 * Phi1 / (pi / 2) .* z1 / h;
%!                                 p2 * Phi2 / (pi / 2) .* z2 / h;
%!                                 repmat(p3, rows(elno_p3), columns(elno_p3))];
%!   if (enable_linear_dist)
%!     p_mid = load_case.pressure.quad8.p(:, 1:3);
%!     load_case.pressure.quad8.p(:, 4:6) = [0.5 * (p_mid(:, 1) + p_mid(:, 2)), 0.5 * (p_mid(:, 2) + p_mid(:, 3)), 0.5 * (p_mid(:, 1) + p_mid(:, 3))];
%!   endif
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.Rlumped] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS, ...
%!                                       FEM_VEC_LOAD_CONSISTENT, ...
%!                                       FEM_VEC_LOAD_LUMPED], ...
%!                                      load_case);

%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat_lumped = fem_sol_static(mesh, dof_map, setfield(mat_ass, "R", mat_ass.Rlumped));
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat_lumped, scale_def / max(norm(sol_stat_lumped.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - lumped pressure load");

%!   X1 = mesh.nodes(unique(elno_p1), 1:3).';
%!   X2 = mesh.nodes(unique(elno_p2), 1:3).';
%!   X3 = mesh.nodes(unique(elno_p3), 1:3).';
%!   dof1 = dof_map.ndof(unique(elno_p1), 1:3);
%!   dof2 = dof_map.ndof(unique(elno_p2), 1:3);
%!   dof3 = dof_map.ndof(unique(elno_p3), 1:3);
%!   F1_con = full(mat_ass.R(dof1)).';
%!   F2_con = full(mat_ass.R(dof2)).';
%!   F3_con = full(mat_ass.R(dof3)).';
%!   M1_con = cross(X1, F1_con);
%!   M2_con = cross(X2, F2_con);
%!   M3_con = cross(X3, F3_con);
%!   F1_lumped = full(mat_ass.Rlumped(dof1)).';
%!   F2_lumped = full(mat_ass.Rlumped(dof2)).';
%!   F3_lumped = full(mat_ass.Rlumped(dof3)).';
%!   M1_lumped = cross(X1, F1_lumped);
%!   M2_lumped = cross(X2, F2_lumped);
%!   M3_lumped = cross(X3, F3_lumped);
%!   Ftot1_con = sum(F1_con, 2);
%!   Ftot2_con = sum(F2_con, 2);
%!   Ftot3_con = sum(F3_con, 2);
%!   Mtot1_con = sum(M1_con, 2);
%!   Mtot2_con = sum(M2_con, 2);
%!   Mtot3_con = sum(M3_con, 2);

%!   Ftot1_lumped = sum(F1_lumped, 2);
%!   Ftot2_lumped = sum(F2_lumped, 2);
%!   Ftot3_lumped = sum(F3_lumped, 2);
%!   Mtot1_lumped = sum(M1_lumped, 2);
%!   Mtot2_lumped = sum(M2_lumped, 2);
%!   Mtot3_lumped = sum(M3_lumped, 2);

%!   F1_an = [(2*(pi/2-1)*((h^2-2*c*h+c^2)/2-c^2/2)*p1*ri)/(pi*h);
%!            (2*((h^2-2*c*h+c^2)/2-c^2/2)*p1*ri)/(pi*h);
%!            0];

%!   F2_an = [-(2*(pi/2-1)*((h^2-2*c*h+c^2)/2-c^2/2)*p2*ro)/(pi*h);
%!            -(2*((h^2-2*c*h+c^2)/2-c^2/2)*p2*ro)/(pi*h);
%!            0];

%!   M1_an = [-(2*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p1*ri)/(pi*h);
%!            (2*(pi/2-1)*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p1*ri)/(pi*h);
%!            0];

%!   M2_an = [(2*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p2*ro)/(pi*h);
%!            -(2*(pi/2-1)*((h^3-3*c*h^2+3*c^2*h-c^3)/3-c^3/3)*p2*ro)/(pi*h);
%!            0];

%!   ys = 2 / 3 * (ro^3 - ri^3) * sin(pi/2) / ((ro^2 - ri^2) * pi / 2);

%!   F3_an = [0;
%!            0;
%!            -p3 * (ro^2 - ri^2) * pi / 4];

%!   M3_an = [-ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            ys * p3 * (ro^2 - ri^2) * pi / 4;
%!            0];

%!   assert(Ftot1_con, F1_an, 1e-4 * norm(F1_an));
%!   assert(Ftot2_con, F2_an, 1e-4 * norm(F2_an));
%!   assert(Ftot1_lumped, F1_an, 2e-3 * norm(F1_an));
%!   assert(Ftot2_lumped, F2_an, 2e-3 * norm(F2_an));

%!   assert(Mtot1_con, M1_an, 1e-4 * norm(M1_an));
%!   assert(Mtot2_con, M2_an, 1e-4 * norm(M2_an));
%!   assert(Mtot1_lumped, M1_an, 5e-3 * norm(M1_an));
%!   assert(Mtot2_lumped, M2_an, 5e-3 * norm(M2_an));

%!   assert(Ftot3_con, F3_an, eps^0.2 * norm(F3_an));
%!   assert(Ftot3_lumped, F3_an, eps^0.2* norm(F3_an));
%!   assert(Mtot3_con, M3_an, eps^0.2 * norm(M3_an));
%!   assert(Mtot3_lumped, M3_an, eps^0.2 * norm(M3_an));

%!   fprintf(stderr, "err(F1_con)=%e\n", norm(Ftot1_con - F1_an) / norm(F1_an));
%!   fprintf(stderr, "err(F2_con)=%e\n", norm(Ftot2_con - F2_an) / norm(F2_an));
%!   fprintf(stderr, "err(F3_con)=%e\n", norm(Ftot3_con - F3_an) / norm(F3_an));
%!   fprintf(stderr, "err(M1_con)=%e\n", norm(Mtot1_con - M1_an) / norm(M1_an));
%!   fprintf(stderr, "err(M2_con)=%e\n", norm(Mtot2_con - M2_an) / norm(M2_an));
%!   fprintf(stderr, "err(M3_con)=%e\n", norm(Mtot3_con - M3_an) / norm(M3_an));
%!   fprintf(stderr, "err(F1_lumped)=%e\n", norm(Ftot1_lumped - F1_an) / norm(F1_an));
%!   fprintf(stderr, "err(F2_lumped)=%e\n", norm(Ftot2_lumped - F2_an) / norm(F2_an));
%!   fprintf(stderr, "err(F3_lumped)=%e\n", norm(Ftot3_lumped - F3_an) / norm(F3_an));
%!   fprintf(stderr, "err(M1_lumped)=%e\n", norm(Mtot1_lumped - M1_an) / norm(M1_an));
%!   fprintf(stderr, "err(M2_lumped)=%e\n", norm(Mtot2_lumped - M2_an) / norm(M2_an));
%!   fprintf(stderr, "err(M3_lumped)=%e\n", norm(Mtot3_lumped - M3_an) / norm(M3_an));

%!   figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect


%!test
%! ## TEST21
%! close all;
%! number_of_modes = 3;
%! scale_eig = 10e-3;
%! tol = 1e-2;
%! plot_modes = true;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 30e-3;
%!   b = 20e-3;
%!   c = 10e-3;
%!   d = -5e-3;
%!   e = 35e-3;
%!   h = 2e-3;
%!   alpha = 1e-6;
%!   beta = 1e-4;
%!   fprintf(fd, "a=%g;\n", a);
%!   fprintf(fd, "b=%g;\n", b);
%!   fprintf(fd, "c=%g;\n", c);
%!   fprintf(fd, "h = %g;\n", h);
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!   fputs(fd, "Point(2) = {a,0.0,0.0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0.0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!   fputs(fd, "Line(1) = {4,3};\n");
%!   fputs(fd, "Line(2) = {3,2};\n");
%!   fputs(fd, "Line(3) = {2,1};\n");
%!   fputs(fd, "Line(4) = {1,4};\n");
%!   fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!   fputs(fd, "  Surface{6}; Layers{Ceil(c/h)}; Recombine;\n");
%!   fputs(fd, "};\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"modal\",1) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"interfaces\",2) = {tmp[2]};\n");
%!   fputs(fd, "SetFactory(\"Built-in\");\n");
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
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   cms_opt.modes.number = int32(6);
%!   cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!   cms_opt.nodes.interfaces.number = int32(rows(mesh.nodes) + 2);
%!   cms_opt.invariants = false;
%!   cms_opt.algorithm = "shift-invert";
%!   mesh.nodes(cms_opt.nodes.modal.number, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%!   mesh.nodes([cms_opt.nodes.interfaces.number], :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                    [1, 2], ...
%!                                                    [cms_opt.nodes.modal.number, ...
%!                                                     cms_opt.nodes.interfaces.number], "quad8");
%!   [mesh, mat_ass_cms, dof_map_cms, sol_eig_cms] = fem_cms_create(mesh, load_case, cms_opt);
%!   mat_ass_cms.Dred = alpha * mat_ass_cms.Mred + beta * mat_ass_cms.Kred;
%!   
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
%!       opt_plot.elem_types = {"quad8"};
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
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 22: patch test of iso20
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 8e-3;
%!   b = 15e-3;
%!   c = 12e-3;
%!   px = 25.79e6;
%!   py = 7.83e6;
%!   pz = 1.3758e6;
%!   mesh_size = 7e-3;
%!   scale_def = 10e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a = %g;\n", a);
%!   fprintf(fd, "b = %g;\n", b);
%!   fprintf(fd, "c = %g;\n", c);
%!   fprintf(fd, "h = %g;\n", mesh_size);
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = {0,0,0,h};\n");
%!   fputs(fd, "Point(2) = {a,0,0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; Layers{Ceil(c / h)}; Recombine; };\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!   fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!   fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.quad8].id] == 4);
%!   grp_id_clamp_y = find([[mesh.groups.quad8].id] == 1);
%!   grp_id_clamp_z = find([[mesh.groups.quad8].id] == 5);
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp_x).nodes, 1) = true;
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp_y).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.quad8(grp_id_clamp_z).nodes, 3) = true;
%!   grp_id_px = find([[mesh.groups.quad8].id] == 2);
%!   grp_id_py = find([[mesh.groups.quad8].id] == 3);
%!   grp_id_pz = find([[mesh.groups.quad8].id] == 6);
%!   elem_id_px = mesh.groups.quad8(grp_id_px).elements;
%!   elem_id_py = mesh.groups.quad8(grp_id_py).elements;
%!   elem_id_pz = mesh.groups.quad8(grp_id_pz).elements;
%!   elno_px = mesh.elements.quad8(elem_id_px, :);
%!   elno_py = mesh.elements.quad8(elem_id_py, :);
%!   elno_pz = mesh.elements.quad8(elem_id_pz, :);
%!   load_case.pressure.quad8.elements = [elno_px; elno_py; elno_pz];
%!   load_case.pressure.quad8.p = [repmat(px, rows(elno_px), columns(elno_px));
%!                                 repmat(py, rows(elno_py), columns(elno_py));
%!                                 repmat(pz, rows(elno_pz), columns(elno_pz))];
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
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
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure_list();
%!   tol = eps^0.7;
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso20(:, :, 1) / -px - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso20(:, :, 2) / -py - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso20(:, :, 3) / -pz - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso20(:, :, 4:6) / max([px,py,pz]))))) < tol);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 23: patch test of iso8
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 8e-3;
%!   b = 15e-3;
%!   c = 12e-3;
%!   px = 25.79e6;
%!   py = 7.83e6;
%!   pz = 1.3758e6;
%!   mesh_size = 7e-3;
%!   scale_def = 10e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a = %g;\n", a);
%!   fprintf(fd, "b = %g;\n", b);
%!   fprintf(fd, "c = %g;\n", c);
%!   fprintf(fd, "h = %g;\n", mesh_size);
%!   fputs(fd, "Point(1) = {0,0,0,h};\n");
%!   fputs(fd, "Point(2) = {a,0,0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; Layers{Ceil(c / h)}; Recombine; };\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!   fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!   fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.iso4].id] == 4);
%!   grp_id_clamp_y = find([[mesh.groups.iso4].id] == 1);
%!   grp_id_clamp_z = find([[mesh.groups.iso4].id] == 5);
%!   load_case.locked_dof(mesh.groups.iso4(grp_id_clamp_x).nodes, 1) = true;
%!   load_case.locked_dof(mesh.groups.iso4(grp_id_clamp_y).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.iso4(grp_id_clamp_z).nodes, 3) = true;
%!   grp_id_px = find([[mesh.groups.iso4].id] == 2);
%!   grp_id_py = find([[mesh.groups.iso4].id] == 3);
%!   grp_id_pz = find([[mesh.groups.iso4].id] == 6);
%!   elem_id_px = mesh.groups.iso4(grp_id_px).elements;
%!   elem_id_py = mesh.groups.iso4(grp_id_py).elements;
%!   elem_id_pz = mesh.groups.iso4(grp_id_pz).elements;
%!   elno_px = mesh.elements.iso4(elem_id_px, :);
%!   elno_py = mesh.elements.iso4(elem_id_py, :);
%!   elno_pz = mesh.elements.iso4(elem_id_pz, :);
%!   load_case.pressure.iso4.elements = [elno_px; elno_py; elno_pz];
%!   load_case.pressure.iso4.p = [repmat(px, rows(elno_px), columns(elno_px));
%!                                 repmat(py, rows(elno_py), columns(elno_py));
%!                                 repmat(pz, rows(elno_pz), columns(elno_pz))];
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
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
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure_list();
%!   tol = eps^0.7;
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso8(:, :, 1) / -px - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso8(:, :, 2) / -py - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso8(:, :, 3) / -pz - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.iso8(:, :, 4:6) / max([px,py,pz]))))) < tol);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 24: patch test of tet10
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 8e-3;
%!   b = 15e-3;
%!   c = 12e-3;
%!   px = 25.79e6;
%!   py = 7.83e6;
%!   pz = 1.3758e6;
%!   mesh_size = 7e-3;
%!   scale_def = 10e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a = %g;\n", a);
%!   fprintf(fd, "b = %g;\n", b);
%!   fprintf(fd, "c = %g;\n", c);
%!   fprintf(fd, "h = %g;\n", mesh_size);
%!   fputs(fd, "Point(1) = {0,0,0,h};\n");
%!   fputs(fd, "Point(2) = {a,0,0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; };\n");
%!   fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!   fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!   fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   grp_id_clamp_x = find([[mesh.groups.tria6].id] == 4);
%!   grp_id_clamp_y = find([[mesh.groups.tria6].id] == 1);
%!   grp_id_clamp_z = find([[mesh.groups.tria6].id] == 5);
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp_x).nodes, 1) = true;
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp_y).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.tria6(grp_id_clamp_z).nodes, 3) = true;
%!   grp_id_px = find([[mesh.groups.tria6].id] == 2);
%!   grp_id_py = find([[mesh.groups.tria6].id] == 3);
%!   grp_id_pz = find([[mesh.groups.tria6].id] == 6);
%!   elem_id_px = mesh.groups.tria6(grp_id_px).elements;
%!   elem_id_py = mesh.groups.tria6(grp_id_py).elements;
%!   elem_id_pz = mesh.groups.tria6(grp_id_pz).elements;
%!   elno_px = mesh.elements.tria6(elem_id_px, :);
%!   elno_py = mesh.elements.tria6(elem_id_py, :);
%!   elno_pz = mesh.elements.tria6(elem_id_pz, :);
%!   load_case.pressure.tria6.elements = [elno_px; elno_py; elno_pz];
%!   load_case.pressure.tria6.p = [repmat(px, rows(elno_px), columns(elno_px));
%!                                 repmat(py, rows(elno_py), columns(elno_py));
%!                                 repmat(pz, rows(elno_pz), columns(elno_pz))];
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
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
%!   figure("visible", "off");
%!   fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   title("static deflection - consistent pressure load");
%!   figure_list();
%!   tol = eps^0.7;
%!   assert(max(max(max(abs(sol_stat.stress.tau.tet10(:, :, 1) / -px - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.tet10(:, :, 2) / -py - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.tet10(:, :, 3) / -pz - 1)))) < tol);
%!   assert(max(max(max(abs(sol_stat.stress.tau.tet10(:, :, 4:6) / max([px,py,pz]))))) < tol);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 25: patch test of iso20
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 8e-3;
%!   b = 15e-3;
%!   c = 12e-3;
%!   shift = 1e-4;
%!   mesh_size = 7e-3;
%!   scale_def = 1e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a = %g;\n", a);
%!   fprintf(fd, "b = %g;\n", b);
%!   fprintf(fd, "c = %g;\n", c);
%!   fprintf(fd, "h = %g;\n", mesh_size);
%!   fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   fputs(fd, "Point(1) = {0,0,0,h};\n");
%!   fputs(fd, "Point(2) = {a,0,0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; Layers{Ceil(c / h)}; Recombine; };\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!   fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!   fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_MASS], ...
%!                                load_case);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 10, shift);
%!   sol_eig.def *= scale_def / max(max(max(abs(sol_eig.def))));
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig);
%!   tolt = eps^0.5;
%!   assert(max(max(max(max(abs(sol_eig.stress.tau.iso20(:, :, :, 1:6)))))) < tolt * max(max(max(max(abs(sol_eig.stress.tau.iso20(:, :, :, 7:end)))))));
%!   tolf = eps^0.4;
%!   assert(all(sol_eig.f(1:6) < tolf * max(sol_eig.f(7:10))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 26: patch test of iso8
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 8e-3;
%!   b = 15e-3;
%!   c = 12e-3;
%!   shift = 1e-4;
%!   mesh_size = 7e-3;
%!   scale_def = 1e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a = %g;\n", a);
%!   fprintf(fd, "b = %g;\n", b);
%!   fprintf(fd, "c = %g;\n", c);
%!   fprintf(fd, "h = %g;\n", mesh_size);
%!   fputs(fd, "Point(1) = {0,0,0,h};\n");
%!   fputs(fd, "Point(2) = {a,0,0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; Layers{Ceil(c / h)}; Recombine; };\n");
%!   fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!   fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!   fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!   fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_MASS], ...
%!                                load_case);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 10, shift);
%!   sol_eig.def *= scale_def / max(max(max(abs(sol_eig.def))));
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig);
%!   tolt = eps^0.5;
%!   assert(max(max(max(max(abs(sol_eig.stress.tau.iso8(:, :, :, 1:6)))))) < tolt * max(max(max(max(abs(sol_eig.stress.tau.iso8(:, :, :, 7:end)))))));
%!   tolf = eps^0.4;
%!   assert(all(sol_eig.f(1:6) < tolf * max(sol_eig.f(7:10))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 27: patch test of tet10
%! close all;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!   [fd, msg] = fopen([filename, ".geo"], "wt");
%!   if (fd == -1)
%!     error("failed to open file \"%s.geo\"", filename);
%!   endif
%!   a = 8e-3;
%!   b = 15e-3;
%!   c = 12e-3;
%!   shift = 1e-4;
%!   mesh_size = 7e-3;
%!   scale_def = 1e-3;
%!   fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!   fprintf(fd, "a = %g;\n", a);
%!   fprintf(fd, "b = %g;\n", b);
%!   fprintf(fd, "c = %g;\n", c);
%!   fprintf(fd, "h = %g;\n", mesh_size);
%!   fputs(fd, "Point(1) = {0,0,0,h};\n");
%!   fputs(fd, "Point(2) = {a,0,0,h};\n");
%!   fputs(fd, "Point(3) = {a,b,0,h};\n");
%!   fputs(fd, "Point(4) = {0,b,0,h};\n");
%!   fputs(fd, "Line(1) = {1,2};\n");
%!   fputs(fd, "Line(2) = {2,3};\n");
%!   fputs(fd, "Line(3) = {3,4};\n");
%!   fputs(fd, "Line(4) = {4,1};\n");
%!   fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!   fputs(fd, "Plane Surface(6) = {5};\n");
%!   fputs(fd, "tmp[] = Extrude {0, 0, c}{ Surface{6}; };\n");
%!   fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-y\", 1) = {tmp[2]};\n");
%!   fputs(fd, "Physical Surface(\"load-x\", 2) = {tmp[3]};\n");
%!   fputs(fd, "Physical Surface(\"load-y\", 3) = {tmp[4]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-x\", 4) = {tmp[5]};\n");
%!   fputs(fd, "Physical Surface(\"clamp-z\", 5) = {6};\n");
%!   fputs(fd, "Physical Surface(\"load-z\", 6) = {tmp[0]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!   unlink([filename, ".msh"]);
%!   load_case.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.material_data.E = 210000e6;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 7850;
%!   mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.M] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_MASS], ...
%!                                load_case);
%!   sol_eig = fem_sol_modal(mesh, dof_map, mat_ass, 10, shift);
%!   sol_eig.def *= scale_def / max(max(max(abs(sol_eig.def))));
%!   [sol_eig.stress] = fem_ass_matrix(mesh, ...
%!                                     dof_map, ...
%!                                     [FEM_VEC_STRESS_CAUCH], ...
%!                                     load_case, ...
%!                                     sol_eig);
%!   tolf = eps^0.4;
%!   tolt = eps^0.5;
%!   assert(max(max(max(max(abs(sol_eig.stress.tau.tet10(:, :, :, 1:6)))))) < tolt * max(max(max(max(abs(sol_eig.stress.tau.tet10(:, :, :, 7:end)))))));
%!   assert(all(sol_eig.f(1:6) < tolf * max(sol_eig.f(7:10))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
