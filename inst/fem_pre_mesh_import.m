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
## @deftypefn {Function File} [@var{mesh}, @var{load_case}] = fem_pre_mesh_import(@var{filename})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_import(@dots{}, @var{format}, @var{options})
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
    options.promote_elem = {"tet4", "tria3", "penta6", "pyra5"};
  endif

  if (~isfield(options, "elem_type"))
    options.elem_type = {"tet10", "tria6", "tet20", "tria10", "tet4", "penta6", "tria3", "iso8", "iso4", "iso20", "iso27", "quad8", "quad9", "penta15", "pyra5"};
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
    [fd, msg] = fopen(filename, "r");

    if (fd == -1)
      error("failed to open file \"%s\": %s", filename, msg);
    endif

    nodes = zeros(0, 4);
    elements = zeros(0, 0, "int32");
    p_name = {};
    p_dim = [];
    p_id = [];

    while (true)
      line = fgetl(fd);

      if (~ischar(line))
        break
      endif

      switch (line)
        case "$MeshFormat"
          [version_number, file_type, data_size, count, msg] = fscanf(fd, "%f %d %d\n", "C");

          if (count ~= 3)
            error("invalid header in Gmsh file");
          endif

          switch (version_number)
            case 2.2
            otherwise
              error("invalid file format: only Gmsh file format version 2 is supported");
          endswitch

          switch (file_type)
            case 0
            otherwise
              error("invalid file format: only Gmsh ASCII format is supported");
          endswitch
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
          elements = zeros(num_elements, 20, "int32");

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

          eltype = fem_pre_mesh_elem_type();

          use_elem_type = false(1, numel(eltype));
          for i=1:numel(eltype)
            switch (eltype(i).name)
              case options.elem_type
                use_elem_type(i) = true;
            endswitch
          endfor

          idx_elem = find(use_elem_type);
          idx_promote(idx_elem) = 1:numel(idx_elem);

          eltype = eltype(use_elem_type);

          for i=1:numel(eltype)
            if (eltype(i).promote > 0)
              eltype(i).promote = idx_promote(eltype(i).promote);
            endif
          endfor

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
    [fd, msg] = fopen(filename, "r");

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
    [fd, msg] = fopen(filename, "r");

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    nodes = zeros(4, 0);
    node_id = zeros(1, 0, "int32");
    material = zeros(4, 0);
    material_id = zeros(1, 0, "int32");

    mesh.nodes = zeros(0, 6);
    mesh.material_data = struct("E", [], "nu", [], "rho", [])([]);

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

            if (~isfield(mesh, "elements"))
              mesh.elements = struct();
            endif

            if (~isfield(mesh, "materials"))
              mesh.materials = struct();
            endif

            switch (elem1(3))
              case 8
                if (~isfield(mesh.materials, "iso4"))
                  mesh.materials.iso4 = [];
                endif
                if (~isfield(mesh.elements, "iso4"))
                  mesh.elements.iso4 = [];
                endif
              case {11, 12}
                if (~isfield(mesh.materials, "iso8"))
                  mesh.materials.iso8 = [];
                endif
                if (~isfield(mesh.elements, "iso8"))
                  mesh.elements.iso8 = [];
                endif
            endswitch

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

