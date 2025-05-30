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
  eltypes = fem_pre_mesh_elem_type();

  gmsh_elem_types = {eltypes.name};

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

  options.elem_types = {options.elem_types{valid_elem_type}};

  fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

  inumelem = int32(0);
  elem_types = fieldnames(mesh.elements);

  for i=1:numel(elem_types)
    switch (elem_types{i})
      case "beam2"
        elname = "line2";
      case "beam3"
        elname = "line3";
      otherwise
        elname = elem_types{i};
    endswitch

    switch (elname)
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
    switch (group_types{i})
      case options.elem_types
      otherwise
        continue;
    endswitch

    groups = getfield(mesh.groups, group_types{i});

    idx_elem_type = fem_pre_mesh_elem_type_index({eltypes.name}, group_types{i});

    if (isempty(idx_elem_type))
      continue;
    endif

    group_dim = eltypes(idx_elem_type).dim;

    for j=1:numel(groups)
      grp_name = groups(j).name;
      grp_name(isspace(grp_name)) = "_"; ## spaces in group names are not allowed for gmsh format
      fprintf(fd, "%d %d \"%s\"\n", group_dim, groups(j).id, grp_name);
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
      case "beam2"
        elname = "line2";
      case "beam3"
        elname = "line3";
      otherwise
        elname = elem_types{i};
    endswitch

    idx_elem_type = fem_pre_mesh_elem_type_index({eltypes.name}, elname);

    if (isempty(idx_elem_type))
      continue;
    endif

    if (~isscalar(idx_elem_type))
      error("element type \"%s\" is not unique", elem_types{i});
    endif

    elem_norder = eltypes(idx_elem_type).nordernonp;

    if (isempty(elem_norder))
      elem_norder = eltypes(idx_elem_type).norder;
    endif

    elem_node_order = zeros(1, numel(elem_norder));
    elem_node_order(elem_norder) = 1:numel(elem_node_order);
    elem_type_id = eltypes(idx_elem_type).id;

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
