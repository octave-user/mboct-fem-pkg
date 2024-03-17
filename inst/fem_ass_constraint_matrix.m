## Copyright (C) 2019(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{C1}, @var{C2}, @var{mesh_constr}, @var{dof_map_constr}, @var{mat_ass_constr}] = fem_ass_constraint_matrix(@var{mesh}, @var{elem_type_master}, @var{grp_id_master}, @var{elem_type_slave}, @var{grp_id_slave}, @var{maxdist}, @var{constrtype})
## Assemble a sparse matrix containing the constraint equations of a single group of surface to node constraints
## @end deftypefn

function [C1, C2, mesh_constr, dof_map_constr, mat_ass_constr, joints] = fem_ass_constraint_matrix(mesh, elem_type_master, grp_id_master, elem_type_slave, grp_id_slave, maxdist, constrtype)
  grp_type_master = getfield(mesh.groups, elem_type_master);
  grp_type_slave = getfield(mesh.groups, elem_type_slave);

  grp_idx_master = find([grp_type_master.id] == grp_id_master);
  grp_idx_slave = find([grp_type_slave.id] == grp_id_slave);

  grp_master = getfield(mesh.groups, elem_type_master)(grp_idx_master);
  grp_slave = getfield(mesh.groups, elem_type_slave)(grp_idx_slave);

  elem_master = getfield(mesh.elements, elem_type_master)(grp_master.elements, :);
  elem_slave = getfield(mesh.elements, elem_type_master)(grp_slave.elements, :);

  node_idx_master = grp_master.nodes(:);
  node_idx_slave = grp_slave.nodes(:);

  node_idx_trans = zeros(rows(mesh.nodes), 1, "int32");

  node_idx_trans(node_idx_master) = 1:rows(node_idx_master);
  node_idx_trans(node_idx_slave) = rows(node_idx_master) + (1:rows(node_idx_slave));



  mesh_constr.nodes = mesh.nodes([node_idx_master; node_idx_slave], :);
  mesh_constr.elements = struct();

  sfncon.slave = node_idx_trans(node_idx_slave);
  sfncon.master = node_idx_trans(elem_master);
  sfncon.maxdist = maxdist;
  sfncon.constraint = constrtype;
  elem_type_map = {"tria6", "sfncon6";
                   "tria6h", "sfncon6h";
                   "iso4", "sfncon4";
                   "quad8", "sfncon8";
                   "quad9", "sfncon9";
                   "tria10", "sfncon10"};

  elem_type_sfncon = [];

  for i=1:rows(elem_type_map)
    switch (elem_type_master)
      case elem_type_map{i, 1}
        elem_type_sfncon = elem_type_map{i, 2};
        break
    endswitch
  endfor

  if (isempty(elem_type_sfncon))
    error("element type \"%s\" is not yet supported", elem_type_master);
  endif

  mesh_constr.elements = setfield(mesh_constr.elements, elem_type_sfncon, sfncon);

  mesh_constr.materials = struct();
  mesh_constr.material_data = struct()([]);

  dof_map_constr = [];
  mat_ass_constr = [];
  
  ## dof_in_use = load_case_constr.locked_dof = false(rows(mesh_constr.nodes), 6);
  ## dof_in_use(:, 1:3) = true;
  ## dof_map_constr = fem_ass_dof_map(mesh_constr, load_case_constr, dof_in_use);

  joints = fem_pre_mesh_constr_surf_to_node(mesh_constr.nodes, mesh_constr.elements, FEM_DO_STRUCTURAL);

  ## [mat_ass_constr.K, ...
  ##  mat_ass_constr.mat_info, ...
  ##  mat_ass_constr.mesh_info] = fem_ass_matrix(mesh_constr, dof_map_constr, FEM_MAT_STIFFNESS, load_case_constr);

  ## C1 = mat_ass_constr.K(dof_map_constr.edof.joints, dof_map_constr.ndof(1:rows(node_idx_master), 1:3)(:));

  ## if (nargout > 1)
  ##   C2 = mat_ass_constr.K(dof_map_constr.edof.joints, dof_map_constr.ndof((rows(node_idx_master) + 1):end, 1:3)(:));
  ## endif

  ## dof_map_inv = zeros(rows(mesh.nodes), 1, "int32");

  ## for j=1:3
  ##   dof_map_inv(dof_map_constr.ndof(:, j)) = 1:rows(dof_map_constr.ndof);
  ## endfor

  ## idx_node_dof = find(dof_map_inv);

  ##  assert(dof_map_inv(dof_map_constr.ndof(:, 1:3)),repmat(int32(1:rows(dof_map_constr.ndof))(:),1,3));

  ## eq_idx_slave = zeros(rows(node_idx_slave), 3, "int32");

  ## [ridx, cidx] = find(C2);

  ## idx_node_slave = dof_map_inv(cidx + dof_map_constr.ndof(rows(node_idx_master), 3)) - rows(node_idx_master);

  C1 = zeros(rows(node_idx_slave), rows(node_idx_master) * 3);

  for i=1:numel(joints)
    for j=2:columns(joints(i).nodes)
      C1(i, (joints(i).nodes(j) - 1) * 3 + (1:3)) = joints(i).C(1, (j - 1) * 6 + (1:3));
    endfor
  endfor

  C2 = zeros(rows(node_idx_slave), rows(node_idx_slave) * 3);

  for i=1:numel(joints)
    C2(i, (joints(i).nodes(1) - rows(node_idx_master) - 1) * 3 + (1:3)) = joints(i).C(1, 1:3);
  endfor
  ## [ridx, cidx] = find(C1);

  ## idx_node_slave = dof_map_inv(cidx);

  ## idx_dof_last_master_no = dof_map_constr.ndof(rows(node_idx_master), 3);

  ## for i=1:rows(node_idx_slave)
  ##   for j=1:3
  ##     idx = find(cidx == dof_map_constr.ndof(i + rows(node_idx_master), j) - idx_dof_last_master_no);
  ##     if (isempty(idx))
  ##       continue;
  ##     endif
  ##     eq_idx_slave(i, j) = ridx(idx);
  ##   endfor
  ## endfor
  ## C1 * U1 + C2 * U2 = 0
  ## C2 * U2 = -C1 * U1
  ## U2n = -C1 * U1
endfunction

%!test
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
%!     a = 10e-3;
%!     b = 12e-3;
%!     c = 14e-3;
%!     h = 2e-3;
%!     F1 = 1e10;
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
%!     fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",2) = {tmp[2]};\n");
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
%!   mesh1 = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh1.materials.tet10 = ones(rows(mesh1.elements.tet10), 1, "int32");
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh1.material_data.rho = 7850;
%!   mesh1.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   data(1).mesh = mesh1;
%!   data(2).mesh = mesh1;
%!   data(2).mesh.nodes(:, 1) += a;
%!   for j=1:numel(data(2).mesh.groups.tria6)
%!     data(2).mesh.groups.tria6(j).id += 100;
%!   endfor
%!   mesh = fem_post_mesh_merge(data);
%!   maxdist = 1e-8 * max(abs([a, b, c]));
%!   constrtype = FEM_CT_SLIDING;
%!   [C1, ...
%!    C2, ...
%!    mesh_constr, ...
%!    dof_map_constr, ...
%!    mat_ass_constr, ...
%!    joints] = fem_ass_constraint_matrix(mesh, ...
%!                                                "tria6", 2, ...
%!                                                "tria6", 101, ...
%!                                                maxdist, constrtype);
%!   mesh.elements.sfncon6.slave = mesh.groups.tria6(3).nodes(:);
%!   mesh.elements.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(2).elements, :);
%!   mesh.elements.sfncon6.maxdist = maxdist;
%!   mesh.elements.sfncon6.constraint = FEM_CT_FIXED;
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(1).elements) = 1;
%!   mesh.materials.tet10(mesh.groups.tet10(2).elements) = 2;
%!   mesh.nodes(end + 1, 1:3) = [3 * a, 0.5 * b, 0.5 * c];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, mesh.groups.tria6(4).id, rows(mesh.nodes), "tria6");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(1).nodes, 1:3) = true;
%!   load_case.loaded_nodes = int32(rows(mesh.nodes));
%!   load_case.loads = [0, 0, F1, 0, 0, 0];
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   [sol_stat, U] = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   U1 = U(dof_map.ndof(mesh.groups.tria6(2).nodes, 1:3).'(:));
%!   U2 = U(dof_map.ndof(mesh.groups.tria6(3).nodes, 1:3).'(:));
%!   f = norm(C1 * U1 + C2 * U2) / max(norm(C1 * U1), norm(C2 * U1));
%!   tol = eps^0.9;
%!   assert_simple(f < tol);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
