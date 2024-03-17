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

function [C1, C2, mesh_constr] = fem_ass_constraint_matrix(mesh, elem_type_master, grp_id_master, elem_type_slave, grp_id_slave, maxdist)
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
  sfncon.constraint = FEM_CT_SLIDING;
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

  mesh_constr.joints = fem_pre_mesh_constr_surf_to_node(mesh_constr.nodes, mesh_constr.elements, FEM_DO_STRUCTURAL);

  nnz_C2 = nnz_C1 = int32(0);

  for i=1:numel(mesh_constr.joints)
    nnz_C1 += (columns(mesh_constr.joints(i).C) - 6) / 2;
    nnz_C2 += 3;
  endfor

  ridx_C1 = zeros(nnz_C1, 1, "int32");
  ridx_C2 = zeros(nnz_C2, 1, "int32");
  cidx_C1 = zeros(nnz_C1, 1, "int32");
  cidx_C2 = zeros(nnz_C2, 1, "int32");
  data_C1 = zeros(nnz_C1, 1);
  data_C2 = zeros(nnz_C2, 1);

  idx_C2 = idx_C1 = int32(0);

  nnz_C2 = nnz_C1 = int32(0);

  for i=1:numel(mesh_constr.joints)
    for j=2:columns(mesh_constr.joints(i).nodes)
      ridx_C1(nnz_C1 + (1:3)) = i;
      cidx_C1(nnz_C1 + (1:3)) = (mesh_constr.joints(i).nodes(j) - 1) * 3 + (1:3);
      data_C1(nnz_C1 + (1:3)) = mesh_constr.joints(i).C(1, (j - 1) * 6 + (1:3));
      nnz_C1 += 3;
    endfor
    ridx_C2(nnz_C2 + (1:3)) = i;
    cidx_C2(nnz_C2 + (1:3)) = (mesh_constr.joints(i).nodes(1) - rows(node_idx_master) - 1) * 3 + (1:3);
    data_C2(nnz_C2 + (1:3)) = mesh_constr.joints(i).C(1, 1:3);
    nnz_C2 += 3;
  endfor

  C1 = sparse(ridx_C1, cidx_C1, data_C1, rows(node_idx_slave), rows(node_idx_master) * 3);
  C2 = sparse(ridx_C2, cidx_C2, data_C2, rows(node_idx_slave), rows(node_idx_slave) * 3);

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
%!     b = 10e-3;
%!     c = 10e-3;
%!     h = 10e-3;
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
%!   [C1, ...
%!    C2, ...
%!    mesh_constr] = fem_ass_constraint_matrix(mesh, ...
%!                                                "tria6", 2, ...
%!                                                "tria6", 101, ...
%!                                                maxdist);
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
%!   f = norm(C1 * U1 + C2 * U2) / (norm(C1 * U1) + norm(C2 * U2));
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
