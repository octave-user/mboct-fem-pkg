## Copyright Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{rbe3}] = fem_pre_mesh_rbe3_from_surf(@var{mesh}, @var{group_id}, @var{master_node_idx})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_rbe3_from_surf(@dots{}, @var{elem_type})
##
## Builds rbe3 elements from specified groups of tria6 elements and assigns weighting factors proportional to equivalent surface area
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{group_id} @dots{} array of group id-numbers used to identify surface elements with slave nodes
##
## @var{master_node_idx} @dots{} array of master node indices for rbe3 elements
##
## @var{elem_type} @dots{} the element type addressed by @var{group_id} (e.g. "tria6" or "iso4")
##
## @end deftypefn

function rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, group_id, master_node_idx, elem_type)
  if (nargin < 3 || nargin > 4 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    elem_type = "tria6";
  endif

  empty_cell = cell(1, numel(group_id));
  
  rbe3 = struct("nodes", empty_cell, "weight", empty_cell, "elements", empty_cell);

  if (numel(group_id) ~= numel(master_node_idx))
    error("numel(group_id) does not match numel(master_node_idx)");
  endif
  
  group_idx = cell(1, numel(group_id));

  msh_groups = getfield(mesh.groups, elem_type);
  
  for i=1:numel(group_id)
    group_idx_i = find([msh_groups.id] == group_id(i));

    if (numel(group_idx_i) == 0)
      error("no elements found with id %d", group_id(i));
    endif

    group_idx{i} = group_idx_i;
  endfor
  
  for i=1:numel(group_idx)
    inode = [msh_groups(group_idx{i}).nodes];
    inode_tr = zeros(rows(mesh.nodes), 1, "int32");
    inode_tr(inode) = 1:numel(inode);
    ielno = getfield(mesh.elements, elem_type)([msh_groups(group_idx{i}).elements], :);
    press_elem.elements = inode_tr(ielno);
    press_elem.p = ones(rows(press_elem.elements), columns(press_elem.elements));
    load_case_i.pressure = setfield(struct(), elem_type, press_elem);
    dof_map_i.ndof = [reshape(1:(numel(inode) * 3), numel(inode), 3), zeros(numel(inode), 3)];
    dof_map_i.totdof = numel(inode) * 3;
    mesh_i.nodes = mesh.nodes(inode, :);
    mesh_i.elements = struct();
    mesh_i.materials = struct();
    mesh_i.material_data = struct("C", [], "rho", [])([]);
    mat_ass_i.R = fem_ass_matrix(mesh_i, ...
                                 dof_map_i, ...
                                 [FEM_VEC_LOAD_CONSISTENT], ...
                                 load_case_i);

    F = norm(mat_ass_i.R(dof_map_i.ndof(:, 1:3)), "rows");
    rbe3(i).nodes = int32([master_node_idx(i), inode]);
    rbe3(i).weight = F.';
    rbe3(i).elements = setfield(struct(), elem_type, [msh_groups(group_idx{i}).elements]);
    
    if (any(rbe3(i).weight < 0))
      error("negative weight factor detected");
    endif
  endfor
endfunction

%!test ##demo
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!  filename(filename == "\\") = "/";
%! endif
%! fd = -1;
%! unwind_protect
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
%! fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
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
%! fputs(fd, "  Surface{6}; Layers{Ceil(c / h)}; Recombine;\n");
%! fputs(fd, "};\n");
%! fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"modal\",1) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"interfaces\",2) = {tmp[2]};\n");
%! unwind_protect_cleanup
%! if (fd ~= -1)
%! fclose(fd);
%! endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", "-optimize_ho", "-ho_min", "0.95", "-ho_max", "1.05", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! node_id_ground = rows(mesh.nodes) + 1;
%! node_id_load = rows(mesh.nodes) + 2;
%! mesh.nodes(node_id_ground, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%! mesh.nodes(node_id_load, :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];

%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);

%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                [1, 2], ...
%!                                                [node_id_ground, ...
%!                                                 node_id_load], "iso4");
%! mesh.elements.joints(1).nodes = node_id_ground;
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.joints(2).nodes = node_id_load;
%! mesh.elements.joints(2).C = [zeros(3, 3), eye(3)];
%! load_case.joints(1).U = zeros(6, 1);
%! load_case.joints(2).U = zeros(3, 1);
%! load_case.loaded_nodes = int32(node_id_load);
%! load_case.loads = [1000, 100, 100, 0, 0, 0];
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! figure("visible", "off");
%! opts.elem_types = {"iso8"};
%! fem_post_sol_plot(mesh, sol_stat, 20e-3 / max(norm(sol_stat.def(:, 1:3), "rows")), 1, opts);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("deformed mesh");
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test ##demo
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!  filename(filename == "\\") = "/";
%! endif
%! fd = -1;
%! unwind_protect
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
%! fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
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
%! unwind_protect_cleanup
%! if (fd ~= -1)
%! fclose(fd);
%! endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", "-optimize_ho", "-ho_min", "0.95", "-ho_max", "1.05", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! node_id_ground = rows(mesh.nodes) + 1;
%! node_id_load = rows(mesh.nodes) + 2;
%! mesh.nodes(node_id_ground, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%! mesh.nodes(node_id_load, :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];

%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);

%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                [1, 2], ...
%!                                                [node_id_ground, ...
%!                                                 node_id_load], "iso4");
%! mesh.elements.joints(1).nodes = node_id_ground;
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.joints(2).nodes = node_id_load;
%! mesh.elements.joints(2).C = [0, 0, 0, 0, 1, 0;];
%! load_case.joints(1).U = zeros(6, 1);
%! load_case.joints(2).U = [1];
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! figure("visible", "off");
%! opts.elem_types = {"iso8"};
%! fem_post_sol_plot(mesh, sol_stat, 20e-3 / max(norm(sol_stat.def(:, 1:3), "rows")), 1, opts);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("deformed mesh");
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect

%!test ##demo
%! close all;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!  filename(filename == "\\") = "/";
%! endif
%! fd = -1;
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! a = 30e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 0e-3;
%! e = 30e-3;
%! h = 4e-3;
%! fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
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
%! unwind_protect_cleanup
%! if (fd ~= -1)
%! fclose(fd);
%! endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-ho_min", "0.95", "-ho_max", "1.05", [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  warning("gmsh failed with status %d", status);
%! endif
%! unlink([filename, ".geo"]);
%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! node_id_ground = rows(mesh.nodes) + 1;
%! node_id_load = rows(mesh.nodes) + 2;
%! mesh.nodes(node_id_ground, :) = [d, 0.5 * b, 0.5 * c, 0, 0, 0];
%! mesh.nodes(node_id_load, :) = [e, 0.5 * b, 0.5 * c, 0, 0, 0];

%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);

%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, ...
%!                                                [1, 2], ...
%!                                                [node_id_ground, ...
%!                                                 node_id_load], "tria6");
%! mesh.elements.joints(1).nodes = node_id_ground;
%! mesh.elements.joints(1).C = eye(6);
%! mesh.elements.joints(2).nodes = node_id_load;
%! mesh.elements.joints(2).C = eye(6);
%! load_case.joints(1).U = zeros(6, 1);
%! load_case.joints(2).U = [0; 0; 0; 30 * pi / 180; 0; 0];
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, ...
%!  mat_ass.R] = fem_ass_matrix(mesh, ...
%!                              dof_map, ...
%!                              [FEM_MAT_STIFFNESS, ...
%!                               FEM_VEC_LOAD_CONSISTENT], ...
%!                              load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%! figure("visible", "off");
%! opts.elem_types = {"tet10"};
%! scale_def = 5e-3;
%! fem_post_sol_plot(mesh, sol_stat, scale_def / max(norm(sol_stat.def(:, 1:3), "rows")), 1, opts);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("deformed mesh");
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
