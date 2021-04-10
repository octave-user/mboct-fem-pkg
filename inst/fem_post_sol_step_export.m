## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_post_sol_step_export(@var{filename}, @var{sol}, @var{idx_sol}, @var{idx_t}, @var{t})
## @deftypefnx {} fem_post_sol_step_export(@dots{}, @var{scale})
## Export stresses and displacements of a single time step or mode shape to a file in Gmsh format.
##
## @var{filename} @dots{} Output filename
##
## @var{sol} @dots{} Finite element solution data structure.
##
## @var{idx_sol} @dots{} Physical index in the solution array to be exported. This is the array index for Octave.
##
## @var{idx_t} @dots{} Logical index of the time step to be exported. This is the index Gmsh will use.
##
## @var{t} @dots{} Time value related to @var{idx_t}. This is the value Gmsh will print to screen for this time step.
##
## @var{scale} @dots{} Multiply displacements by @var{scale} before exporting it. The actual display scale may be changed in Gmsh.
##
## @seealso{fem_post_sol_export, fem_post_sol_external}
## @end deftypefn

function fem_post_sol_step_export(filename, sol, idx_sol, idx_t, t, scale)
  if (nargin < 5 || nargin > 6 || nargout > 0)
    print_usage();
  endif

  if (nargin < 6)
    scale = 1;
  endif

  fd = -1;
  
  unwind_protect
    [fd, msg] = fopen(filename, "w");

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    nodal_fields = struct("name", {"def", "theta", "p"}, ...
                          "title", {"nodal deformation", "nodal temperature", "acoustic pressure"}, ...
                          "dim", {3, 1, 1}, ...
                          "value", {@(x) x(:, 1:3, idx_sol), @(x) x(:, idx_sol), @(x) x(:, idx_sol)});

    for i=1:numel(nodal_fields)
      if (isfield(sol, nodal_fields(i).name))
        x = nodal_fields(i).value(getfield(sol, nodal_fields(i).name));
        fprintf(fd, "$NodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n%d\n%d\n", nodal_fields(i).title, t, idx_t - 1, nodal_fields(i).dim, size(x, 1));
        fprintf(fd, "%d %e\n", [1:size(x, 1); x.']);
        fputs(fd, "$EndNodeData\n");
      endif
    endfor
    
    field_type = {"stress", "strain"};

    idxtens = int32([1, 4, 6, 4, 2, 5, 6, 5, 3]);
    stress_type = {"discontinuous stress tensor", "continuous stress tensor", "van Mises stress", "discontinuous strain", "continuous strain"};
    stress_field = {"tau", "taum", "vmis", "epsilon", "epsilonm"};
    stress_comp = int32([9, 9, 1, 9, 9]);
    
    for n=1:numel(field_type)
      if (isfield(sol, field_type{n}))        
        for l=1:numel(stress_field)
          if (isfield(getfield(sol, field_type{n}), stress_field{l}))
            elem_stress = {"iso4", "quad8", "iso8", "iso20", "tet10", "penta15", "tet10h"};

            inumelem = int32(0);
            
            for i=1:numel(elem_stress)
              if (isfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i}))
                tau = getfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i});
                inumelem += rows(tau);
              endif
            endfor

            fprintf(fd, "$ElementNodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n%d\n%d\n", stress_type{l}, t, idx_t - 1, stress_comp(l), inumelem);

            inumelem = int32(0);

            for i=1:numel(elem_stress)
              if (isfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i}))
                switch (elem_stress{i})
                  case "iso3"
                    idxnode = int32([1:3]);
                  case "iso4"
                    idxnode = int32([1:4]);
		  case "quad8"
		    idxnode = int32([1:8]);
                  case "iso8"
                    idxnode = int32([5:8, 1:4]);
		  case "iso20"
		    idxnode = int32([5:8, 1:4, 17, 19, 20, 18, 9, 12, 14, 10, 11, 13, 15, 16]);
		  case "penta15"
		    idxnode =  int32([1, 2, 3, 4, 5, 6, 7, 10, 8, 13, 15, 14, 9, 11, 12]);
                  case {"tet10", "tet10h"}
                    idxnode = int32([1:8, 10, 9]);
                endswitch
                
                tau = getfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i});

                tauout = zeros(2 + numel(idxnode) * stress_comp(l), rows(tau));
                tauout(1, :) = inumelem + (1:rows(tau));
                tauout(2, :) = numel(idxnode);
                
                for k=1:numel(idxnode)
                  if (stress_comp(l) == 1)
                    taukl = tau(:, idxnode(k), idx_sol);
                  else
                    taukl = tau(:, idxnode(k), idxtens, idx_sol);
                  endif
                  
                  tauout(2 + (k - 1) * stress_comp(l) + (1:stress_comp(l)), :) = reshape(taukl, rows(tau), stress_comp(l)).';
                endfor

                inumelem += rows(tau);
                
                format = ["%d %d", repmat(" %e", 1, numel(idxnode) * stress_comp(l)), "\n"];

                fprintf(fd, format, tauout);
              endif
            endfor

            fputs(fd, "$EndElementNodeData\n");
          endif
        endfor
      endif
    endfor
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

%!demo
%! close all;
%! E1 = 70000e6;
%! nu1 = 0.3;
%! rho1 = 2700;
%! E2 = 60000e6;
%! nu2 = 0.3;
%! rho2 = 5000;
%! F1 = 100;
%! h = 5e-3;
%! param.a = 50e-3;
%! param.b = 20e-3;
%! param.c = 10e-3;
%! param.o = 25e-3;
%! param.g = 0.5e-3;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     geo_file = [filename, ".geo"];
%!     fd = fopen(geo_file, "w");
%!     if (fd == -1)
%!       error("failed to open file %s", geo_file);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0.0, 0.0, 0.0};\n");
%!     fputs(fd, "Point(2) = {  a, 0.0, 0.0};\n");
%!     fputs(fd, "Point(3) = {  a,   b, 0.0};\n");
%!     fputs(fd, "Point(4) = {0.0,   b, 0.0};\n");
%!     fputs(fd, "Point(5) = {    o, 0.0, c + g};\n");
%!     fputs(fd, "Point(6) = {o + a, 0.0, c + g};\n");
%!     fputs(fd, "Point(7) = {o + a,   b, c + g};\n");
%!     fputs(fd, "Point(8) = {    o,   b, c + g};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line(5) = {5,6};\n");
%!     fputs(fd, "Line(6) = {6,7};\n");
%!     fputs(fd, "Line(7) = {7,8};\n");
%!     fputs(fd, "Line(8) = {8,5};\n");
%!     fputs(fd, "Line Loop(9) = {1,2,3,4};\n");
%!     fputs(fd, "Line Loop(10) = {5,6,7,8};\n");
%!     fputs(fd, "Plane Surface(11) = {9};\n");
%!     fputs(fd, "Plane Surface(12) = {10};\n");
%!     fputs(fd, "v1[] = Extrude {0,0.0,c} {\n");
%!     fputs(fd, "  Surface{11};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "v2[] = Extrude {0,0.0,c} {\n");
%!     fputs(fd, "  Surface{12};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"v1\") = {v1[1]};\n");
%!     fputs(fd, "Physical Volume(\"v2\") = {v2[1]};\n");
%!     fputs(fd, "Physical Surface(\"s1\") = {v1[0]};\n");
%!     fputs(fd, "Physical Surface(\"s2\") = {12};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v1[5]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v2[3]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.element_size = h;
%!   opt.mesh.jacobian_range = [0.5, 1.5];
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, param, opt);
%!   mesh.material_data(1).rho = rho1;
%!   mesh.material_data(1).C = fem_pre_mat_isotropic(E1, nu1);
%!   mesh.material_data(2).rho = rho2;
%!   mesh.material_data(2).C = fem_pre_mat_isotropic(E2, nu2);
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(1).elements) = 1;
%!   mesh.materials.tet10(mesh.groups.tet10(2).elements) = 2;
%!   mesh.nodes(end + 1, 1:3) = [param.o + param.a, 0.5 * param.b, 1.5 * param.c + param.g];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, mesh.groups.tria6(4).id, rows(mesh.nodes), "tria6");
%!   load_case = fem_pre_load_case_create_empty(2);
%!   load_case(1).locked_dof = false(size(mesh.nodes));
%!   load_case(1).locked_dof(mesh.groups.tria6(3).nodes, 1:3) = true;
%!   load_case(1).loaded_nodes = int32(rows(mesh.nodes));
%!   load_case(1).loads = [0, 0, F1, 0, 0, 0];
%!   load_case(2).loaded_nodes = int32(rows(mesh.nodes));
%!   load_case(2).loads = [F1, 0, 0, 0, 0, 0];
%!   elem.sfncon6.slave = mesh.groups.tria6(1).nodes(:);
%!   elem.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(2).elements, :);
%!   elem.sfncon6.maxdist = param.g * (1 + sqrt(eps));
%!   mesh.elements.joints = fem_pre_mesh_constr_surf_to_node(mesh.nodes, elem);
%!   dof_map = fem_ass_dof_map(mesh, load_case(1));
%!   [mat_ass.K, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   for i=1:size(sol_stat.def, 3)
%!     fem_post_sol_step_export(sprintf("%s_post_pro_%02d.msh", filename, i),  sol_stat, i, i, i);
%!   endfor
%!   opts.format = "gmsh";
%!   fem_post_mesh_export([filename,  "_post_pro_00.msh"], mesh, opts);
%!   fn = dir([filename, "_post_pro_*.msh"]);
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([filename, "_post_pro.geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [filename, "_post_pro.geo"]);
%!     endif
%!     for i=1:numel(fn)
%!       fprintf(fd, "Merge \"%s\";\n", fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!     fputs(fd, "Mesh.SurfaceEdges = 0;\n");
%!     fputs(fd, "Mesh.VolumeFaces = 0;\n");
%!     fputs(fd, "Mesh.SurfaceFaces = 0;\n");
%!     fputs(fd, "Mesh.VolumeEdges = 0;\n");
%!     fputs(fd, "General.Trackball = 0;\n");
%!     fputs(fd, "General.RotationX = -90;\n");
%!     fputs(fd, "View[0].Visible = 1;\n");
%!     fputs(fd, "View[1].Visible = 0;\n");
%!     fputs(fd, "View[0].VectorType = 5;\n");
%!     fputs(fd, "View[0].DisplacementFactor = 500;\n");
%!     fputs(fd, "View[0].ExternalView = 1;\n");
%!     fputs(fd, "View[0].ShowElement = 1;\n");
%!     fputs(fd, "View[0].IntervalsType = 3;\n");
%!     fputs(fd, "View[0].ShowTime = 1;\n");
%!     for i=1:size(sol_stat.def, 3)
%!       fprintf(fd, "View[0].TimeStep = %d;\n", i - 1);
%!       fprintf(fd, "Print \"%s_post_pro_%02d.jpg\";\n", filename, i);
%!     endfor
%!     fputs(fd, "Exit;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   status = spawn_wait(spawn("gmsh", {[filename, "_post_pro.geo"]}));
%!   if (status ~= 0)
%!     error("gmsh failed with status %d", status);
%!   endif
%!   fn = dir([filename, "_post_pro_*.jpg"]);
%!   for i=1:numel(fn)
%!     [img, alpha, map] = imread(fullfile(fn(i).folder, fn(i).name));
%!     figure("visible", "off");
%!     imshow(img, map);
%!     title(sprintf("deformed mesh load case %d", i));
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
