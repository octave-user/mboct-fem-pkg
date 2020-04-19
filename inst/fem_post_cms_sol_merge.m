## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{sol} = fem_post_cms_sol_merge(@var{mesh}, @var{dof_map}, @var{sol_tot})
##
## Merge the solutions (displacements, stresses) of several modal bodies
##
## @var{mesh} @dots{} the combined mesh returned from fem_post_mesh_merge
##
## @var{dof_map} @dots{} the combined degree of freedom mapping returned from fem_post_mesh_merge
##
## @var{sol_tot} @dots{} nodal solutions of several flexible bodies returned from fem_post_cms_expand
## @seealso{fem_post_cms_sol_import}
## @end deftypefn

function [sol] = fem_post_cms_sol_merge(mesh, dof_map, sol_tot)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  num_nodes = int32(0);
  num_steps = int32(-1);
  num_dof = int32(-1);

  for i=1:numel(sol_tot.bodies)
    num_nodes += size(sol_tot.bodies(i).def, 1);

    if (num_dof == -1)
      num_dof = size(sol_tot.bodies(i).def, 2);
    elseif (num_dof ~= size(sol_tot.bodies(i).def, 2))
      error("invalid number of dof detected");
    endif

    if (num_steps == -1)
      num_steps = size(sol_tot.bodies(i).def, 3);
    elseif (num_steps ~= size(sol_tot.bodies(i).def, 3))
      error("invalid number of load steps detected");
    endif
  endfor

  sol.def = zeros(num_nodes, num_dof, num_steps);

  for i=1:numel(sol_tot.bodies)
    nnodes = size(sol_tot.bodies(i).def, 1);
    sol.def(dof_map.submesh.offset.nodes(i) + (1:nnodes), :, :) = sol_tot.bodies(i).def;
  endfor

  if (isfield(sol_tot.bodies(i), "stress"))
    elem_type = {"iso4", "iso8", "tet10"};
    stress_type = {"tau", "taum", "vmis"};

    sol.stress = struct();
    
    for k=1:numel(stress_type)
      for j=1:numel(elem_type)
        num_elem = int32(0);
        num_nodes = int32(-1);
        num_comp = int32(-1);
        num_steps = int32(-1);

        stress_elem_type_m = struct();
        
        for i=1:numel(sol_tot.bodies)
          if (isfield(sol_tot.bodies(i).stress, stress_type{k}))
            stress_type_k = getfield(sol_tot.bodies(i).stress, stress_type{k});

            if (isfield(stress_type_k, elem_type{j}))
              stress_elem = getfield(stress_type_k, elem_type{j});
              num_elem += size(stress_elem, 1);

              if (num_nodes == -1)
                num_nodes = size(stress_elem, 2);
              elseif (num_nodes ~= size(stress_elem, 2))
                error("invalid number of element nodes detected");
              endif

              if (num_comp == -1)
                num_comp = size(stress_elem, 3);
              elseif (num_comp ~= size(stress_elem, 3))
                error("invalid number of stress components detected");
              endif

              if (num_steps == -1)
                num_steps = size(stress_elem, 4);
              elseif (num_steps ~= size(stress_elem, 4))
                error("invalid number of load steps detected");
              endif
            endif
          endif
        endfor

        if (num_elem > 0)
          stress_elem_m = zeros(num_elem, num_nodes, num_comp, num_steps);
          offset_elem = getfield(dof_map.submesh.offset.elements, elem_type{j});

          for i=1:numel(sol_tot.bodies)
            if (isfield(sol_tot.bodies(i).stress, stress_type{k}))
              stress_type_b = getfield(sol_tot.bodies(i).stress, stress_type{k});
              
              if (isfield(stress_type_b, elem_type{j}))
                stress_elem_b = getfield(stress_type_b, elem_type{j});
                stress_elem_m(offset_elem(i) + (1:size(stress_elem_b, 1)), :, :, :) = stress_elem_b;
              endif
            endif
          endfor

          stress_elem_type_m = setfield(stress_elem_type_m, elem_type{j}, stress_elem_m);
        endif
      endfor

      if (numel(fieldnames(stress_elem_type_m)))
        sol.stress = setfield(sol.stress, stress_type{k}, stress_elem_type_m);
      endif
    endfor
  endif
endfunction

%!demo
%! ## DEMO
%! close all;
%! number_of_modes = 10;
%! a = [30e-3, 100e-3];
%! b = [20e-3, 20e-3];
%! c = [10e-3, 5e-3];
%! d = [0, 40e-3];
%! h = 3.5e-3;
%! p = [25e6, 25e6];
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! for i=1:numel(a)
%!   unwind_protect
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "wt");
%!       if fd == -1
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       unwind_protect
%!         fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!         fprintf(fd, "a=%g;\n", a(i));
%!         fprintf(fd, "b=%g;\n", b(i));
%!         fprintf(fd, "c=%g;\n", c(i));
%!         fprintf(fd, "d=%g;\n", d(i));
%!         fprintf(fd, "h = %g;\n", h);
%!         fputs(fd, "Point(1) = {d,0.0,0.0,h};\n");
%!         fputs(fd, "Point(2) = {a+d,0.0,0.0,h};\n");
%!         fputs(fd, "Point(3) = {a+d,b,0.0,h};\n");
%!         fputs(fd, "Point(4) = {d,b,0.0,h};\n");
%!         fputs(fd, "Line(1) = {4,3};\n");
%!         fputs(fd, "Line(2) = {3,2};\n");
%!         fputs(fd, "Line(3) = {2,1};\n");
%!         fputs(fd, "Line(4) = {1,4};\n");
%!         fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!         fputs(fd, "Plane Surface(6) = {5};\n");
%!         fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!         fputs(fd, "  Surface{6};\n");
%!         fputs(fd, "};\n");
%!         fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!         fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[4]};\n");
%!         fputs(fd, "Physical Surface(\"load\",2) = {tmp[3]};\n");
%!       unwind_protect_cleanup
%!         fclose(fd);
%!       end_unwind_protect
%!       fprintf(stderr, "meshing ...\n");
%!       pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", [filename, ".geo"]});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         warning("gmsh failed with status %d", status);
%!       endif
%!     unwind_protect_cleanup
%!       unlink([filename, ".geo"]);
%!     end_unwind_protect
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh_data(i).mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%!     fprintf(stderr, "%d nodes\n", rows(mesh_data(1).mesh.nodes));
%!   unwind_protect_cleanup
%!     unlink([filename, ".msh"]);
%!   end_unwind_protect
%! endfor
%! for i=1:numel(mesh_data)
%!   fprintf(stderr, "assembling matrices ...\n");
%!   mesh_data(i).load_case.locked_dof = false(rows(mesh_data(i).mesh.nodes), 6);
%!   mesh_data(i).load_case.locked_dof(mesh_data(i).mesh.groups.tria6(find([[mesh_data(i).mesh.groups.tria6].id] == 1)).nodes, :) = true;
%!   mesh_data(i).load_case.pressure.tria6.elements = mesh_data(i).mesh.elements.tria6(mesh_data(i).mesh.groups.tria6(find([mesh_data(i).mesh.groups.tria6.id] == 2)).elements, :);
%!   mesh_data(i).load_case.pressure.tria6.p = repmat(p(i),size(mesh_data(i).load_case.pressure.tria6.elements)); 
%!   mesh_data(i).mesh.materials.tet10 = ones(rows(mesh_data(i).mesh.elements.tet10), 1, "int32");
%!   mesh_data(i).mesh.material_data.E = 210000e6;
%!   mesh_data(i).mesh.material_data.nu = 0.3;
%!   mesh_data(i).mesh.material_data.rho = 7850;
%!   mesh_data(i).mesh.material_data.C = fem_pre_mat_isotropic(mesh_data(i).mesh.material_data.E, mesh_data(i).mesh.material_data.nu);
%!   mesh_data(i).dof_map = fem_ass_dof_map(mesh_data(i).mesh, mesh_data(i).load_case);
%!   [mesh_data(i).mat_ass.K, ...
%!    mesh_data(i).mat_ass.R] = fem_ass_matrix(mesh_data(i).mesh, ...
%!                                             mesh_data(i).dof_map, ...
%!                                             [FEM_MAT_STIFFNESS, ...
%!                                              FEM_VEC_LOAD_CONSISTENT], mesh_data(i).load_case);
%!   sol.bodies(i).def = fem_sol_static(mesh_data(i).mesh, mesh_data(i).dof_map, mesh_data(i).mat_ass).def;
%!   sol.bodies(i).stress = fem_ass_matrix(mesh_data(i).mesh, ...
%!                                         mesh_data(i).dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         mesh_data(i).load_case, ...
%!                                         sol.bodies(i));
%! endfor
%! [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data);
%! [sol_comb] = fem_post_cms_sol_merge(mesh_comb, dof_map_comb, sol);
%! opts.rotation_angle = [0, 0, 0];
%! opts.print_and_exit = 1;
%! opts.print_to_file = filename;
%! opts.skin_only = true;
%! opts.show_element = true;
%! unwind_protect
%!   fem_post_sol_external(mesh_comb, sol_comb, opts);
%!   [img, map, alpha] = imread([opts.print_to_file, "_001.jpg"]);
%!   figure("visible", "off");
%!   imshow(img, map);
%!   title("Gmsh - stress and deformation");
%! unwind_protect_cleanup
%!   unlink([opts.print_to_file, "_001.jpg"]);
%! end_unwind_protect
%! figure_list();
%! unwind_protect_cleanup
%! if (numel(filename))
%!   fn = dir([filename, "*"]);
%!   for i=1:numel(fn)
%!     unlink(fullfile(fn(i).folder, fn(i).name));
%!   endfor
%! endif
%! end_unwind_protect
