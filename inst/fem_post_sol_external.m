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
## @deftypefn{Function File} fem_post_sol_external(@var{mesh}, @var{sol}, @var{varargin})
## Do postprocessing of a finite element solution within Gmsh.
##
## @var{mesh} @dots{} Finite element mesh data structure.
##
## @var{sol} @dots{} Finite element solution data structure.
##
## @var{varargin} @dots{} Additional arguments passed to fem_post_sol_export
##
## @seealso{fem_post_sol_export}
## @end deftypefn

function fem_post_sol_external(mesh, sol, varargin)
  if (nargin < 2 || nargin > 3 || nargout > 0)
    print_usage();
  endif
  
  prefix = "";

  unwind_protect
    prefix = tempname();

    if (ispc())
      prefix(prefix == "\\") = "/";
    endif
    
    post_pro_geo = fem_post_sol_export(prefix, mesh, sol, varargin{:});

    pid = spawn("gmsh", {post_pro_geo});

    status = spawn_wait(pid);
    
    fprintf(stderr, "gmsh returned with status %d\n", status);

  unwind_protect_cleanup
    if (numel(prefix))
      files = dir([prefix, "*"]);

      for i=1:numel(files)
        [err, msg] = unlink(fullfile(files(i).folder, files(i).name));

        if (err ~= 0)
          warning("failed to remove file \"%s\": %d (%s)", files(i).name, err, msg);
        endif
      endfor
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
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(3).nodes, 1:3) = true;
%!   load_case.loaded_nodes = int32(rows(mesh.nodes));
%!   load_case.loads = [0, 0, F1, 0, 0, 0];
%!   elem.sfncon6.slave = mesh.groups.tria6(1).nodes(:);
%!   elem.sfncon6.master = mesh.elements.tria6(mesh.groups.tria6(2).elements, :);
%!   elem.sfncon6.maxdist = param.g * (1 + sqrt(eps));
%!   mesh.elements.joints = fem_pre_mesh_constr_surf_to_node(mesh.nodes, elem);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
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
%!   opt_post.print_and_exit = true;
%!   opt_post.print_to_file = [filename, "_post_pro_jpg"];
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   opt_post.scale_def = 300;
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   fn = dir([filename, "_post_pro_jpg*.jpg"]);
%!   for i=1:numel(fn)
%!    [img, alpha, map] = imread(fullfile(fn(i).folder, fn(i).name));
%!    figure("visible", "off");
%!    imshow(img, map);
%!    title("van Mises stress - Gmsh");
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
