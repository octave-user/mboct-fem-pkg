## Copyright (C) 2018(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{mesh} = fem_pre_mesh_unstruct_create(@var{geo_file}, @var{param_dim}, @var{options})
## Create an unstructured tetrahedral or hexahedral mesh using Gmsh.
##
## @var{geo_file} @dots{} Gmsh script file which builds a 3D CAD model based on parametric dimensions passed in @var{param_dim}.
##
## @var{param_dim} @dots{} Scalar struct with variables to be passed to the Gmsh script.
##
## @var{options}.verbose @dots{} Enable verbose output.
##
## @var{options}.output_file @dots{} Prefix of a temporary file where the mesh is written to.
##
## @var{options}.mesh.element_size @dots{} Define the global mesh element size.
##
## @var{options}.mesh.jacobian_range @dots{} Define relative limits for the Jacobian determinant.
## @seealso{fem_pre_mesh_struct_create}
## @end deftypefn

function mesh = fem_pre_mesh_unstruct_create(geo_file, param_dim, options)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  if (~isfield(options, "interactive"))
    options.interactive = false;
  endif
  
  f_temp_files = false;
  geo_file_tmp = "";
  mesh_file = "";

  unwind_protect
    if (~isfield(options, "output_file"))
      options.output_file = tempname();

      if (ispc())
        options.output_file(options.output_file == "\\") = "/";
      endif

      f_temp_files = true;
    endif

    geo_file_tmp = [options.output_file, ".geo"];
    mesh_file = [options.output_file, ".msh"];

    [info, err, msg] = stat(geo_file);

    if (err ~= 0)
      error("geometry definition file \"%s\" not found", geo_file);
    endif

    fem_pre_geo_export(geo_file_tmp, geo_file, param_dim);

    if (options.verbose)
      fprintf(stderr, ...
              "generating mesh: input file \"%s\", mesh file \"%s\"\n", ...
              geo_file_tmp, ...
              mesh_file);
    endif

    [info, err, msg] = stat(mesh_file);

    if (err == 0)
      [err, msg] = unlink(mesh_file);
      if (err ~= 0)
        error("failed to access file \"%s\"", mesh_file);
      endif
    endif

    if (options.verbose)
      tic();
    endif

    if (~isfield(options, "mesh"))
      options.mesh = struct();
    endif

    if (~isfield(options.mesh, "order"))
      options.mesh.order = 2;
    endif

    if (~isfield(options.mesh, "dim"))
      options.mesh.dim = 3;
    endif

    cmdline = {"-format", "msh2", ...
               sprintf("-%d", options.mesh.dim), ...
               "-order", sprintf("%d", options.mesh.order)};

    if (~isfield(options.mesh, "optimize"))
      options.mesh.optimize = false;
    endif

    if (options.mesh.optimize)
      if (options.mesh.order > 1)
        cmdline{end + 1} = "-optimize_ho";
      else
        cmdline{end + 1} = "-optimize";
      endif
    endif

    if (options.interactive)
      pid = spawn("gmsh", {geo_file_tmp});
      
      status = spawn_wait(pid);
    endif
    
    if (isfield(options, "mesh"))
      if (isfield(options.mesh, "element_size"))
        cmdline{end + 1} = "-clmin";
        cmdline{end + 1} = sprintf("%g", 0.75 * options.mesh.element_size);
        cmdline{end + 1} = "-clmax";
        cmdline{end + 1} = sprintf("%g", 1.25 * options.mesh.element_size);
      endif

      if (isfield(options.mesh, "jacobian_range"))
        cmdline{end + 1} = "-ho_min";
        cmdline{end + 1} = sprintf("%g", options.mesh.jacobian_range(1));
        cmdline{end + 1} = "-ho_max";
        cmdline{end + 1} = sprintf("%g", options.mesh.jacobian_range(2));
      endif
    endif

    cmdline{end + 1} = geo_file_tmp;
    cmdline{end + 1} = "-o";
    cmdline{end + 1} = mesh_file;

    pid = spawn("gmsh", cmdline);

    status = spawn_wait(pid);

    if (options.verbose)
      toc();
    endif

    if (status ~= 0)
      warning("gmsh returned with exit status %d", status);
    endif

    [info, err, msg] = stat(mesh_file);

    if (err ~= 0)
      error("gmsh failed to generate mesh for input file \"%s\"", geo_file_tmp);
    endif

    if (options.verbose)
      fprintf(stderr, "loading mesh file \"%s\"\n", mesh_file);
      tic();
    endif

    mesh = fem_pre_mesh_import(mesh_file, "gmsh", options.mesh);

    if (options.verbose)
      toc();
      fprintf(stderr, "%d nodes generated\n", rows(mesh.nodes));
    endif
  unwind_protect_cleanup
    if (f_temp_files)
      if (numel(geo_file_tmp))
        [~] = unlink(geo_file_tmp);
      endif
      if (numel(mesh_file))
        [~] = unlink(mesh_file);
      endif
    endif
  end_unwind_protect
endfunction

%!demo
%! ## DEMO1
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.15e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0.01e-3; ## FIXME: di=0 not supported by Gmsh 4.5.5
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.22;
%! B = 1.37;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.D^2 * pi / 4);
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
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(8) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Circle(4) = {4,5,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,8};\n");
%!     fputs(fd, "Line(7) = {8,1};\n");
%!     fputs(fd, "Line Loop(9) = {1,2,3,4,5,6,7};\n");
%!     fputs(fd, "Plane Surface(11) = {9};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{11}; Layers{Ceil(d * Pi / 4 / h0)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{11,v[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {11};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[7]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 1;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.iso4(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.iso4(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.iso4(3).nodes, 1) = true;
%!   load_case.pressure.iso4.elements = mesh.elements.iso4(mesh.groups.iso4(4).elements, :);
%!   load_case.pressure.iso4.p = repmat(p, rows(load_case.pressure.iso4.elements), 4);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   Kt = max(max(sol_stat.stress.vmis.iso8)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.15 * Kt_a);
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
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.15e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0.01e-3; ## FIXME: di=0 not supported by Gmsh 4.5.5
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.62;
%! B = 3.5;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.d^2 * pi / 4);
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
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * L, 0.0, 0.5 * d, h1};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Line(4) = {4,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,1};\n");
%!     fputs(fd, "Line Loop(7) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(8) = {7};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{8}; Layers{Ceil(d * Pi / 4 / h0)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{8,v[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {8};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[6]};\n");
%!   unwind_protect_cleanup
%!     fclose(fd);
%!   end_unwind_protect
%!   opt.mesh.order = 1;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.iso4(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.iso4(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.iso4(3).nodes, 1) = true;
%!   load_case.pressure.iso4.elements = mesh.elements.iso4(mesh.groups.iso4(4).elements, :);
%!   load_case.pressure.iso4.p = repmat(p, rows(load_case.pressure.iso4.elements), 4);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.iso8)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.15 * Kt_a);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO3
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.15e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0e-3;
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.22;
%! B = 1.37;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.D^2 * pi / 4);
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
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(8) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Circle(4) = {4,5,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,8};\n");
%!     fputs(fd, "Line(7) = {8,1};\n");
%!     fputs(fd, "Line Loop(9) = {1,2,3,4,5,6,7};\n");
%!     fputs(fd, "Plane Surface(11) = {9};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{11};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {11};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[7]};\n");
%!   unwind_protect_cleanup
%!     fclose(fd);
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.tria6(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.tria6(3).nodes, 1) = true;
%!   load_case.pressure.tria6.elements = mesh.elements.tria6(mesh.groups.tria6(4).elements, :);
%!   load_case.pressure.tria6.p = repmat(p, rows(load_case.pressure.tria6.elements), 6);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   Kt = max(max(sol_stat.stress.vmis.tet10)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.08 * Kt_a);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO4
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.15e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0e-3;
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.62;
%! B = 3.5;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.d^2 * pi / 4);
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
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * L, 0.0, 0.5 * d, h1};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Line(4) = {4,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,1};\n");
%!     fputs(fd, "Line Loop(7) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(8) = {7};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{8};\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {8};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[6]};\n");
%!   unwind_protect_cleanup
%!     fclose(fd);
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.tria6(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.tria6(3).nodes, 1) = true;
%!   load_case.pressure.tria6.elements = mesh.elements.tria6(mesh.groups.tria6(4).elements, :);
%!   load_case.pressure.tria6.p = repmat(p, rows(load_case.pressure.tria6.elements), 6);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.tet10)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.05 * Kt_a);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO5
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.2e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0.01e-3; ## FIXME: di=0 not supported by Gmsh 4.6.0 (degenerated hexahedrons are created)
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.22;
%! B = 1.37;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.D^2 * pi / 4);
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
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(8) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Circle(4) = {4,5,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,8};\n");
%!     fputs(fd, "Line(7) = {8,1};\n");
%!     fputs(fd, "Line Loop(9) = {1,2,3,4,5,6,7};\n");
%!     fputs(fd, "Plane Surface(11) = {9};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{11}; Layers{Ceil(d * Pi / 4 / h0)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{11,v[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {11};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[7]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.quad8(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.quad8(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.quad8(3).nodes, 1) = true;
%!   load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(4).elements, :);
%!   load_case.pressure.quad8.p = repmat(p, rows(load_case.pressure.quad8.elements), 8);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   Kt = max(max(sol_stat.stress.vmis.iso20)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.15 * Kt_a);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO6
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.2e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0.01e-3; ## FIXME: di=0 not supported by Gmsh 4.5.5
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.62;
%! B = 3.5;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.d^2 * pi / 4);
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
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * L, 0.0, 0.5 * d, h1};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Line(4) = {4,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,1};\n");
%!     fputs(fd, "Line Loop(7) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(8) = {7};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{8}; Layers{Ceil(d * Pi / 4 / h0)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Recombine Surface{8,v[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume\") = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\") = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\") = {8};\n");
%!     fputs(fd, "Physical Surface(\"clamp\") = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\") = {v[6]};\n");
%!   unwind_protect_cleanup
%!     fclose(fd);
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.quad8(1).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.quad8(2).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.quad8(3).nodes, 1) = true;
%!   load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(4).elements, :);
%!   load_case.pressure.quad8.p = repmat(p, rows(load_case.pressure.quad8.elements), 8);
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   Kt = max(max(sol_stat.stress.vmis.iso20)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.15 * Kt_a);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO7
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h0 = 0.2e-3;
%! geo.h1 = 2e-3;
%! geo.D = 15e-3;
%! geo.d = 13e-3;
%! geo.di = 0.01e-3; ## FIXME: di=0 not supported
%! geo.r = 0.5 * (geo.D - geo.d);
%! geo.L = 2 * geo.D;
%! geo.t = 0.5 * (geo.D - geo.d);
%! geo.w = 2 * sqrt(geo.r^2 - (geo.r - geo.t)^2);
%! A = 0.22;
%! B = 1.37;
%! Kt_a = 1 + 1 / sqrt(A * geo.r / geo.t + 2 * B * geo.r / geo.d * (1 + 2 * geo.r / geo.d)^2);
%! tauxx_n = F1 / (geo.d^2 * pi / 4);
%! tauxx_a = tauxx_n * Kt_a;
%! p = -F1 / (geo.D^2 * pi / 4);
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
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(4) = {     0.0, 0.0, 0.5 * d, h0};\n");
%!     fputs(fd, "Point(5) = {     0.0, 0.0, 0.5 * d + r, h0};\n");
%!     fputs(fd, "Point(6) = { 0.5 * w, 0.0, 0.5 * D, h0};\n");
%!     fputs(fd, "Point(7) = { 0.5 * L, 0.0, 0.5 * D, h1};\n");
%!     fputs(fd, "Point(8) = { 0.5 * L, 0.0, 0.5 * di, h1};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Circle(3) = {3,5,4};\n");
%!     fputs(fd, "Circle(4) = {4,5,6};\n");
%!     fputs(fd, "Line(5) = {6,7};\n");
%!     fputs(fd, "Line(6) = {7,8};\n");
%!     fputs(fd, "Line(7) = {8,1};\n");
%!     fputs(fd, "Line Loop(9) = {1,2,3,4,5,6,7};\n");
%!     fputs(fd, "Plane Surface(11) = {9};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{11}; Layers{Ceil(d * Pi / 4 / h0)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",0) = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\",1) = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\",2) = {11};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {v[7]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   elem_types = fieldnames(mesh.elements);
%!   mesh.materials = struct();
%!   for i=1:numel(elem_types)
%!     mesh.materials = setfield(mesh.materials, elem_types{i}, ones(rows(getfield(mesh.elements, elem_types{i})), 1, "int32"));
%!   endfor
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==1)).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==2)).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.quad8(find([mesh.groups.quad8.id]==3)).nodes, 1) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==3)).nodes, 1) = true;
%!   idx = find([mesh.groups.quad8.id]==4);
%!   if (numel(idx))
%!     load_case.pressure.quad8.elements = mesh.elements.quad8(mesh.groups.quad8(idx).elements, :);
%!     load_case.pressure.quad8.p = repmat(p, rows(load_case.pressure.quad8.elements), 8);
%!   endif
%!   idx = find([mesh.groups.tria6.id]==4);
%!   if (numel(idx))
%!     load_case.pressure.tria6.elements = mesh.elements.tria6(mesh.groups.tria6(idx).elements, :);
%!     load_case.pressure.tria6.p = repmat(p, rows(load_case.pressure.tria6.elements), 6);
%!   endif
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_SCA_STRESS_VMIS], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   Kt = max(max(sol_stat.stress.vmis.penta15)) / tauxx_n;
%!   fprintf(stdout, "Kt_a=%.2f\n", Kt_a);
%!   fprintf(stdout, "Kt=%.2f\n", Kt);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
%!   figure_list();
%!   assert(Kt, Kt_a, 0.15 * Kt_a);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO8
%! close all;
%! E = 210000e6;
%! nu = 0.3;
%! rho = 7850;
%! F1 = 120;
%! geo.h = 0.25e-3;
%! geo.Do = 4e-3;
%! geo.Di = 2e-3;
%! geo.L = 4e-3;
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
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Point(1) = {-0.5 * L, 0.0, 0.5 * Di, h};\n");
%!     fputs(fd, "Point(2) = {-0.5 * L, 0.0, 0.5 * Do, h};\n");
%!     fputs(fd, "Point(3) = { 0.5 * L, 0.0, 0.5 * Do, h};\n");
%!     fputs(fd, "Point(4) = { 0.5 * L, 0.0, 0.5 * Di, h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "Line(2) = {2,3};\n");
%!     fputs(fd, "Line(3) = {3,4};\n");
%!     fputs(fd, "Line(4) = {4,1};\n");
%!     fputs(fd, "Line Loop(5) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "v[] = Extrude {{1.0,0.0,0.0},{0.0,0.0,0.0},-Pi/2} {\n");
%!     fputs(fd, "  Surface{6}; Layers{Ceil(Do * Pi / 4 / h)}; Recombine;\n");
%!     fputs(fd, "};\n");
%!     fputs(fd, "Physical Volume(\"volume\",0) = {v[1]};\n");
%!     fputs(fd, "Physical Surface(\"bottom\",1) = {v[0]};\n");
%!     fputs(fd, "Physical Surface(\"front\",2) = {6};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {v[2]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {v[4]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   opt.mesh.elem_type = {"tria6", "quad8", "penta15", "tet10h"};
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   elem_types = fieldnames(mesh.elements);
%!   mesh.materials = struct();
%!   for i=1:numel(elem_types)
%!     mesh.materials = setfield(mesh.materials, elem_types{i}, ones(rows(getfield(mesh.elements, elem_types{i})), 1, "int32"));
%!   endfor
%!   load_case.locked_dof = false(size(mesh.nodes));
%!   node_id_disp = [];
%!   if (isfield(mesh.groups, "tria6"))
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==1)).nodes, 3) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==2)).nodes, 2) = true;
%!   load_case.locked_dof(mesh.groups.tria6(find([mesh.groups.tria6.id]==3)).nodes, 1) = true;
%!   node_id_disp = unique([node_id_disp, mesh.groups.tria6(find([mesh.groups.tria6.id]==4)).nodes]);
%!   endif
%!   if (isfield(mesh.groups, "quad8"))
%!   load_case.locked_dof(mesh.groups.quad8(find([mesh.groups.quad8.id]==3)).nodes, 1) = true;
%!   node_id_disp = unique([node_id_disp, mesh.groups.quad8(find([mesh.groups.quad8.id]==4)).nodes]);
%!   endif
%!   empty_cell = cell(1, numel(node_id_disp));
%!   mesh.elements.joints = struct("C", empty_cell, "nodes", empty_cell);
%!   load_case.joints = struct("U", empty_cell);
%!   for i=1:numel(node_id_disp)
%!     mesh.elements.joints(i).C = [1, 0, 0, 0, 0, 0];
%!     mesh.elements.joints(i).nodes = node_id_disp(i);
%!     load_case.joints(i).U = geo.L / E;
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mtot, ...
%!    mat_info, ...
%!    mesh_info] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_STIFFNESS, ...
%!                                 FEM_VEC_LOAD_CONSISTENT, ...
%!                                 FEM_SCA_TOT_MASS], ...
%!                                load_case);
%!   opt_sol.number_of_threads = int32(4);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass, opt_sol);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   if (mesh_info.detJ.min <= 0)
%!     error("Jacobian is singular");
%!   endif
%!   m_a = (geo.Do^2 - geo.Di^2) * pi / 4 * geo.L / 4 * rho;
%!   assert(mat_ass.mtot, m_a, eps^0.4 * m_a);
%!   opt_post.scale_def = 3000;
%!   opt_post.show_element = true;
%!   opt_post.print_to_file = [filename, "_post"];
%!   opt_post.print_and_exit = true;
%!   opt_post.rotation_angle = [-pi/2, 0, 0];
%!   fem_post_sol_external(mesh, sol_stat, opt_post);
%!   figure("visible", "off");
%!   [img, map, alpha] = imread([opt_post.print_to_file, "_001.jpg"]);
%!   imshow(img, alpha);
%!   title("deformed mesh/van Mises stress");
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
%! ## DEMO9
%! close all;
%! SI_unit_meter = 1;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1;
%! geo.D = 148e-3 / SI_unit_meter;
%! geo.d = 140e-3 / SI_unit_meter;
%! geo.h = 10e-3 / SI_unit_meter;
%! param.rho = 1.33 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.c = 225 / (SI_unit_meter / SI_unit_second);
%! param.eta = 8.3492e-06 / (SI_unit_pascal * SI_unit_second);
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
%!     fputs(fd, "vout = newv;\n");
%!     fputs(fd, "Sphere(vout) = {0, 0, 0, 0.5 * D};\n");
%!     fputs(fd, "vin = newv;\n");
%!     fputs(fd, "Sphere(vin) = {0, 0, 0, 0.5 * d};\n");
%!     fputs(fd, "v = newv;\n");
%!     fputs(fd, "BooleanDifference(v) = {Volume{vout}; Delete; }{ Volume{vin}; Delete; };\n");
%!     fputs(fd, "Physical Volume(\"volume\",1) = {v};\n");
%!     fputs(fd, "Physical Surface(\"surface\", 8) = {2};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v}; } } = h;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   opt.mesh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh = fem_pre_mesh_reorder(mesh);
%!   mesh.material_data.c = 220;
%!   mesh.material_data.rho = 1.33;
%!   load_case.locked_dof = false(rows(mesh.nodes), 1);
%!   load_case.domain = FEM_DO_ACOUSTICS;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.c = param.c;
%!   mesh.material_data.eta = param.eta;
%!   dof_map = fem_ass_dof_map(mesh, load_case);
%!   [mat_ass.Ka, ...
%!    mat_ass.Ma, ...
%!    mat_ass.coll_Ka, ...
%!    mat_ass.coll_Ma] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_ACOUSTICS_RE, ...
%!                                       FEM_MAT_MASS_ACOUSTICS_RE, ...
%!                                       FEM_VEC_COLL_STIFF_ACOUSTICS, ...
%!                                       FEM_VEC_COLL_MASS_ACOUSTICS], ...
%!                                      load_case);
%!   N = 80;
%!   fref = 2 * pi * 1 / (1 / SI_unit_second);
%!   opt_sol.rho = fref^2;
%!   opt_sol.tolerance = sqrt(eps);
%!   opt_sol.algorithm = "shift-invert";
%!   opt_sol.solver = "umfpack";
%!   opt_sol.number_of_threads = int32(6);
%!   opt_sol.refine_max_iter = int32(10);
%!   opt_sol.pre_scaling = true;
%!   opt_so.problem = "acoustic";
%!   [Phi, lambda, err] = fem_sol_eigs(mat_ass.Ka, mat_ass.Ma, N, opt_sol);
%!   sol.p = -Phi(dof_map.ndof, :) * diag(imag(lambda));
%!   sol.f = imag(lambda) / (2 * pi);
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO10
%! close all;
%! SI_unit_meter = 1;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! geo.D = 148e-3 / SI_unit_meter;
%! geo.d = 136e-3 / SI_unit_meter;
%! geo.h = 10e-3 / SI_unit_meter;
%! geo.x1 = 0 / SI_unit_meter;
%! geo.y1 = 0 / SI_unit_meter;
%! geo.z1 = 100e-3 / SI_unit_meter;
%! geo.x2 = 0 / SI_unit_meter;
%! geo.y2 = 0 / SI_unit_meter;
%! geo.z2 = geo.z1;
%! geo.lsuc = 300e-3 / SI_unit_meter;
%! geo.dsuc = 6e-3 / SI_unit_meter;
%! param.fmin = 0;
%! param.fmax = 2000 / (SI_unit_second^-1);
%! param.maxdef = 10e-3 / SI_unit_meter;
%! param.num_freq_modal = 50000;
%! param.rho = 1.33 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.c = 225 / (SI_unit_meter / SI_unit_second);
%! param.eta = 8.3492e-06 / (SI_unit_pascal * SI_unit_second);
%! param.m1 = 2.2 / SI_unit_kilogram;
%! param.J1 = 1e-6 * [2.4427674e+03 -3.2393301e+01 -5.3623318e+02
%!                    -3.2393301e+01  3.7506484e+03  7.9745924e+01
%!                    -5.3623318e+02  7.9745924e+01  4.2943071e+03] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg1 = zeros(3, 1);
%! param.m2 = 2.2978223 / SI_unit_kilogram;
%! param.J2 = 1e-6 * [8.6104517e+03, 5.5814351e+01, -3.0103453e+02;
%!                    5.5814351e+01, 1.1905289e+04,  2.0425595e+02;
%!                    -3.0103453e+02, 2.0425595e+02,  1.2109596e+04] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg2 = 1e-3 * [12.940131;
%!                      -6.0192164e-01;
%!                      5.9683120e+01] / SI_unit_meter;
%! param.offset1 = (6.63e-3 + 2.4e-3) / SI_unit_meter;
%! param.N = int32(60);
%! param.use_impedance = false;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.D = 12.12e-3 / SI_unit_meter + helspr1.d;
%! helspr1.L = 27.7e-3 / SI_unit_meter;
%! helspr1.n = 5.7;
%! helspr1.ni = 2.7;
%! helspr1.ng = 0.75;
%! helspr1.m = 60;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! helspr1.material.nu = helspr1.material.E / (2 * helspr1.material.G) - 1;
%! helspr1.material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! [helspr1.material.alpha, helspr1.material.beta] = fem_pre_mat_rayleigh_damping([1e-2; 0.03e-2], [20, 1000] / (SI_unit_second^-1));
%! helspr1.X = [[45.20e-3,  45.20e-3, -62.5e-3;
%!               43.13e-3, -43.13e-3,   0.0e-3] / SI_unit_meter;
%!              [-helspr1.L, -helspr1.L, -helspr1.L]];
%! section1.A = helspr1.d^2 * pi / 4;
%! section1.Ay = 9 / 10 * section1.A;
%! section1.Az = section1.Ay;
%! section1.Iy = helspr1.d^4 * pi / 64;
%! section1.Iz = section1.Iy;
%! section1.It = section1.Iy + section1.Iz;
%! rubspr1.d = 11.1e-3 / SI_unit_meter;
%! rubspr1.D = 28.6e-3 / SI_unit_meter;
%! rubspr1.L = 10.5e-3 / SI_unit_meter;
%! rubspr1.material.E = 2.6e6 / SI_unit_pascal;
%! rubspr1.material.nu = 0.499;
%! rubspr1.material.rho = 910 / (SI_unit_kilogram / SI_unit_meter^3);
%! [rubspr1.material.alpha, rubspr1.material.beta] = fem_pre_mat_rayleigh_damping([10e-2; 30e-2], [30, 100] / (SI_unit_second^-1));
%! rubspr1.X = [(170e-3 * [0.5, 0.5, -0.5, -0.5] + 8e-3) / SI_unit_meter;
%!              70e-3 * [-0.5, 0.5, -0.5, 0.5] / SI_unit_meter;
%!              -[1, 1, 1, 1] * (helspr1.L + rubspr1.L + param.offset1)];
%! rubspr1.e2 = [1, 0, 0];
%! section2.A = (rubspr1.D^2 - rubspr1.d^2) * pi / 4;
%! section2.Ay = 0.8 * section2.A;
%! section2.Az = section2.Ay;
%! section2.Iy = (rubspr1.D^4 - rubspr1.d^4) * pi / 64;
%! section2.Iz = section2.Iy;
%! section2.It = section2.Iy + section2.Iz;
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
%!     fputs(fd, "vout = newv;\n");
%!     fputs(fd, "Sphere(vout) = {x2, y2, z2, 0.5 * D};\n");
%!     fputs(fd, "vsuc = newv;\n");
%!     fputs(fd, "Cylinder(vsuc) = {x2, y2, z2, lsuc, 0, 0, 0.5 * dsuc};\n");
%!     fputs(fd, "vouts = newv;\n");
%!     fputs(fd, "BooleanUnion(vouts) = { Volume{vout}; Delete; }{ Volume{vsuc}; Delete; };\n")
%!     fputs(fd, "vin = newv;\n");
%!     fputs(fd, "Sphere(vin) = {x1, y1, z1, 0.5 * d};\n");
%!     fputs(fd, "v = newv;\n");
%!     fputs(fd, "BooleanDifference(v) = {Volume{vouts}; Delete; }{ Volume{vin}; Delete; };\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {v};\n");
%!     fputs(fd, "Physical Surface(\"inside\", 2) = {4};\n");
%!     fputs(fd, "Physical Surface(\"outside\", 3) = {1};\n");
%!     fputs(fd, "Physical Surface(\"inlet\", 4) = {3};\n");
%!     fputs(fd, "Physical Surface(\"tube\", 5) = {2};\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v}; } } = h;\n");
%!     fputs(fd, "MeshSize{PointsOf{Surface{2}; } } = dsuc * Pi / 7.;\n");
%!     fputs(fd, "ReorientMesh Volume{v};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   opt.mesh.order = 2;
%!   opt.mesh.elem_type = {"tria6h", "tet10h"};
%!   opt.interactive = false;
%!   mesh = fem_pre_mesh_unstruct_create(geo_file, geo, opt);
%!   mesh = fem_pre_mesh_reorder(mesh);
%!   node_idx_rb1 = int32(rows(mesh.nodes) + 1);
%!   node_idx_rb2 = int32(rows(mesh.nodes) + 2);
%!   grp_idx_fsi1 = int32(find([mesh.groups.tria6h.id] == 2));
%!   grp_idx_fsi2 = int32(find([mesh.groups.tria6h.id] == 3));
%!   grp_idx_inlet = int32(find([mesh.groups.tria6h.id] == 4));
%!   mesh.nodes(node_idx_rb1, :) = [0, 0, rubspr1.L + param.offset1 + helspr1.L, zeros(1, 3)];
%!   mesh.nodes(node_idx_rb2, :) = [0, 0, rubspr1.L, zeros(1, 3)];
%!   node_idx_helspr1 = int32(rows(mesh.nodes) + 1);
%!   helspr1.Phi = linspace(0, 2 * pi * helspr1.n, ceil(helspr1.m * helspr1.n)).' + 2 * pi * helspr1.ni;
%!   node_idx_rubspr1 = node_idx_helspr1 + numel(helspr1.Phi) * columns(helspr1.X);
%!   helspr1.x = 0.5 * helspr1.D * cos(helspr1.Phi);
%!   helspr1.y = 0.5 * helspr1.D * sin(helspr1.Phi);
%!   helspr1.z = (helspr1.L - helspr1.d * (2 * (helspr1.ni - helspr1.ng) + 1)) * linspace(0, 1, numel(helspr1.Phi))(:) + helspr1.d * (helspr1.ni - helspr1.ng + 0.5);
%!   helspr1.e2 = repmat([0, 0, 1], numel(helspr1.Phi) - 1, 1);
%!   idx_beam = int32(0);
%!   for i=1:columns(helspr1.X)
%!     elnodes = int32([(1:numel(helspr1.Phi) - 1).', (2:numel(helspr1.Phi)).']) + node_idx_helspr1 - 1 + (i - 1) * numel(helspr1.Phi);
%!     mesh.nodes(node_idx_helspr1 - 1 + (1:numel(helspr1.Phi)) + (i - 1) * numel(helspr1.Phi), 1:6) = [[helspr1.x, helspr1.y, helspr1.z] + helspr1.X(:, i).', zeros(numel(helspr1.Phi), 3)];
%!     mesh.elements.beam2(idx_beam + (1:numel(helspr1.Phi) - 1)) = struct("nodes", mat2cell(elnodes, ones(numel(helspr1.Phi) - 1, 1, "int32"), 2), ...
%!                                                                         "section", mat2cell(repmat(section1, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32")), ...
%!                                                                         "e2", mat2cell(helspr1.e2, ones(numel(helspr1.Phi) - 1, 1, "int32"), 3));
%!     idx_beam += numel(helspr1.Phi) - 1;
%!   endfor
%!   for i=1:columns(rubspr1.X)
%!     elnodes = int32([1, 2]) + node_idx_rubspr1 - 1 + (i - 1) * 2;
%!     mesh.nodes(node_idx_rubspr1 - 1 + (1:2) + 2 * (i - 1), 1:6) = [[zeros(1, 3); zeros(1, 2), rubspr1.L] + rubspr1.X(:, i).', zeros(2, 3)];
%!     mesh.elements.beam2(idx_beam + 1) = struct("nodes", mat2cell(elnodes, ones(1, 1, "int32"), 2), ...
%!                                                "section", mat2cell(repmat(section2, 1, 1), ones(1, 1, "int32")), ...
%!                                                "e2", mat2cell(rubspr1.e2, ones(1, 1, "int32"), 3));
%!     ++idx_beam;
%!   endfor
%!   empty_cell = cell(1, 2);
%!   mesh.elements.bodies = struct("m", empty_cell, "J", empty_cell, "lcg", empty_cell, "nodes", empty_cell);
%!   mesh.elements.bodies(1).m = param.m1;
%!   mesh.elements.bodies(1).J = param.J1;
%!   mesh.elements.bodies(1).lcg = param.lcg1;
%!   mesh.elements.bodies(1).nodes = node_idx_rb1;
%!   mesh.elements.bodies(2).m = param.m2;
%!   mesh.elements.bodies(2).J = param.J2;
%!   mesh.elements.bodies(2).lcg = param.lcg2;
%!   mesh.elements.bodies(2).nodes = node_idx_rb2;
%!   empty_cell = cell(1, 3);
%!   mesh.material_data = struct("E", empty_cell, "nu", empty_cell, "rho", empty_cell, "c", empty_cell, "eta", empty_cell);
%!   mesh.material_data(1).c = param.c;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).eta = param.eta;
%!   mesh.material_data(2).E = helspr1.material.E;
%!   mesh.material_data(2).nu = helspr1.material.nu;
%!   mesh.material_data(2).beta = helspr1.material.beta;
%!   mesh.material_data(2).alpha = helspr1.material.alpha;
%!   mesh.material_data(2).rho = helspr1.material.rho;
%!   mesh.material_data(3).E = rubspr1.material.E;
%!   mesh.material_data(3).nu = rubspr1.material.nu;
%!   mesh.material_data(3).rho = rubspr1.material.rho;
%!   mesh.material_data(3).beta = rubspr1.material.beta;
%!   mesh.material_data(3).alpha = rubspr1.material.alpha;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 7);
%!   if (~param.use_impedance)
%!     load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_inlet).nodes, 7) = true;
%!   endif
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_fsi1).nodes, 4:6) = true;
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_fsi2).nodes, 4:6) = true;
%!   for i=1:columns(rubspr1.X)
%!     load_case_dof.locked_dof(node_idx_rubspr1 + 2 * (i - 1), 1:6) = true;
%!   endfor
%!   load_case_dof.domain = FEM_DO_FLUID_STRUCT;
%!   mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!   mesh.materials.beam2 = [repmat(int32(2), columns(helspr1.X) * (numel(helspr1.Phi) - 1), 1);
%!                           repmat(int32(3), columns(rubspr1.X), 1)];
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).c = param.c;
%!   mesh.elements.fluid_struct_interface.tria6h = mesh.elements.tria6h([[mesh.groups.tria6h([grp_idx_fsi1, grp_idx_fsi2])].elements], :);
%!   if (param.use_impedance)
%!     mesh.elements.acoustic_impedance.tria6h.nodes = mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_inlet).elements, :);
%!     mesh.elements.acoustic_impedance.tria6h.z = repmat(param.rho * param.c, size(mesh.elements.acoustic_impedance.tria6h.nodes));
%!     mesh.materials.acoustic_impedance.tria6h = ones(rows(mesh.elements.acoustic_impedance.tria6h.nodes), 1, "int32");
%!   endif
%!   empty_cell = cell(1, numel(mesh.groups.tria6h(grp_idx_fsi1).nodes) + numel(mesh.groups.tria6h(grp_idx_fsi2).nodes) + 2 * columns(helspr1.X) + columns(rubspr1.X));
%!   mesh.elements.joints = struct("nodes", empty_cell, "C", empty_cell);
%!   idx_joint = int32(0);
%!   for i=1:numel(mesh.groups.tria6h(grp_idx_fsi1).nodes)
%!     node_idx_slave = mesh.groups.tria6h(grp_idx_fsi1).nodes(i);
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb1, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb1, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [eye(3), -skew(lslave), -eye(3), zeros(3, 3)];
%!   endfor
%!   for i=1:numel(mesh.groups.tria6h(grp_idx_fsi2).nodes)
%!     node_idx_slave = mesh.groups.tria6h(grp_idx_fsi2).nodes(i);
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [eye(3), -skew(lslave), -eye(3), zeros(3, 3)];
%!   endfor
%!   for i=1:columns(helspr1.X)
%!     node_idx_slave = node_idx_helspr1 - 1 + numel(helspr1.Phi) * i;
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb1, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb1, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                               zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%!   endfor
%!   for i=1:columns(helspr1.X)
%!     node_idx_slave = node_idx_helspr1 + numel(helspr1.Phi) * (i - 1);
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                               zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%!   endfor
%!   for i=1:columns(rubspr1.X)
%!     node_idx_slave = node_idx_rubspr1 + 2 * i - 1;
%!     lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!     ++idx_joint;
%!     mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!     mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                               zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%!   endfor
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   empty_cell = cell(1, 6);
%!   load_case = struct("loads", empty_cell, "loaded_nodes", empty_cell);
%!   idx_load = int32(0);
%!   for i=1:6
%!     ++idx_load;
%!     load_case(idx_load).loaded_nodes = node_idx_rb1;
%!     load_case(idx_load).loads = zeros(1, 6);
%!     load_case(idx_load).loads(i) = 1;
%!   endfor
%!   [mat_ass.Mfs_re, ...
%!    mat_ass.Mfs_im, ...
%!    mat_ass.Kfs_re, ...
%!    mat_ass.Kfs_im, ...
%!    mat_ass.Dfs_re, ...
%!    mat_ass.Dfs_im, ...
%!    mat_ass.Rfs, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_MASS_FLUID_STRUCT_IM, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_STIFFNESS_FLUID_STRUCT_IM, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_RE, ...
%!                                         FEM_MAT_DAMPING_FLUID_STRUCT_IM, ...
%!                                         FEM_VEC_LOAD_FLUID_STRUCT], ...
%!                                        load_case);
%!   if (~nnz(mat_ass.Mfs_im))
%!     mat_ass = rmfield(mat_ass, "Mfs_im");
%!   endif
%!   if (~nnz(mat_ass.Dfs_im))
%!     mat_ass = rmfield(mat_ass, "Dfs_im");
%!   endif
%!   if (~nnz(mat_ass.Kfs_im))
%!     mat_ass = rmfield(mat_ass, "Kfs_im");
%!   endif
%!   opt_linsol.solver = "umfpack";
%!   opt_linsol.pre_scaling = true;
%!   opt_linsol.number_of_threads = int32(4);
%!   opt_linsol.refine_max_iter = int32(20);
%!   opt_linsol.epsilon_refinement = 1e-10;
%!   opt_linsol.verbose = int32(0);
%!   opt_linsol.weighted_matching = true;
%!   opt_linsol.scaling = true;
%!   opt_eigs.maxit = int32(100);
%!   opt_eigs.disp = int32(0);
%!   opt_eigs.p = 5 * param.N;
%!   [sol_eig, Phi] = fem_sol_modal_fsi(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eigs);
%!   [Phi, h] = fem_sol_modes_scale(mat_ass.Mfs_re, mat_ass.Kfs_re, sol_eig.lambda, Phi, mat_ass.Rfs);
%!   omega = linspace(2 * pi * param.fmin, 2 * pi * param.fmax, param.num_freq_modal);
%!   U2 = fem_sol_harmonic_modal(h, sol_eig.lambda, Phi(dof_map.ndof(node_idx_rb2, 1:6), :), omega);
%!   for i=1:3
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     hold on;
%!     plot(omega / (2*pi) * (SI_unit_second^-1), 20 * log10(abs(-omega.^2 .* U2(i, :, i)) * (SI_unit_meter / SI_unit_second^2)), "-;modal;r");
%!     xlabel("f [Hz]");
%!     ylabel("a [dB/(1 m/s^2/N)]");
%!     ylim([-80, 10]);
%!     yticks(-100:5:10);
%!     grid on;
%!     grid minor on;
%!     title("frequency response magnitude");
%!     subplot(2, 1, 2);
%!     hold on;
%!     plot(omega / (2*pi) * (SI_unit_second^-1), 180 / pi * arg(-omega.^2 .* U2(i, :, i)), "-;modal;r");
%!     xlabel("f [Hz]");
%!     ylabel("arg(a) [deg]");
%!     ylim([-180,180]);
%!     yticks(-180:30:180);
%!     grid on;
%!     grid minor on;
%!     title("frequency response phase");
%!   endfor
%!   idx_freq = find(sol_eig.f >= 0);
%!   sol_eig.f = sol_eig.f(idx_freq);
%!   sol_eig.D = sol_eig.D(idx_freq);
%!   sol_eig.lambda = sol_eig.lambda(idx_freq);
%!   sol_eig.def = sol_eig.def(:, :, idx_freq);
%!   sol_eig.p = sol_eig.p(:, idx_freq);
%!   for i=1:size(sol_eig.def, 3)
%!     sre = max(max(abs(real(sol_eig.def(:, 1:3, i)))));
%!     sim = max(max(abs(imag(sol_eig.def(:, 1:3, i)))));
%!     if (sim > sre)
%!       sol_eig.def(:, :, i) = imag(sol_eig.def(:, :, i)) * param.maxdef / sim;
%!       sol_eig.p(:, i) = imag(sol_eig.p(:, i)) * param.maxdef / sim;
%!     else
%!       sol_eig.def(:, :, i) = real(sol_eig.def(:, :, i)) * param.maxdef / sre;
%!       sol_eig.p(:, i) = real(sol_eig.p(:, i)) * param.maxdef / sre;
%!     endif
%!   endfor
%!   opt_post.elem_types={"tria6h","beam2"};
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ## DEMO11
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.m1 = 2.2 / SI_unit_kilogram;
%! param.J1 = 1e-6 * [2.4427674e+03 -3.2393301e+01 -5.3623318e+02
%!                    -3.2393301e+01  3.7506484e+03  7.9745924e+01
%!                    -5.3623318e+02  7.9745924e+01  4.2943071e+03] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg1 = zeros(3, 1);
%! param.m2 = 2.2978223 / SI_unit_kilogram;
%! param.J2 = 1e-6 * [8.6104517e+03, 5.5814351e+01, -3.0103453e+02;
%!                    5.5814351e+01, 1.1905289e+04,  2.0425595e+02;
%!                    -3.0103453e+02, 2.0425595e+02,  1.2109596e+04] / (SI_unit_kilogram * SI_unit_meter^2);
%! param.lcg2 = 1e-3 * [12.940131;
%!                      -6.0192164e-01;
%!                      5.9683120e+01] / SI_unit_meter;
%! param.offset1 = (6.63e-3 + 2.4e-3) / SI_unit_meter;
%! param.N = 200;
%! param.maxdef = 10e-3 / SI_unit_meter;
%! param.fmin = 0;
%! param.fmax = 2000 / (SI_unit_second^-1);
%! param.num_freq_modal = 10000;
%! param.num_freq_nodal = 400;
%! helspr1.d = 1.3e-3 / SI_unit_meter;
%! helspr1.D = 12.12e-3 / SI_unit_meter + helspr1.d;
%! helspr1.L = 27.7e-3 / SI_unit_meter;
%! helspr1.n = 5.7;
%! helspr1.ni = 2.7;
%! helspr1.ng = 0.75;
%! helspr1.m = 40;
%! helspr1.material.E = 206000e6 / SI_unit_pascal;
%! helspr1.material.G = 81500e6 / SI_unit_pascal;
%! [helspr1.material.alpha, helspr1.material.beta] = fem_pre_mat_rayleigh_damping([1e-2; 0.03e-2], [20, 1000] / (SI_unit_second^-1));
%! helspr1.material.nu = helspr1.material.E / (2 * helspr1.material.G) - 1;
%! helspr1.material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! helspr1.X = [[45.20e-3,  45.20e-3, -62.5e-3;
%!               43.13e-3, -43.13e-3,   0.0e-3] / SI_unit_meter;
%!              [-helspr1.L, -helspr1.L, -helspr1.L]];
%! section1.A = helspr1.d^2 * pi / 4;
%! section1.Ay = 9 / 10 * section1.A;
%! section1.Az = section1.Ay;
%! section1.Iy = helspr1.d^4 * pi / 64;
%! section1.Iz = section1.Iy;
%! section1.It = section1.Iy + section1.Iz;
%! rubspr1.d = 11.1e-3 / SI_unit_meter;
%! rubspr1.D = 28.6e-3 / SI_unit_meter;
%! rubspr1.L = 10.5e-3 / SI_unit_meter;
%! rubspr1.material.E = 2.6e6 / SI_unit_pascal;
%! rubspr1.material.nu = 0.499;
%! [rubspr1.material.alpha, rubspr1.material.beta] = fem_pre_mat_rayleigh_damping([10e-2; 30e-2], [30, 100] / (SI_unit_second^-1));
%! rubspr1.material.rho = 910 / (SI_unit_kilogram / SI_unit_meter^3);
%! rubspr1.X = [(170e-3 * [0.5, 0.5, -0.5, -0.5] + 8e-3) / SI_unit_meter;
%!              70e-3 * [-0.5, 0.5, -0.5, 0.5] / SI_unit_meter;
%!              -[1, 1, 1, 1] * (helspr1.L + rubspr1.L + param.offset1)];
%! rubspr1.e2 = [1, 0, 0];
%! section2.A = (rubspr1.D^2 - rubspr1.d^2) * pi / 4;
%! section2.Ay = 0.8 * section2.A;
%! section2.Az = section2.Ay;
%! section2.Iy = (rubspr1.D^4 - rubspr1.d^4) * pi / 64;
%! section2.Iz = section2.Iy;
%! section2.It = section2.Iy + section2.Iz;
%! node_idx_rb1 = int32(1);
%! node_idx_rb2 = int32(2);
%! mesh.nodes(node_idx_rb1, :) = [0, 0, rubspr1.L + param.offset1 + helspr1.L, zeros(1, 3)];
%! mesh.nodes(node_idx_rb2, :) = [0, 0, rubspr1.L, zeros(1, 3)];
%! node_idx_helspr1 = int32(rows(mesh.nodes) + 1);
%! helspr1.Phi = linspace(0, 2 * pi * helspr1.n, ceil(helspr1.m * helspr1.n)).' + 2 * pi * helspr1.ni;
%! node_idx_rubspr1 = node_idx_helspr1 + numel(helspr1.Phi) * columns(helspr1.X);
%! helspr1.x = 0.5 * helspr1.D * cos(helspr1.Phi);
%! helspr1.y = 0.5 * helspr1.D * sin(helspr1.Phi);
%! helspr1.z = (helspr1.L - helspr1.d * (2 * (helspr1.ni - helspr1.ng) + 1)) * linspace(0, 1, numel(helspr1.Phi))(:) + helspr1.d * (helspr1.ni - helspr1.ng + 0.5);
%! helspr1.e2 = repmat([0, 0, 1], numel(helspr1.Phi) - 1, 1);
%! idx_beam = int32(0);
%! for i=1:columns(helspr1.X)
%!   elnodes = int32([(1:numel(helspr1.Phi) - 1).', (2:numel(helspr1.Phi)).']) + node_idx_helspr1 - 1 + (i - 1) * numel(helspr1.Phi);
%!   mesh.nodes(node_idx_helspr1 - 1 + (1:numel(helspr1.Phi)) + (i - 1) * numel(helspr1.Phi), 1:6) = [[helspr1.x, helspr1.y, helspr1.z] + helspr1.X(:, i).', zeros(numel(helspr1.Phi), 3)];
%!   mesh.elements.beam2(idx_beam + (1:numel(helspr1.Phi) - 1)) = struct("nodes", mat2cell(elnodes, ones(numel(helspr1.Phi) - 1, 1, "int32"), 2), ...
%!                                                                       "section", mat2cell(repmat(section1, numel(helspr1.Phi) - 1, 1), ones(numel(helspr1.Phi) - 1, 1, "int32")), ...
%!                                                                       "e2", mat2cell(helspr1.e2, ones(numel(helspr1.Phi) - 1, 1, "int32"), 3));
%!   idx_beam += numel(helspr1.Phi) - 1;
%! endfor
%! for i=1:columns(rubspr1.X)
%!   elnodes = int32([1, 2]) + node_idx_rubspr1 - 1 + (i - 1) * 2;
%!   mesh.nodes(node_idx_rubspr1 - 1 + (1:2) + 2 * (i - 1), 1:6) = [[zeros(1, 3); zeros(1, 2), rubspr1.L] + rubspr1.X(:, i).', zeros(2, 3)];
%!   mesh.elements.beam2(idx_beam + 1) = struct("nodes", mat2cell(elnodes, ones(1, 1, "int32"), 2), ...
%!                                              "section", mat2cell(repmat(section2, 1, 1), ones(1, 1, "int32")), ...
%!                                              "e2", mat2cell(rubspr1.e2, ones(1, 1, "int32"), 3));
%!   ++idx_beam;
%! endfor
%! empty_cell = cell(1, 2);
%! mesh.elements.bodies = struct("m", empty_cell, "J", empty_cell, "lcg", empty_cell, "nodes", empty_cell);
%! mesh.elements.bodies(1).m = param.m1;
%! mesh.elements.bodies(1).J = param.J1;
%! mesh.elements.bodies(1).lcg = param.lcg1;
%! mesh.elements.bodies(1).nodes = node_idx_rb1;
%! mesh.elements.bodies(2).m = param.m2;
%! mesh.elements.bodies(2).J = param.J2;
%! mesh.elements.bodies(2).lcg = param.lcg2;
%! mesh.elements.bodies(2).nodes = node_idx_rb2;
%! empty_cell = cell(1, 2);
%! mesh.material_data = struct("E", empty_cell, "nu", empty_cell, "rho", empty_cell, "c", empty_cell);
%! mesh.material_data(1).E = helspr1.material.E;
%! mesh.material_data(1).nu = helspr1.material.nu;
%! mesh.material_data(1).beta = helspr1.material.beta;
%! mesh.material_data(1).alpha = helspr1.material.alpha;
%! mesh.material_data(1).rho = helspr1.material.rho;
%! mesh.material_data(2).E = rubspr1.material.E;
%! mesh.material_data(2).nu = rubspr1.material.nu;
%! mesh.material_data(2).beta = rubspr1.material.beta;
%! mesh.material_data(2).alpha = rubspr1.material.alpha;
%! mesh.material_data(2).rho = rubspr1.material.rho;
%! load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%! for i=1:columns(rubspr1.X)
%!   load_case_dof.locked_dof(node_idx_rubspr1 + 2 * (i - 1), 1:6) = true;
%! endfor
%! load_case_dof.domain = FEM_DO_STRUCTURAL;
%! mesh.materials.beam2 = [repmat(int32(1), columns(helspr1.X) * (numel(helspr1.Phi) - 1), 1);
%!                         repmat(int32(2), columns(rubspr1.X), 1)];
%! empty_cell = cell(1, numel(2 * columns(helspr1.X) + columns(rubspr1.X)));
%! mesh.elements.joints = struct("nodes", empty_cell, "C", empty_cell);
%! idx_joint = int32(0);
%! for i=1:columns(helspr1.X)
%!   node_idx_slave = node_idx_helspr1 - 1 + numel(helspr1.Phi) * i;
%!   lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb1, 1:3);
%!   ++idx_joint;
%!   mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb1, node_idx_slave]);
%!   mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                             zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%! endfor
%! for i=1:columns(helspr1.X)
%!   node_idx_slave = node_idx_helspr1 + numel(helspr1.Phi) * (i - 1);
%!   lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!   ++idx_joint;
%!   mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!   mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                             zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%! endfor
%! for i=1:columns(rubspr1.X)
%!   node_idx_slave = node_idx_rubspr1 + 2 * i - 1;
%!   lslave = mesh.nodes(node_idx_slave, 1:3) - mesh.nodes(node_idx_rb2, 1:3);
%!   ++idx_joint;
%!   mesh.elements.joints(idx_joint).nodes = int32([node_idx_rb2, node_idx_slave]);
%!   mesh.elements.joints(idx_joint).C = [     eye(3), -skew(lslave),     -eye(3), zeros(3, 3);
%!                                             zeros(3, 3),        eye(3), zeros(3, 3),     -eye(3)];
%! endfor
%! empty_cell = cell(1, 6);
%! load_case = struct("nodes", empty_cell, "loaded_nodes", empty_cell);
%! for i=1:numel(load_case)
%!   load_case(i).loaded_nodes = int32(node_idx_rb1);
%!   load_case(i).loads = zeros(1, 6);
%!   load_case(i).loads(i) = 1;
%! endfor
%! dof_map = fem_ass_dof_map(mesh, load_case_dof);
%! [mat_ass.M, ...
%!  mat_ass.K, ...
%!  mat_ass.D, ...
%!  mat_ass.R, ...
%!  mat_ass.mat_info, ...
%!  mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_MASS, ...
%!                                       FEM_MAT_STIFFNESS, ...
%!                                       FEM_MAT_DAMPING, ...
%!                                       FEM_VEC_LOAD_CONSISTENT], ...
%!                                      load_case);
%! opt_linsol.solver = "pastix";
%! opt_linsol.pre_scaling = true;
%! opt_linsol.number_of_threads = int32(4);
%! opt_linsol.refine_max_iter = int32(50);
%! opt_linsol.epsilon_refinement = 1e-10;
%! opt_linsol.verbose = int32(0);
%! opt_eig.p = 5 * param.N;
%! opt_eig.maxit = int32(100);
%! [sol_eig, Phi] = fem_sol_modal_damped(mesh, dof_map, mat_ass, param.N, opt_linsol, opt_eig);
%! [Phi, h] = fem_sol_modes_scale(mat_ass.M, mat_ass.K, sol_eig.lambda, Phi, mat_ass.R);
%! omega = linspace(2 * pi * param.fmin, 2 * pi * param.fmax, param.num_freq_modal);
%! U2 = fem_sol_harmonic_modal(h, sol_eig.lambda, Phi(dof_map.ndof(node_idx_rb2, :), :), omega);
%! for i=1:size(sol_eig.def, 3)
%!   sre = max(max(abs(real(sol_eig.def(:, 1:3, i)))));
%!   sim = max(max(abs(imag(sol_eig.def(:, 1:3, i)))));
%!   if (sim > sre)
%!     sol_eig.def(:, :, i) = imag(sol_eig.def(:, :, i)) * param.maxdef / sim;
%!   else
%!     sol_eig.def(:, :, i) = real(sol_eig.def(:, :, i)) * param.maxdef / sre;
%!   endif
%! endfor
%! if (param.num_freq_nodal > 0)
%! omegan = linspace(2 * pi * param.fmin, 2 * pi * param.fmax, param.num_freq_nodal);
%! U2n = zeros(6, numel(omegan), columns(mat_ass.R));
%! for i=1:numel(omegan)
%!  Un = fem_sol_factor(-omegan(i)^2 * mat_ass.M + 1j * omegan(i) * mat_ass.D + mat_ass.K, opt_linsol) \ mat_ass.R;
%!  for j=1:columns(Un)
%!    U2n(:, i, j) = Un(dof_map.ndof(node_idx_rb2, :), j);
%!  endfor
%! endfor
%! endif
%! figure("visible", "off");
%! hold on;
%! semilogy(sol_eig.f * SI_unit_second^-1, 100 * sol_eig.D, "-;D(f);r");
%! xlabel("f [Hz]");
%! ylabel("D [%]");
%! title("modal damping ratio");
%! grid on;
%! grid minor on;
%! for i=1:3
%!   figure("visible", "off");
%!   subplot(2, 1, 1);
%!   hold on;
%!   plot(omega / (2*pi) * (SI_unit_second^-1), 20 * log10(abs(-omega.^2 .* U2(i, :, i)) * (SI_unit_meter / SI_unit_second^2)), "-;modal;r");
%!   if (param.num_freq_nodal > 0)
%!     plot(omegan / (2*pi) * (SI_unit_second^-1), 20 * log10(abs(-omegan.^2 .* U2n(i, :, i)) * (SI_unit_meter / SI_unit_second^2)), "-;nodal;b");
%!   endif
%!   xlabel("f [Hz]");
%!   ylabel("a [dB/(1 m/s^2/N)]");
%!   ylim([-80, 10]);
%!   yticks(-100:5:10);
%!   grid on;
%!   grid minor on;
%!   title("frequency response magnitude");
%!   subplot(2, 1, 2);
%!   hold on;
%!   plot(omega / (2*pi) * (SI_unit_second^-1), 180 / pi * arg(-omega.^2 .* U2(i, :, i)), "-;modal;r");
%!   if (param.num_freq_nodal > 0)
%!     plot(omegan / (2*pi) * (SI_unit_second^-1), 180 / pi * arg(-omegan.^2 .* U2n(i, :, i)), "-;nodal;b");
%!   endif
%!   xlabel("f [Hz]");
%!   ylabel("arg(a) [deg]");
%!   ylim([-180,180]);
%!   yticks(-180:30:180);
%!   grid on;
%!   grid minor on;
%!   title("frequency response phase");
%! endfor
%! opt_post.elem_types = {"beam2"};
