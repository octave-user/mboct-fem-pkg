## Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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
        unlink(geo_file_tmp);
      endif
      if (numel(mesh_file))
        unlink(mesh_file);
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
