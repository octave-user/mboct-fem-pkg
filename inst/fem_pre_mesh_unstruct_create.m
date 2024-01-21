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

