## Copyright (C) 2016(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_ehd_pre_comp_mat_plot(@var{comp_mat}, @var{options})
## Create contour plots of a compliance matrix for debugging purposes.
##
## @var{comp_mat} @dots{} Return value from fem_ehd_comp_mat_struct or fem_ehd_comp_mat_unstruct.
##
## @var{options}.nodal_plot @dots{} Create a contour plot of the radial deflection at each node.
##
## @var{options}.modal_plot @dots{} Create a contour plot of the radial deformation of each mode shape.
##
## @var{options}.plot_mesh @dots{} Create a plot of the deformed mesh.
##
## @var{options}.scale_modes @dots{} Scale factor for plots of the deformed mesh.
##
## @var{options}.plot_load_cases @dots{} Create a contour plot of the pressure distribution.
##
## @var{options}.merge_files @dots{} Create a single .pdf file for all figures.
##
## @var{options}.plot_matrix @dots{} Create a contour plot of each column of the compliance matrix.
##
## @var{options}.plot_const_pressure @dots{} Create a contour plot of the radial deflection of the bearing surface under a constant pressure load.
##
## @var{options}.output_file @dots{} Output filename of a .pdf file containing all the graphs.
##
## @end deftypefn

function fem_ehd_pre_comp_mat_plot(comp_mat, options)
  if (nargin ~= 2 || nargout > 0)
    print_usage();
  endif

  if (~isfield(options, "nodal_plot"))
    options.nodal_plot = true;
  endif

  if (~isfield(options, "modal_plot"))
    options.modal_plot = true;
  endif

  if (~isfield(options, "scale_modes"))
    options.scale_modes = 1;
  endif
  
  if (~isfield(options, "plot_mesh"))
    options.plot_mesh = true;
  endif

  if (~isfield(options, "plot_load_cases"))
    options.plot_load_cases = true;
  endif

  if (~isfield(options, "merge_files"))
    options.merge_files = true;
  endif

  if (~isfield(options, "plot_step"))
    options.plot_step = struct();
  endif

  if (~isfield(options.plot_step, "x"))
    options.plot_step.x = 1;
  endif

  if (~isfield(options.plot_step, "z"))
    options.plot_step.z = 1;
  endif

  if (isfield(options, "output_file"))
    unlink(options.output_file);
  endif

  if (~isfield(options, "contour_levels"))
    options.contour_levels = int32(50);
  endif

  is_struct_mesh = true;

  for i=1:numel(comp_mat)
    if (~(isfield(comp_mat(i), "mesh") && isfield(comp_mat(i).mesh, "structured")))
      is_struct_mesh = false;
      break;
    endif   
  endfor

  if (~isfield(options, "plot_matrix"))
    options.plot_matrix = ~is_struct_mesh;
  endif

  if (~isfield(options, "plot_const_pressure"))
    options.plot_const_pressure = ~is_struct_mesh;
  endif

  output_files = {};
  
  for i=1:numel(comp_mat)
    output_files = fem_ehd_comp_mat_plot_generic(comp_mat(i), output_files, options);
    
    if (is_struct_mesh)
      output_files = fem_ehd_comp_mat_plot_struct_mesh(comp_mat(i), output_files, options);
    endif
  endfor
  
  if (options.merge_files && isfield(options, "output_file"))
    pdf_merge(output_files, options.output_file);
  endif
endfunction

function output_files = fem_ehd_comp_mat_plot_struct_mesh(comp_mat, output_files, options)
  narginchk(3, 3);

  r = comp_mat.mesh.structured.r;
  s = comp_mat.mesh.structured.s;
  t = comp_mat.mesh.structured.t;
  Ur = zeros(numel(s), numel(t));
  Phi = zeros(1, numel(s));
  z = zeros(1, numel(t));
  X = zeros(3, numel(s));
  
  for j=1:numel(s)
    inode = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(1).position.r), j, find(t == comp_mat.loads(1).position.t));
    X(:, j) = comp_mat.mesh.nodes(inode, 1:3);
  endfor

  [rc, Xc] = fem_ehd_center_of_circle(X);
  
  Phi = atan2(X(2, 1) - Xc(2), X(1, 1) - Xc(1)) + comp_mat.bearing_surf.grid_x / rc;
  
  for k=1:numel(t)
    inode = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(1).position.r), 1, k);
    z(k) = comp_mat.mesh.nodes(inode, 3);
  endfor
  
  hnd = -1;
  
  if (options.nodal_plot && isfield(comp_mat, "sol_stat") && options.plot_load_cases)
    for i=1:size(comp_mat.sol_stat.def, 3)
      for j=1:numel(s)
        for k=1:numel(t)
          inode = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(i).position.r), j, k);
          n = comp_mat.mesh.nodes(inode, 1:2).' - Xc(1:2);
          Ur(j, k) = comp_mat.sol_stat.def(inode, 1:2, i) * (n / norm(n));
        endfor
      endfor
      inode_load = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(i).position.r), ...
                                                      comp_mat.loads(i).index.i, ...
                                                      comp_mat.loads(i).index.j ...
                                                      + comp_mat.bearing_surf.offset_t);

      hnd = figure("visible", "off");
      colormap("jet");
      contourf(180/pi*Phi, 1e3 * z, 1e6 * Ur.');
      colorbar();
      for j=1:numel(s)
        for k=1:numel(t)
          inode = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(i).position.r), j, k);
          text(180 / pi * Phi(j), 1e3 * z(k), sprintf("N%d", inode), "rotation", 90);
        endfor
      endfor
      xlabel("Phi [deg]");
      ylabel("z [mm]");
      grid on;
      grid minor on;
      title(sprintf("pressure at node %d = %g - radial deformation (%g:%g) [um]", ...
                    inode_load, ...
                    comp_mat.reference_pressure, ...
                    1e6 * min(min(Ur)), ...
                    1e6 * max(max(Ur))));

      if (isfield(options, "output_file"))
        output_files{end + 1} = gen_pdf_file(hnd, options, output_files);
      endif
    endfor
  endif
  
  if (options.modal_plot)
    if (isfield(comp_mat, "modal") && options.plot_mesh)
      for i=1:size(comp_mat.modal.def, 3)
        hnd = figure("visible", "off");
        fem_post_sol_plot(comp_mat.mesh, comp_mat.modal, options.scale_modes, i, options);
        xlabel("x [m]");
        ylabel("y [m]");
        zlabel("z [m]");
        grid on;
        grid minor on;
        title(sprintf("mode %d: %g", i, comp_mat.lambda(i)));

        if (isfield(options, "output_file"))
          output_files{end + 1} = gen_pdf_file(hnd, options, output_files);
        endif
      endfor
    endif
    
    if (isfield(comp_mat, "modal_sol") && options.plot_load_cases)
      U = comp_mat.modal_sol.def;
      
      for i=1:size(U, 3)
        inode_load = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(i).position.r), ...
                                                        comp_mat.loads(i).index.i, ...
                                                        comp_mat.loads(i).index.j ...
                                                        + comp_mat.bearing_surf.offset_t);
        
        for j=1:numel(s)
          for k=1:numel(t)
            inode = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(i).position.r), j, k);
            n = comp_mat.mesh.nodes(inode, 1:2).' - Xc(1:2);
            Ur(j, k) = U(inode, 1:2, i) * (n / norm(n));
          endfor
        endfor
        

        hnd = figure("visible", "off");
        contourf(180/pi * Phi, 1e3 * z, 1e6 * Ur.');
        colorbar();
        for j=1:numel(s)
          for k=1:numel(z)
            inode = comp_mat.mesh.structured.inode_idx(find(r == comp_mat.loads(i).position.r), j, k);
            text(180 / pi * Phi(j), ...
                 1e3 * z(k), ...
                 sprintf("N%d", inode), "rotation", 90);
          endfor
        endfor
        xlabel("Phi [deg]");
        ylabel("z [mm]");
        grid on;
        grid minor on;
        title(sprintf("pressure at node %d = %g - modal radial deformation (%g:%g) [um]", ...
                      inode_load, ...
                      comp_mat.reference_pressure, ...
                      1e6 * min(min(Ur)), ...
                      1e6 * max(max(Ur))));
        
        if (isfield(options, "output_file"))
          output_files{end + 1} = gen_pdf_file(hnd, options, output_files);
        endif
      endfor
    endif
  endif
endfunction

function output_files = fem_ehd_comp_mat_plot_generic(comp_mat, output_files, options)
  hnd = -1;
  
  if (options.plot_const_pressure && isfield(comp_mat, "C"))
    w = zeros(numel(comp_mat.bearing_surf.grid_x), numel(comp_mat.bearing_surf.grid_z));

    for l=1:rows(w)
      w(l, :) = sum(comp_mat.C((l - 1) * columns(w) + (1:columns(w)), 1:end - columns(w)), 2);
    endfor
    
    zrange = 1e6 * linspace(min(min(w)), max(max(w)), options.contour_levels);
    hnd = figure("visible", "off");
    colormap jet;
    contourf(1e3 * comp_mat.bearing_surf.grid_x, 1e3 * comp_mat.bearing_surf.grid_z, 1e6 * w.', zrange);
    xlabel("x [mm]");
    ylabel("z [mm]");
    grid on;
    grid minor on;
    colorbar();
    title(sprintf("constant pressure load case min(w)=%.2fum, max(w)=%.2fum, w [um]", ...
                  1e6 * min(min(w)), ...
                  1e6 * max(max(w))));
    
    if (isfield(options, "output_file"))
      output_files{end + 1} = gen_pdf_file(hnd, options, output_files);
    endif
  endif
  
  if (options.plot_matrix && isfield(comp_mat, "C"))
    for j=1:options.plot_step.x:numel(comp_mat.bearing_surf.grid_x)
      for k=1:options.plot_step.z:numel(comp_mat.bearing_surf.grid_z)
        w = zeros(numel(comp_mat.bearing_surf.grid_x), numel(comp_mat.bearing_surf.grid_z));
        for l=1:rows(w)
          w(l, :) = comp_mat.C((l - 1) * columns(w) + (1:columns(w)), (j - 1) * columns(w) + k);
        endfor
        zrange = 1e6 * linspace(min(min(w)), max(max(w)), options.contour_levels);
        hnd = figure("visible", "off");
        colormap jet;
        contourf(1e3 * comp_mat.bearing_surf.grid_x, 1e3 * comp_mat.bearing_surf.grid_z, 1e6 * w.', zrange);
        text(1e3 * comp_mat.bearing_surf.grid_x(j), 1e3 * comp_mat.bearing_surf.grid_z(k), sprintf("p(%d,%d)", j, k));
        xlabel("x [mm]");
        ylabel("z [mm]");
        grid on;
        grid minor on;
        colorbar();
        title(sprintf("load case (%d, %d) x=%.2fmm z=%.2fmm, min(w)=%.2fum, max(w)=%.2fum, w [um]", ...
                      j, ...
                      k, ...
                      1e3 * comp_mat.bearing_surf.grid_x(j), ...
                      1e3 * comp_mat.bearing_surf.grid_z(k), ...
                      1e6 * min(min(w)), ...
                      1e6 * max(max(w))));
        
        if (isfield(options, "output_file"))
          output_files{end + 1} = gen_pdf_file(hnd, options, output_files);
        endif
      endfor
    endfor

    if (isfield(comp_mat, "D"))
      for j=1:columns(comp_mat.D)
        w = zeros(numel(comp_mat.bearing_surf.grid_x), numel(comp_mat.bearing_surf.grid_z));
        for l=1:rows(w)
          w(l, :) = comp_mat.D((l - 1) * columns(w) + (1:columns(w)), j);
        endfor
        zrange = 1e6 * linspace(min(min(w)), max(max(w)), options.contour_levels);
        hnd = figure("visible", "off");
        colormap jet;
        contourf(1e3 * comp_mat.bearing_surf.grid_x, 1e3 * comp_mat.bearing_surf.grid_z, 1e6 * w.', zrange);
        xlabel("x [mm]");
        ylabel("z [mm]");
        grid on;
        grid minor on;
        colorbar();
        title(sprintf("contribution of mode (%d) min(w)=%.2fum, max(w)=%.2fum, w [um]", ...
                      j, ...
                      1e6 * min(min(w)), ...
                      1e6 * max(max(w))));
        
        if (isfield(options, "output_file"))
          output_files{end + 1} = gen_pdf_file(hnd, options, output_files);
        endif
      endfor
    endif
  endif
endfunction

function output_file = gen_pdf_file(hnd, options, output_files)
  [out_dir, out_name, out_ext] = fileparts(options.output_file);
  out_ext = ".pdf";

  output_file  = fullfile(out_dir, sprintf("%s_%03d%s", out_name, numel(output_files) + 1, out_ext));
  
  unlink(output_file);
  print("-dpdf", "-color", "-landscape", "-S1024,768", output_file);
  [info, err] = stat(output_file);
  
  if (err ~= 0 || info.size == 0)
    error("failed to create file \"%s\"", output_file);
  endif
  
  close(hnd);
endfunction

function [r, Xc, Rc] = fem_ehd_center_of_circle(X)
  X1 = X(:, 1);
  X2 = X(:, 2);
  X3 = X(:, 3);

  R1 = [0, -1, 0;
        1,  0, 0;
        0,  0, 1];

  n12 = (X2 - X1) / norm(X2 - X1);
  n23 = (X3 - X2) / norm(X3 - X2);

  nc = cross(n12, n23);

  nc /= norm(nc);

  Rc = [n12, n23, nc];

  n12 = Rc * R1 * Rc.' * n12;
  n23 = Rc * R1 * Rc.' * n23;
  
  X12 = 0.5 * (X1 + X2);
  X23 = 0.5 * (X2 + X3);

  t = line_intersection([X12, X23], [n12, n23]);
  
  Xc = t(1) * n12 + X12;

  r = norm(Xc - X1);
endfunction

function [t] = line_intersection(X0, n)
  A = [-n(:, 1).' * n(:, 1),  n(:, 1).' * n(:, 2);
        n(:, 2).' * n(:, 1), -n(:, 2).' * n(:, 2)];

  b = [ n(:, 1).' * (X0(:, 1) - X0(:, 2));
       -n(:, 2).' * (X0(:, 1) - X0(:, 2))];
  
  t = A \ b;
endfunction

%!demo
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! mesh_size = 1.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  error("gmsh failed with status %d", status);
%! endif

%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! node_idx_modal = rows(mesh.nodes) + 1;
%! node_idx_itf1 = rows(mesh.nodes) + 2;
%! node_idx_itf2 = rows(mesh.nodes) + 3;
%! mesh.nodes(node_idx_itf2, 1:3) = [0, 0, c + 0.5 * b];
%! mesh.nodes(node_idx_itf1, 1:3) = [0, 0, c + 0.5 * b];
%! mesh.nodes(node_idx_modal, :) = 0;
%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [grp_id_p1, grp_id_p2], [node_idx_itf1, node_idx_itf2], "tria6");
%! bearing_surf(1).group_idx = grp_id_p1;
%! bearing_surf(1).options.reference_pressure = 1e9;
%! bearing_surf(1).options.mesh_size = 1.5 * mesh_size;
%! bearing_surf(1).options.bearing_type = "shell";
%! bearing_surf(1).options.matrix_type = "nodal substruct";
%! bearing_surf(1).master_node_no = node_idx_itf1;
%! bearing_surf(1).nodes = mesh.groups.tria6(grp_id_p1).nodes;
%! bearing_surf(1).r = ri;
%! bearing_surf(1).w = b;
%! bearing_surf(1).X0 = [0; 0; b/2 + c];
%! bearing_surf(1).R = eye(3);
%! bearing_surf(1).relative_tolerance = 0;
%! bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! bearing_surf(2).group_idx = grp_id_p2;
%! bearing_surf(2).options.reference_pressure = 1e9;
%! bearing_surf(2).options.mesh_size = 1.5 * mesh_size;
%! bearing_surf(2).options.bearing_type = "journal";
%! bearing_surf(2).options.matrix_type = "nodal substruct";
%! bearing_surf(2).master_node_no = node_idx_itf2;
%! bearing_surf(2).nodes = mesh.groups.tria6(grp_id_p2).nodes;
%! bearing_surf(2).r = ro;
%! bearing_surf(2).w = b;
%! bearing_surf(2).X0 = [0; 0; b/2 + c];
%! bearing_surf(2).R = eye(3);
%! bearing_surf(2).relative_tolerance = 0;
%! bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! [load_case_press, bearing_surf] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! load_case = fem_pre_load_case_merge(load_case, load_case_press);
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! E = 210000e6;
%! nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! cms_opt.nodes.modal.number = node_idx_modal;
%! cms_opt.nodes.modal.name = "node_id_modal";
%! cms_opt.nodes.interfaces(1).number = node_idx_itf1;
%! cms_opt.nodes.interfaces(1).name = "node_id_itf1";
%! cms_opt.nodes.interfaces(2).number = node_idx_itf2;
%! cms_opt.nodes.interfaces(2).name = "node_id_itf2";
%! cms_opt.number_of_threads = int32(4);
%! [mesh, mat_ass, dof_map, sol_eig, cms_opt] = fem_cms_create(mesh, load_case, cms_opt);
%! comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%! opt_plot.plot_nodal = true;
%! opt_plot.contour_levels = int32(15);
%! opt_plot.plot_step.x = int32(2);
%! opt_plot.plot_step.z = int32(2);
%! fem_ehd_pre_comp_mat_plot(comp_mat, opt_plot);
%! figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
