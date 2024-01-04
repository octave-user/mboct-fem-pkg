## Copyright (C) 2018(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn{Function File} [@var{mesh_split}] = fem_pre_mesh_elem_split(@var{mesh})
## @deftypefnx{} [@dots{}] = fem_pre_mesh_elem_split(@var{mesh}, @var{options})
## Convert eight node hexahedra to six node pentahedra by splitting,
## if the maximum corner angle is above certain limit.
##
## @var{mesh} @dots{} Original mesh with distorted elements.
##
## @var{options}.iso8.max_angle @dots{} Limiting value for corner angles which causes splitting.
##
## @var{mesh_split} @dots{} New mesh with reduced corner angles.
## @seealso{fem_pre_mesh_elem_angle}
## @end deftypefn

function mesh_split = fem_pre_mesh_elem_split(mesh, options)
  if (nargin < 1 || nargin > 2 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "iso8"))
    options.iso8 = struct();
  endif

  if (~isfield(options.iso8, "max_angle"))
    options.iso8.max_angle = 175 * pi / 180;
  endif

  mesh_split = mesh;
  angle_data = fem_pre_mesh_elem_angle(mesh);

  cond_split = any(mean(angle_data.iso8, 3) > options.iso8.max_angle, 2);
  idx_elem_keep = find(~cond_split);
  idx_elem_split = find(cond_split);

  inum_elem = numel(idx_elem_keep);
  mesh_split.elements.iso8 = zeros(numel(idx_elem_keep) + 2 * numel(idx_elem_split), columns(mesh.elements.iso8), "int32");
  mesh_split.materials.iso8 = zeros(numel(idx_elem_keep) + 2 * numel(idx_elem_split), 1, "int32");

  if (numel(idx_elem_keep))
    mesh_split.elements.iso8(1:inum_elem, :) = mesh.elements.iso8(idx_elem_keep, :);
    mesh_split.materials.iso8(1:inum_elem, :) = mesh.materials.iso8(idx_elem_keep);
  endif

  for i=1:numel(idx_elem_split)
    Phi = mean(angle_data.iso8(idx_elem_split(i), :, :), 3);

    idx_max_angle = find(Phi == max(Phi));

    if (numel(idx_max_angle))
      switch (idx_max_angle(1))
        case {1, 3}
          idx_nodes = [1, 2, 3, 3, 5, 6, 7, 7;
                       1, 3, 4, 4, 5, 7, 8, 8];
        case {2, 4}
          idx_nodes = [2, 3, 4, 4, 6, 7, 8, 8;
                       1, 2, 4, 4, 5, 6, 8, 8];
        case {5, 7}
          idx_nodes = [5, 6, 2, 2, 8, 7, 3, 3;
                       5, 2, 1, 1, 8, 3, 4, 4];
        case {6, 8}
          idx_nodes = [6, 2, 1, 1, 7, 3, 4, 4;
                       5, 6, 1, 1, 8, 7, 4, 4];
        case {9, 11}
          idx_nodes = [5, 1, 4, 4, 6, 2, 3, 3;
                       5, 4, 8, 8, 6, 3, 7, 7];
        case {10, 12}
          idx_nodes = [1, 4, 8, 8, 2, 3, 7, 7;
                       1, 8, 5, 5, 2, 7, 6, 6];
        otherwise
          error("invalid node index %d detected", idx_max_angle(1));
      endswitch

      for j=1:rows(idx_nodes)
        mesh_split.elements.iso8(++inum_elem, :) = mesh.elements.iso8(idx_elem_split(i), idx_nodes(j, :));
        mesh_split.materials.iso8(inum_elem) = mesh.materials.iso8(idx_elem_split(i));
      endfor
      continue
    endif
  endfor
endfunction

%!test
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 120e-3;
%! rho = 7850;
%! m = 3 / 8 * a * b * c * rho;
%! tol = eps^0.9;
%!
%! X = [0.5 * a, 0.5 * b, c;
%!            0, 0.5 * b, c;
%!            0,       0, c;
%!            a,       0, c;
%!      0.5 * a, 0.5 * b, 0;
%!            0, 0.5 * b, 0;
%!            0,       0, 0;
%!            a,       0, 0];
%! Phi1 = [0, 20, 45, 330] * pi / 180;
%! Phi2 = [0, 60, 270, 285] * pi / 180;
%! Phi3 = [0, -15, 30, -195] * pi / 180;
%! options.iso8.max_angle = 95 * pi / 180;
%! rand("seed", 0);
%! idx_elem = int32([1:8;
%!                   5, 1, 4, 8, 6, 2, 3, 7;
%!                   4, 3, 7, 8, 1, 2, 6, 5;
%!                   3, 2, 6, 7, 4, 1, 5, 8;
%!                   1, 4, 8, 5, 2, 3, 7, 6;
%!                   1, 5, 6, 2, 4, 8, 7, 3;
%!                   3, 2, 6, 7, 4, 1, 5, 8;
%!                   8, 7, 6, 5, 4, 3, 2, 1;
%!                   7, 6, 5, 8, 3, 2, 1, 4;
%!                   5, 8, 7, 6, 1, 4, 3, 2]);
%! idx_mesh = int32(0);
%! for i=1:rows(idx_elem)
%! for j=1:length(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   mesh(++idx_mesh).nodes = [X * R1.', zeros(rows(X), 3)];
%!   mesh(idx_mesh).elements.iso8 = idx_elem(i, :);
%!   mesh(idx_mesh).materials.iso8 = int32([1]);
%!   load_case(idx_mesh).locked_dof = false(rows(mesh(idx_mesh).nodes), 6);
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh(idx_mesh).material_data.rho = rho;
%!   mesh(idx_mesh).material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map(idx_mesh) = fem_ass_dof_map(mesh(idx_mesh), load_case(idx_mesh));
%!   [rndval, idx_node] = sort(rand(rows(mesh(idx_mesh).nodes), 1));
%!   idx_node_inv(idx_node) = 1:rows(mesh(idx_mesh).nodes);
%!   for k=1:rows(mesh(j).elements.iso8)
%!     mesh(idx_mesh).elements.iso8(k, :) = idx_node_inv(mesh(idx_mesh).elements.iso8(k, :));
%!   endfor
%!   mesh(idx_mesh).nodes = mesh(idx_mesh).nodes(idx_node, :);
%!   mesh_split(idx_mesh) = fem_pre_mesh_elem_split(mesh(idx_mesh), options);
%!   mesh_split(idx_mesh).nodes(:, 3) += d;
%!   mat_ass(idx_mesh).m = fem_ass_matrix(mesh(idx_mesh), ...
%!                                        dof_map(idx_mesh), ...
%!                                        [FEM_SCA_TOT_MASS], ...
%!                                        load_case(idx_mesh));
%!   mat_ass_split(idx_mesh).m = fem_ass_matrix(mesh_split(idx_mesh), ...
%!                                              dof_map(idx_mesh), ...
%!                                              [FEM_SCA_TOT_MASS], ...
%!                                              load_case(idx_mesh));
%!   if (do_plot)
%!   figure("visible","off");
%!   fem_post_sol_plot(mesh(idx_mesh));
%!   fem_post_sol_plot(mesh_split(idx_mesh));
%!   title(sprintf("mesh %d: elements before %d, elements after %d", j, rows(mesh(idx_mesh).elements.iso8), rows(mesh_split(idx_mesh).elements.iso8)));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   endif
%! endfor
%! endfor
%! assert_simple(all(abs([mat_ass.m] - m) < tol * m));
%! assert_simple(all(abs([mat_ass_split.m] - m) < tol * m));
%! if (do_plot)
%! figure_list();
%! endif

%!test
%! do_plot = false;
%! if (do_plot)
%! close all;
%! endif
%! a = 70e-3;
%! b = 20e-3;
%! c = 10e-3;
%! d = 120e-3;
%! rho = 7850;
%! m = 3 / 8 * a * b * c * rho;
%! tol = eps^0.9;
%!
%! X = [      a,       b, c;
%!            0,       b, c;
%!            0, 0.5 * b, c;
%!      0.5 * a, 0.5 * b, c;
%!            a,       b, 0;
%!            0,       b, 0;
%!            0, 0.5 * b, 0;
%!      0.5 * a, 0.5 * b, 0];
%!
%! Phi1 = [0, 20, 45, 330] * pi / 180;
%! Phi2 = [0, 60, 270, 285] * pi / 180;
%! Phi3 = [0, -15, 30, -195] * pi / 180;
%! options.iso8.max_angle = 95 * pi / 180;
%! rand("seed", 0);
%! idx_elem = int32([1:8;
%!                   5, 1, 4, 8, 6, 2, 3, 7;
%!                   4, 3, 7, 8, 1, 2, 6, 5;
%!                   3, 2, 6, 7, 4, 1, 5, 8;
%!                   1, 4, 8, 5, 2, 3, 7, 6;
%!                   1, 5, 6, 2, 4, 8, 7, 3;
%!                   3, 2, 6, 7, 4, 1, 5, 8;
%!                   8, 7, 6, 5, 4, 3, 2, 1;
%!                   7, 6, 5, 8, 3, 2, 1, 4;
%!                   5, 8, 7, 6, 1, 4, 3, 2]);
%! idx_mesh = int32(0);
%! for i=1:rows(idx_elem)
%! for j=1:length(Phi1)
%!   R1 = euler123_to_rotation_matrix([Phi1(j); Phi2(j); Phi3(j)]);
%!   mesh(++idx_mesh).nodes = [X * R1.', zeros(rows(X), 3)];
%!   mesh(idx_mesh).elements.iso8 = idx_elem(i, :);
%!   mesh(idx_mesh).materials.iso8 = int32([1]);
%!   load_case(idx_mesh).locked_dof = false(rows(mesh(idx_mesh).nodes), 6);
%!   E = 210000e6;
%!   nu = 0.3;
%!   mesh(idx_mesh).material_data.rho = rho;
%!   mesh(idx_mesh).material_data.C = fem_pre_mat_isotropic(E, nu);
%!   dof_map(idx_mesh) = fem_ass_dof_map(mesh(idx_mesh), load_case(idx_mesh));
%!   [rndval, idx_node] = sort(rand(rows(mesh(idx_mesh).nodes), 1));
%!   idx_node_inv(idx_node) = 1:rows(mesh(idx_mesh).nodes);
%!   for k=1:rows(mesh(j).elements.iso8)
%!     mesh(idx_mesh).elements.iso8(k, :) = idx_node_inv(mesh(idx_mesh).elements.iso8(k, :));
%!   endfor
%!   mesh(idx_mesh).nodes = mesh(idx_mesh).nodes(idx_node, :);
%!   mesh_split(idx_mesh) = fem_pre_mesh_elem_split(mesh(idx_mesh), options);
%!   mesh_split(idx_mesh).nodes(:, 3) += d;
%!   mat_ass(idx_mesh).m = fem_ass_matrix(mesh(idx_mesh), ...
%!                                        dof_map(idx_mesh), ...
%!                                        [FEM_SCA_TOT_MASS], ...
%!                                        load_case(idx_mesh));
%!   mat_ass_split(idx_mesh).m = fem_ass_matrix(mesh_split(idx_mesh), ...
%!                                              dof_map(idx_mesh), ...
%!                                              [FEM_SCA_TOT_MASS], ...
%!                                              load_case(idx_mesh));
%!   if (do_plot)
%!   figure("visible","off");
%!   fem_post_sol_plot(mesh(idx_mesh));
%!   fem_post_sol_plot(mesh_split(idx_mesh));
%!   title(sprintf("mesh %d: elements before %d, elements after %d", j, rows(mesh(idx_mesh).elements.iso8), rows(mesh_split(idx_mesh).elements.iso8)));
%!   xlabel("x [m]");
%!   ylabel("y [m]");
%!   zlabel("z [m]");
%!   grid on;
%!   grid minor on;
%!   endif
%! endfor
%! endfor
%! assert_simple(all(abs([mat_ass.m] - m) < tol * m));
%! assert_simple(all(abs([mat_ass_split.m] - m) < tol * m));
%! if (do_plot)
%! figure_list();
%! endif

%!test
%! do_plot = false;
%! if (do_plot)
%! close all;
%! endif
%! a = 20e-3;
%! b = 15e-3;
%! c = 10e-3;
%! rho = 7850;
%! m = 3 / 8 * a * b * c * rho;
%! tol_m = eps^0.9;
%! N = 10;
%! r = 1;
%! X = [0.5 * a, 0.5 * b, c;
%!            0, 0.5 * b, c;
%!            0,       0, c;
%!            a,       0, c;
%!      0.5 * a, 0.5 * b, 0;
%!            0, 0.5 * b, 0;
%!            0,       0, 0;
%!            a,       0, 0];
%!
%! mesh_data(1).mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh_data(1).mesh.elements.iso8 = int32([1:8]);
%! mesh_data(1).mesh.materials.iso8 = int32([1]);
%! mesh_data(1).load_case.locked_dof = false(rows(mesh_data(1).mesh.nodes), 6);
%! E = 210000e6;
%! nu = 0.3;
%! mesh_data(1).mesh.material_data.rho = rho;
%! mesh_data(1).mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%! options.iso8.max_angle = 95 * pi / 180;
%! mesh_data(2).mesh = fem_pre_mesh_elem_split(mesh_data(1).mesh, options);
%! mesh_data(2).load_case = mesh_data(1).load_case;
%! for i=1:numel(mesh_data)
%!   mesh_data(i).dof_map = fem_ass_dof_map(mesh_data(i).mesh, mesh_data(i).load_case);
%!   [mesh_data(i).mat_ass.M, ...
%!    mesh_data(i).mat_ass.K, ...
%!    mesh_data(i).mat_ass.m] = fem_ass_matrix(mesh_data(i).mesh, ...
%!                                             mesh_data(i).dof_map, ...
%!                                             [FEM_MAT_MASS, ...
%!                                              FEM_MAT_STIFFNESS, ...
%!                                              FEM_SCA_TOT_MASS], ...
%!                                             mesh_data(i).load_case);
%!   mesh_data(i).sol_eig = fem_sol_modal(mesh_data(i).mesh, mesh_data(i).dof_map, mesh_data(i).mat_ass, N, r);
%! endfor
%! opts.print_and_exit = true;
%! opts.print_to_file = "";
%! unwind_protect
%!   opts.print_to_file = tempname();
%!   for i=1:numel(mesh_data)
%!     for j=8:9
%!       opts.scale_def = 2.5e-3 / max(max(abs(mesh_data(i).sol_eig.def(:, 1:3, j))));
%!       opts.output_step_idx = j;
%!       if (do_plot)
%!         fem_post_sol_external(mesh_data(i).mesh, mesh_data(i).sol_eig, opts);
%!         [img, map, alpha] = imread(sprintf("%s_%03d.jpg", opts.print_to_file, 1));
%!         figure("visible", "off");
%!         imshow(img, map);
%!         title(sprintf("mode %d f=%.0fHz", j, mesh_data(i).sol_eig.f(j)));
%!       endif
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (numel(opts.print_to_file))
%!     fn = dir([opts.print_to_file, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! for i=1:numel(mesh_data)
%!   assert_simple(mesh_data(i).mat_ass.m, m, tol_m * m);
%! endfor
%! if (do_plot)
%!   figure_list();
%! endif
