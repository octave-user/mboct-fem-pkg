## Copyright (C) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_post_sol_plot(@var{mesh})
## @deftypefnx {} fem_post_sol_plot(@var{mesh}, @var{sol})
## @deftypefnx {} fem_post_sol_plot(@var{mesh}, @var{sol}, @var{scale})
## @deftypefnx {} fem_post_sol_plot(@var{mesh}, @var{sol}, @var{scale}, @var{idx_sol})
## @deftypefnx {} fem_post_sol_plot(@var{mesh}, @var{sol}, @var{scale}, @var{idx_sol}, @var{options})
##
## Plot a mesh in deformed or undeformed state.
##
## @var{mesh} @dots{} Mesh data structure returned from fem_pre_mesh_struct_create or fem_pre_mesh_unstruct_create.
##
## @var{sol} @dots{} Solution data structure returned from fem_sol_static.
##
## @var{scale} @dots{} Scale factor for deformations.
##
## @var{idx_sol} @dots{} Index of the solution if more than one solution exists in <@var{sol}>.
##
## @var{options}.node_labels @dots{} Plot node labels.
##
## @var{options}.elem_groups @dots{} Plot only selected groups of elements.
##
## @var{options}.elem_types @dots{} Plot only selected types of elements.
##
## @end deftypefn

function fem_post_sol_plot(mesh, sol, scale, idx_sol, options)
  if (nargin < 1 || nargin > 5 || nargout > 0)
    print_usage();
  endif

  if (nargin < 3)
    scale = 1;
  endif

  if (nargin < 4)
    idx_sol = 1;
  endif

  if (nargin < 5)
    options = struct();
  endif

  if (~isfield(options, "node_labels"))
    options.node_labels = false;
  endif

  if (~isfield(options, "elem_types"))
    options.elem_types = {"iso8", "iso20", "tet10", "rbe3", "tria6", "iso4", "quad8", "tria3"};
  endif

  if (~isfield(options, "elem_groups"))
    options.elem_groups = struct();
  endif

  nodes = mesh.nodes(:, 1:3);

  if (nargin >= 2 && numel(sol))
    nodes += scale * sol.def(:, 1:3, idx_sol);
  endif

  inumfaces = int32(0);

  for etype_idx=1:numel(options.elem_types)
    if (~isfield(mesh.elements, options.elem_types{etype_idx}))
      continue
    endif

    elements = getfield(mesh.elements, options.elem_types{etype_idx});

    switch (options.elem_types{etype_idx})
      case "iso8"
        inumfaces += 12 * rows(elements);
      case "iso20"
        inumfaces += 36 * rows(elements);	
      case "iso4"
        inumfaces += 2 * rows(elements);
      case "quad8"
	inumfaces += 6 * rows(elements);
      case "tet10"
        inumfaces += 16 * rows(elements);
      case "tria6"
        inumfaces += 4 * rows(elements);
      case "iso4"
        inumfaces += 2 * rows(elements);
      case "tria3"
        inumfaces += rows(elements);
      case "rbe3"
      otherwise
        error("unknown element type with %d nodes", columns(elements));
    endswitch
  endfor

  face_data = zeros(inumfaces, 3, "int32");

  inumfaces = int32(0);

  hold on;

  for etype_idx=1:length(options.elem_types)
    if (~isfield(mesh.elements, options.elem_types{etype_idx}))
      continue
    endif

    elements = getfield(mesh.elements, options.elem_types{etype_idx});

    switch (options.elem_types{etype_idx})
      case "iso8"
        faces = [1,2,3;
                 1,4,3;
                 1,2,6;
                 1,5,6;
                 1,4,5;
                 4,8,5;
                 4,3,8;
                 3,7,8;
                 3,2,6,
                 3,7,6,
                 6,7,8,
                 6,5,8];
      case "iso20"
	faces = [1, 9, 12;
		 2, 10, 9;
		 3, 11, 10;
		 4, 12, 11;
		 10, 11, 12;
		 9, 10, 12;
		 5, 13, 16;
		 6, 14, 13;
		 7, 15, 14;
		 8, 16, 15;
		 13, 14, 16;
		 14, 15, 16;
		 5, 13, 17;
		 6, 18, 13;
		 2, 9, 18;
		 1, 17, 9;
		 13, 18, 17;
		 9, 17, 18;
		 3, 11, 19;
		 4, 20, 11;
		 8, 15, 20;
		 7, 19, 15;
		 11, 20, 15;
		 11, 15, 19;
		 1, 12, 17;
		 4, 20, 12;
		 8, 16, 20;
		 5, 17, 16;
		 12, 16, 17;
		 12, 20, 16;
		 2, 10, 18;
		 3, 19, 10;
		 7, 14, 19;
		 6, 18, 14;
		 10, 19, 14;
		 10, 14, 18];
      case "iso4"
        faces = [1, 2, 3;
                 3, 4, 1];
      case "quad8"
	faces = [1, 5, 8;
		 2, 6, 5;
		 3, 7, 6;
		 4, 8, 7;
		 5, 6, 7;
		 5, 7, 8];
      case "tet10"
        faces = [1, 5, 8; ## 1
                 5, 2, 8;
                 2, 9, 8;
                 8, 9, 4;
                 2, 6, 9; ## 2
                 9, 10, 4;
                 6, 3, 10;
                 6, 9, 10;
                 3, 7, 10; ## 3
                 7, 8, 10;
                 7, 1, 8;
                 8, 4, 10;
                 1, 5, 7; ## 4
                 5, 2, 6;
                 5, 6, 7;
                 6, 3, 7];
      case "tria6"
        faces = [1, 4, 5;
                 2, 5, 4;
                 1, 5, 6;
                 3, 6, 5];
      case "tria3"
        faces = [1, 2, 3];
      case "rbe3"
        for elem_idx=1:numel(elements)
          for node_idx=2:numel(elements(elem_idx).nodes)
            xi_b = nodes(elements(elem_idx).nodes([1, node_idx]), 1:3);
            hnd = plot3(xi_b(:, 1), xi_b(:, 2), xi_b(:, 3));
            set(hnd, "color", zeros(1, 3));
            set(hnd, "linewidth", 3);
          endfor
        endfor
        continue
      otherwise
        error("unknown element type with %d nodes", columns(elements));
    endswitch

    if (isfield(options.elem_groups, options.elem_types{etype_idx}))
      elem_groups = getfield(mesh.groups, options.elem_types{etype_idx});
            
      for group_id=getfield(options.elem_groups, options.elem_types{etype_idx})
        for elem_idx=[[elem_groups(find([elem_groups.id] == group_id))].elements]
          face_data(inumfaces + (1:rows(faces)), 1:3) = elements(elem_idx, :)(faces);
          inumfaces += rows(faces);
        endfor
      endfor
    else
      for elem_idx=1:rows(elements)
        face_data(inumfaces + (1:rows(faces)), 1:3) = elements(elem_idx, :)(faces);
        inumfaces += rows(faces);
      endfor
    endif
  endfor

  if (nargin >= 2 && numel(sol))
    color_table = jet(1024);
    norm_U = norm(sol.def(:, 1:3, idx_sol), "rows");
    max_U = max(norm_U);
    grid_U = linspace(0, max_U, rows(color_table)).';

    for color_idx=1:columns(color_table)
      color_data(:, color_idx) = interp1(grid_U, color_table(:, color_idx), norm_U, "linear");
    endfor
  else
    color_data = [0, 0, 1];
  endif

  if (inumfaces)
    patch('Faces', face_data(1:inumfaces, :), 'Vertices', nodes, 'FaceVertexCData', color_data, 'FaceColor', 'interp');
  endif
  
  if (options.node_labels)
    for node_idx=1:rows(nodes)
      text(nodes(node_idx, 1), nodes(node_idx, 2), nodes(node_idx, 3), sprintf('n%d',node_idx));
    endfor
  endif

  daspect([1, 1, 1]);
endfunction

%!test ##demo
%! ## DEMO 1
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! h = 10e-3;
%! geometry.l = 60e-3;
%! geometry.w = 20e-3;
%! geometry.h = 50e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [1; 0; 0];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! figure("visible", "off");
%! fem_post_sol_plot(mesh);
%! view(30, 30);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! view(15, 15);
%! grid on;
%! grid minor on;
%! title("undeformed cube mesh");
%! figure_list();

%!test ##demo
%! ## DEMO 2
%! close all;
%! material.E = 210000e6;
%! material.nu = 0.3;
%! material.rho = 7850;
%! h = 5e-3;
%! scale = 20e-3;
%! geometry.l = 60e-3;
%! geometry.w = 20e-3;
%! geometry.h = 50e-3;
%! mesh_size.num_elem_l = ceil(geometry.l / h);
%! mesh_size.num_elem_w = ceil(geometry.w / h);
%! mesh_size.num_elem_h = ceil(geometry.h / h);
%! f = [0; 0; -1000];
%! [mesh, load_case] = fem_pre_mesh_cube_create(geometry, mesh_size, material, f);
%! dof_map = fem_ass_dof_map(mesh, load_case);
%! [mat_ass.K, mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_MAT_STIFFNESS, ...
%!                                          FEM_VEC_LOAD_CONSISTENT], ...
%!                                         load_case);
%! sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!
%! figure("visible", "off");
%! fem_post_sol_plot(mesh, sol_stat, scale / max(max(abs(sol_stat.def))));
%! view(15, 15);
%! xlabel("x [m]");
%! ylabel("y [m]");
%! zlabel("z [m]");
%! grid on;
%! grid minor on;
%! title("deformed cube mesh");
%! figure_list();
