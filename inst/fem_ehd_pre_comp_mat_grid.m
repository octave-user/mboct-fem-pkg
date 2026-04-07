function [bearing_surf, idx] = fem_ehd_pre_comp_mat_grid(mesh, bearing_surf, options)
  if (~isfield(options, "elem_type"))
    options.elem_type = "tria6";
  endif

  if (~isfield(options, "interpolate_interface"))
    options.interpolate_interface = false(size(bearing_surf));
  endif
  
  for i=1:numel(bearing_surf)
    nodes = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).nodes;
    elements = getfield(mesh.elements, options.elem_type)(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements, :);

    if (isfield(bearing_surf(i), "relative_tolerance") && isfield(bearing_surf(i), "absolute_tolerance"))
      z = bearing_surf(i).R(:, 3).' * (mesh.nodes(nodes, 1:3).' - bearing_surf(i).X0);

      dz = abs([0.5 * bearing_surf(i).w - max(z);
                0.5 * bearing_surf(i).w + min(z)]);

      tol = bearing_surf(i).relative_tolerance * bearing_surf(i).w + bearing_surf(i).absolute_tolerance;

      if (~options.interpolate_interface(i) && any(dz > tol))
        error("z-coordinate of nodes at bearing %d does not cover the complete bearing width", i);
      endif
    endif

    clear z dz tol elements nodes;
  endfor

  idx = int32(0);

  for i=1:numel(bearing_surf)
    if (isfield(bearing_surf(i).options, "mesh_size"))
      dx = bearing_surf(i).options.mesh_size;
      Nx = round(2 * pi * bearing_surf(i).r / dx);
      Nz = round(bearing_surf(i).w / dx);
    else
      Nx = bearing_surf(i).options.number_of_nodes_x;
      Nz = bearing_surf(i).options.number_of_nodes_z;
    endif

    N = [max(4, Nx), max(3, Nz)];

    if (isfield(bearing_surf(i).options, "bearing_model") && ischar(bearing_surf(i).options.bearing_model))
      switch (bearing_surf(i).options.bearing_model)
        case "EHD/FE"
          if (mod(N(1), 2) == 0)
            ++N(1);
          endif
          if (mod(N(2), 2) == 0)
            ++N(2);
          endif
      endswitch
    endif

    bearing_surf(i).grid_x = linspace(0, 2 * pi * bearing_surf(i).r, N(1));
    bearing_surf(i).grid_z = linspace(-0.5 * bearing_surf(i).w, 0.5 * bearing_surf(i).w, N(2));
    idx += (numel(bearing_surf(i).grid_x) - 1) * numel(bearing_surf(i).grid_z);
    clear dx N Nx Nz;
  endfor
endfunction
