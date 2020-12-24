function [bearing_surf, idx] = fem_ehd_pre_comp_mat_grid(mesh, bearing_surf)
  for i=1:numel(bearing_surf)
    nodes = mesh.groups.tria6(bearing_surf(i).group_idx).nodes;
    elements = mesh.elements.tria6(mesh.groups.tria6(bearing_surf(i).group_idx).elements, :);

    if (isfield(bearing_surf(i), "relative_tolerance") && isfield(bearing_surf(i), "absolute_tolerance"))
      z = bearing_surf(i).R(:, 3).' * (mesh.nodes(nodes, 1:3).' - bearing_surf(i).X0);

      dz = abs([0.5 * bearing_surf(i).w - max(z);
	        0.5 * bearing_surf(i).w + min(z)]);

      tol = bearing_surf(i).relative_tolerance * bearing_surf(i).w + bearing_surf(i).absolute_tolerance;
      
      if (any(dz > tol))
        error("z-coordinate of nodes at bearing %d does not cover the complete bearing width", i);
      endif
    endif

    clear z dz tol elements nodes;
  endfor
  
  idx = int32(0);
  
  for i=1:numel(bearing_surf)
    dx = bearing_surf(i).options.mesh_size;
    N = [max(4, round(2 * pi * bearing_surf(i).r / dx)), max(3, round(bearing_surf(i).w / dx))];

    if (isfield(bearing_surf(i).options, "bearing_model") && ischar(bearing_surf(i).options.bearing_model))
      switch (bearing_surf(i).options.bearing_model)
        case "EHD/FE"
          for j=1:numel(N)
            if (mod(N(j), 2) == 0)
              ++N(j);
            endif
          endfor
      endswitch
    endif
    
    bearing_surf(i).grid_x = linspace(0, 2 * pi * bearing_surf(i).r, N(1));
    bearing_surf(i).grid_z = linspace(-0.5 * bearing_surf(i).w, 0.5 * bearing_surf(i).w, N(2));
    idx += (numel(bearing_surf(i).grid_x) - 1) * numel(bearing_surf(i).grid_z);
    clear dx N;
  endfor
endfunction
