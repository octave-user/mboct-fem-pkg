function comp_mat = fem_ehd_pre_comp_mat_load(filename)
  comp_mat.D = comp_mat.E = [];

  fd = -1;
  unwind_protect
    [fd, msg] = fopen(filename);

    if (fd == -1)
      error("failed to open file \"%s\"", filename)
    endif

    while (true)
      line = fgets(fd);

      if (~ischar(line))
        break;
      endif

      [d, count] = sscanf(line, "bearing diameter: %g", "C");

      if (count == 1)
        comp_mat.bearing_dimensions.bearing_diameter = d;
        continue;
      endif

      [w, count] = sscanf(line, "bearing width: %g", "C");

      if (count == 1)
        comp_mat.bearing_dimensions.bearing_width = w;
        continue;
      endif

      [n, count] = sscanf(line, "circumferential grid: %d", "C");

      if (count == 1)
        x = fscanf(fd, "%g", [n, 2]);
        comp_mat.bearing_surf.grid_x = x(:, 2);
        continue;
      endif

      [n, count] = sscanf(line, "axial grid: %d", "C");

      if (count == 1)
        z = fscanf(fd, "%g", [n, 2]);
        comp_mat.bearing_surf.grid_z = z(:, 2);
        continue;
      endif

      [rows, cols, count] = sscanf(line, "substruct total contrib matrix: %d %d", "C");

      if (count == 2)
        comp_mat.D = fscanf(fd, "%g", [rows, cols]);
        continue;
      endif

      [rows, cols, count] = sscanf(line, "substruct total residual matrix: %d %d", "C");

      if (count == 2)
        comp_mat.E = fscanf(fd, "%g", [rows, cols]);
        continue;
      endif
    endwhile
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
