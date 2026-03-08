function [D, E] = fem_ehd_pre_comp_mat_load(filename)
  D = E = [];
  
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

      [rows, cols, count] = sscanf(line, "substruct total contrib matrix: %d %d", "C");

      if (count == 2)
        D = fscanf(fd, "%g", [rows, cols]);
        continue;
      endif

      [rows, cols, count] = sscanf(line, "substruct total residual matrix: %d %d", "C");

      if (count == 2)
        E = fscanf(fd, "%g", [rows, cols]);
        continue;
      endif
    endwhile
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
