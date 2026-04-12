## Copyright (C) 2026(-2026) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{comp_mat}] = fem_ehd_pre_comp_load(@var{filename})
## Load a compliance matrix from an input file, suitable for MBDyn's module-hydrodynamic_plain_bearing2
##
## @var{filename} @dots{} Path of the file to be loaded
##
## @var{comp_mat} @dots{} Compliance matrix data
##
## @seealso{fem_ehd_pre_comp_mat_export}
## @end deftypefn

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
