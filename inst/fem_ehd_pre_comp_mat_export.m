## Copyright (C) 2016(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_ehd_pre_comp_mat_export(@var{comp_mat}, @var{options}, @var{output_file})
## @deftypefnx {} fem_ehd_pre_comp_mat_export(@var{comp_mat}, @var{options}, @var{output_filename})
## Print a compliance matrix to a file which can be loaded by MBDyn's module-hydrodynamic_plain_bearing2.
##
## @var{comp_mat} @dots{} Return value from fem_ehd_pre_comp_mat_struct or fem_ehd_pre_comp_mat_unstruct.
##
## @var{options}.matrix_type @dots{} One of "modal", "nodal", "nodal substruct", "nodal substruct total".
##
## @var{output_file} @dots{} Open file descriptor where the compliance matrices are written to.
##
## @var{output_filename} @dots{} Output filename where the compliance matrices are written to.
##
## @seealso{fem_ehd_pre_comp_mat_struct, fem_ehd_pre_comp_mat_unstruct}
## @end deftypefn

function fem_ehd_pre_comp_mat_export(comp_mat, options, varargin)
  if (nargin < 1 || nargin > 3 || nargout > 0)
    print_usage();
  endif

  if (~(isstruct(comp_mat) && isscalar(comp_mat)))
    error("argument comp_mat must be a scalar struct");
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "matrix_type"))
    options.matrix_type = "default";
  endif

  switch (options.matrix_type)
    case "default"
      if (isfield(comp_mat, "C"))
        options.matrix_type = "nodal";
      else
        options.matrix_type = "modal";
      endif
  endswitch

  if (~isfield(options, "modal_threshold"))
    options.modal_threshold = 1e-12;
  endif

  owns_fd = false;
  fd = -1;

  unwind_protect
    if (nargin < 3)
      fd = stdout;
    elseif (ischar(varargin{1}))
      owns_fd = true;
      [fd, msg] = fopen(varargin{1}, "w");
      if (fd == -1)
        error("failed to open file \"%s\": %s", varargin{1}, msg);
      endif
    else
      fd = varargin{1};
    endif

    file_format = 12;

    d = comp_mat.bearing_dimensions.bearing_diameter;
    w = comp_mat.bearing_dimensions.bearing_width;

    x = comp_mat.bearing_surf.grid_x;
    z = comp_mat.bearing_surf.grid_z;
    z -= (max(z) + min(z)) / 2;

    N0 = numel(x) * numel(z);
    N1 = (numel(x) - 1) * numel(z);

    if (numel(x) < 3)
      error("numel(comp_mat.bearing_surf.grid_x must be at least three")
    endif

    if (abs((x(end) - x(1)) / (d * pi) - 1) >  sqrt(eps))
      error("invalid range for comp_mat.bearing_surf.grid_x");
    endif

    ic = zeros(numel(x), numel(z), "int32");

    for i=1:numel(x)
      for j=1:numel(z)
        ic(i, j) = (i - 1) * length(z) + j;
      endfor
    endfor

    ix = [(1:numel(x) + 1) + 2, 1:2];
    iz = 1:numel(z);
    x = [x([end - (2:-1:1)]) - d * pi, x(1:end - 1), x(1:2) + d * pi];
    ic = ic([1:end - 1, 1:2, end - 2:end - 1], :);

    fprintf(fd, "file format:\t%d\n", file_format);
    fprintf(fd, "\nbearing diameter:\t%.16g\n", d);
    fprintf(fd, "\nbearing width:\t%.16g\n", w);
    fprintf(fd, "\ncircumferential grid:\t%d\n", numel(ix));

    for i=1:numel(x)
      fprintf(fd, "%d\t%.16g\n", i, x(i));
    endfor

    fprintf(fd, "\naxial grid:\t%d\n", numel(z));

    for i=1:numel(z)
      fprintf(fd, "%d\t%.16g\n", i, z(i));
    endfor

    fprintf(fd, "\nnodes:\t%d\t%d\n", numel(ic), N1);
    idx = int32(0);

    for i=1:rows(ic)
      for j=1:columns(ic)
        fprintf(fd, "%d\t%d\t%d\t%d\n", ++idx, ic(i, j), ix(i), iz(j));
      endfor
    endfor

    fprintf(fd, "\nreference pressure:\t%e\n", comp_mat.reference_pressure);

    warning("error", "Octave:singular-matrix", "local");
    warning("error", "Octave:nearly-singular-matrix", "local");

    fprintf(fd, "\nbearing center:\t");
    fprintf(fd, "%.16g\t", comp_mat.dX);
    fprintf(fd, "\n");

    fprintf(fd, "\nbearing orientation:\t");

    for i=1:3
      fprintf(fd, "%.16g\t", comp_mat.dR(i, :));
      fprintf(fd, "\n\t\t\t");
    endfor

    switch (options.matrix_type)
      case "nodal"
        if (~isfield(comp_mat, "C"))
          error("missing field \"C\"");
        endif

        if (rows(comp_mat.C) ~= N0 || columns(comp_mat.C) ~= N0)
          error("invalid size for comp_mat.C");
        endif

        C = comp_mat.C(1:end - numel(z), 1:end - numel(z));

        print_matrix(fd, C, "compliance matrix");

        fprintf(fd, "\n\n## cond(C)=%e\n", cond(C));
      case "nodal substruct"
        if (~isfield(comp_mat, "C"))
          error("missing field \"C\"");
        endif

        if (~isfield(comp_mat, "D"))
          error("missing field \"D\"");
        endif

        if (~isfield(comp_mat, "E"))
          error("missing field \"E\"");
        endif

        if (rows(comp_mat.C) ~= N0 || columns(comp_mat.C) ~= N0)
          error("invalid size for comp_mat.C");
        endif

        if (rows(comp_mat.D) ~= N0 )
          error("invalid size for comp_mat.D");
        endif

        if (columns(comp_mat.E) ~= N0)
          error("invalid size for comp_mat.E");
        endif

        C = comp_mat.C(1:end - numel(z), 1:end - numel(z));

        [mode_index, D, E] = extract_modes(comp_mat.D, comp_mat.E, z, options.modal_threshold);

        print_index_vector(fd, mode_index, "modal subset vector");

        print_matrix(fd, C, "compliance matrix substruct");

        fprintf(fd, "\n\n## cond(C)=%e\n", cond(C));

        print_matrix(fd, D, "substruct contrib matrix");
        print_matrix(fd, E, "substruct residual matrix");
      case {"nodal substruct total", "modal substruct total"}
        if (~isfield(comp_mat, "D"))
          error("missing field \"D\"");
        endif

        if (~isfield(comp_mat, "E"))
          error("missing field \"E\"");
        endif

        if (rows(comp_mat.D) ~= N0 )
          error("invalid size for comp_mat.D");
        endif

        if (columns(comp_mat.E) ~= N0)
          error("invalid size for comp_mat.E");
        endif

        [mode_index, D, E] = extract_modes(comp_mat.D, comp_mat.E, z, options.modal_threshold);

        print_index_vector(fd, mode_index, "modal subset vector");
        print_matrix(fd, D, "substruct total contrib matrix");
        print_matrix(fd, E, "substruct total residual matrix");
      case "modal"
        if (~isfield(comp_mat, "KPhi"))
          error("missing field \"KPhi\"");
        endif

        if (~isfield(comp_mat, "RPhi"))
          error("missing field \"RPhi\"");
        endif

        if (~isfield(comp_mat, "Phin"))
          error("missing field \"Phin\"");
        endif

        if (rows(comp_mat.Phin) ~= N0)
          error("invalid size for comp_mat.Phin");
        endif

        if (columns(comp_mat.RPhi) ~= N0)
          error("invalid size for comp_mat.RPhi");
        endif

        Phin = comp_mat.Phin(1:end - numel(z), :);
        KPhi = comp_mat.KPhi;
        RPhi = comp_mat.RPhi(:, 1:end - numel(z));

        print_matrix(fd, Phin, "mode shapes");
        print_matrix(fd, KPhi \ RPhi, "modal load");

        fprintf(fd, "\n\n## cond(KPhi)=%e\n", cond(KPhi));
      otherwise
        error("unknown option matrix_type=\"%s\"", options.matrix_type);
    endswitch

  unwind_protect_cleanup
    if (owns_fd && fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

function print_matrix(fd, C, C_name)
  fprintf(fd, "\n%s:\t%d\t%d\n", C_name, rows(C), columns(C));

  for i=1:rows(C)
    fprintf(fd, "%.16e\t", C(i, :));
    fprintf(fd, "\n");
  endfor
endfunction

function print_index_vector(fd, mode_index, name)
  fprintf(fd, "\n%s:\t%d\n", name, numel(mode_index));
  fprintf(fd, "%d\n", mode_index);
endfunction

function [mode_index, D, E] = extract_modes(D, E, z, threshold)
  mode_index = find((norm(D, "cols")(:) / norm(D) >= threshold) | (norm(E, "rows") / norm(E) >= threshold));
  D = D(1:end - numel(z), mode_index);
  E = E(mode_index, 1:end - numel(z));
endfunction
