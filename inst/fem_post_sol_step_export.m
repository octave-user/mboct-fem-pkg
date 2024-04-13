## Copyright (C) 2019(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_post_sol_step_export(@var{filename}, @var{sol}, @var{idx_sol}, @var{idx_t}, @var{t})
## @deftypefnx {} fem_post_sol_step_export(@dots{}, @var{scale})
## Export stresses and displacements of a single time step or mode shape to a file in Gmsh format.
##
## @var{filename} @dots{} Output filename
##
## @var{sol} @dots{} Finite element solution data structure.
##
## @var{idx_sol} @dots{} Physical index in the solution array to be exported. This is the array index for Octave.
##
## @var{idx_t} @dots{} Logical index of the time step to be exported. This is the index Gmsh will use.
##
## @var{t} @dots{} Time value related to @var{idx_t}. This is the value Gmsh will print to screen for this time step.
##
## @var{scale} @dots{} Multiply displacements by @var{scale} before exporting it. The actual display scale may be changed in Gmsh.
##
## @seealso{fem_post_sol_export, fem_post_sol_external}
## @end deftypefn

function fem_post_sol_step_export(filename, sol, idx_sol, idx_t, t, scale)
  if (nargin < 5 || nargin > 6 || nargout > 0)
    print_usage();
  endif

  if (nargin < 6)
    scale = 1;
  endif

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(filename, "w");

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    nodal_fields = struct("name", {"def", "vel", "theta", "p"}, ...
                          "title", {"nodal deformation", "nodal velocity", "nodal temperature", "acoustic pressure"}, ...
                          "dim", {3, 3, 1, 1}, ...
                          "value", {@(x) x(:, 1:3, idx_sol), @(x) x(:, 1:3, idx_sol), @(x) x(:, idx_sol), @(x) x(:, idx_sol)});

    for i=1:numel(nodal_fields)
      if (isfield(sol, nodal_fields(i).name))
        x = nodal_fields(i).value(getfield(sol, nodal_fields(i).name));
        if (isreal(x))
          fprintf(fd, "$NodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n%d\n%d\n", nodal_fields(i).title, t, idx_t - 1, nodal_fields(i).dim, size(x, 1));
          fprintf(fd, "%d %.16e\n", [1:size(x, 1); x.']);
          fputs(fd, "$EndNodeData\n");
        else
          fprintf(fd, "$NodeData\n1\n\"real(%s)\"\n1\n%e\n3\n%d\n%d\n%d\n", nodal_fields(i).title, t, idx_t - 1, nodal_fields(i).dim, size(x, 1));
          fprintf(fd, "%d %.16e\n", [1:size(x, 1); real(x.')]);
          fputs(fd, "$EndNodeData\n");
          fprintf(fd, "$NodeData\n1\n\"imag(%s)\"\n1\n%e\n3\n%d\n%d\n%d\n", nodal_fields(i).title, t, idx_t - 1, nodal_fields(i).dim, size(x, 1));
          fprintf(fd, "%d %.16e\n", [1:size(x, 1); imag(x.')]);
          fputs(fd, "$EndNodeData\n");
        endif
      endif
    endfor

    field_type = {"stress", "strain", "particle_velocity"};

    idxtens = int32([1, 4, 6, 4, 2, 5, 6, 5, 3]);
    stress_type = {"discontinuous stress tensor", "continuous stress tensor", ...
                   "van Mises stress", "discontinuous strain", "continuous strain", ...
                   "particle velocity vector", "particle velocity normal"};
    stress_field = {"tau", "taum", "vmis", "epsilon", "epsilonm", "v", "vn"};
    elem_stress = {"iso4", "quad8", "quad9", "iso8", "iso8f", "iso8upc", "iso20", "iso20upc", "iso27", "iso27upc", "tet10", "penta15", "tet10h", "tet10upc", "tet20", "iso20r", "iso20fr", "quad8r", "iso20upcr", "penta15upc"};
    stress_comp = int32([9, 9, 1, 9, 9, 3, 1]);

    for n=1:numel(field_type)
      if (isfield(sol, field_type{n}))
        for l=1:numel(stress_field)
          if (isfield(getfield(sol, field_type{n}), stress_field{l}))
            inumelem = int32(0);

            for i=1:numel(elem_stress)
              if (isfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i}))
                tau = getfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i});
                inumelem += rows(tau);
              endif
            endfor

            fprintf(fd, "$ElementNodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n%d\n%d\n", stress_type{l}, t, idx_t - 1, stress_comp(l), inumelem);

            inumelem = int32(0);

            for i=1:numel(elem_stress)
              if (isfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i}))
                switch (elem_stress{i})
                  case "iso3"
                    idxnode = int32([1:3]);
                  case "iso4"
                    idxnode = int32([1:4]);
                  case {"tria6", "tria6h"}
                    idxnode = int32([1:6]);
                  case "tria10"
                    idxnode = int32([1:10]);
                  case "quad8"
                    idxnode = int32([1:8]);
                  case "quad9"
                    idxnode = int32([1:9]);                    
                  case "quad8r"
                    idxnode = int32([3,4,1,2,7,8,5,6]);
                  case {"iso8", "iso8f", "iso8upc"}
                    idxnode = int32([5:8, 1:4]);
                  case {"iso20", "iso20upc"}
                    idxnode = int32([5:8, 1:4, 17, 19, 20, 18, 9, 12, 14, 10, 11, 13, 15, 16]);
                  case {"iso27", "iso27upc"}
                    idxnode = int32([1,2,3,4,5,6,7,8,9,12,14,10,11,13,15,16,17,19,20,18,21,22,24,25,23,26,27]);
                  case {"iso20r", "iso20fr", "iso20upcr"}
                    idxnode = int32([1,2,3,4,5,6,7,8,9,12,14,10,17,19,20,18,11,13,15,16]);
                  case {"penta15", "penta15upc"}
                    idxnode =  int32([1, 2, 3, 4, 5, 6, 7, 10, 8, 13, 15, 14, 9, 11, 12]);
                  case {"tet10", "tet10h", "tet10upc"}
                    idxnode = int32([1:8, 10, 9]);
                  case "tet20"
                    idxnode([1,5,6,2,7,8,3,9,10,17,16,20,14,19,12,18,15,13,11,4]) = int32(1:20);
                  otherwise
                    error("unknown element type \"%s\"", elem_stress{i});
                endswitch

                tau = getfield(getfield(getfield(sol, field_type{n}), stress_field{l}), elem_stress{i});

                tauout = zeros(2 + numel(idxnode) * stress_comp(l), rows(tau));
                tauout(1, :) = inumelem + (1:rows(tau));
                tauout(2, :) = numel(idxnode);

                for k=1:numel(idxnode)
                  switch (stress_comp(l))
                    case 1
                      taukl = tau(:, idxnode(k), idx_sol);
                    case 9
                      taukl = tau(:, idxnode(k), idxtens, idx_sol);
                    otherwise
                      taukl = tau(:, idxnode(k), :, idx_sol);
                  endswitch

                  tauout(2 + (k - 1) * stress_comp(l) + (1:stress_comp(l)), :) = reshape(taukl, rows(tau), stress_comp(l)).';
                endfor

                inumelem += rows(tau);

                format = ["%d %d", repmat(" %.16e", 1, numel(idxnode) * stress_comp(l)), "\n"];

                fprintf(fd, format, tauout);
              endif
            endfor

            fputs(fd, "$EndElementNodeData\n");
          endif
        endfor
      endif
    endfor
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

