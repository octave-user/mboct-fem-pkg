## Copyright (C) 2019(-2020) Reinhard <octave-user@a1.net>
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
## @var{idx_sol} @dots{} Data index in the solution array to be exported
##
## @var{idx_t} @dots{} Index of the time step to be exported
##
## @var{t} @dots{} Time value related to @var{idx_t}
##
## @var{scale} @dots{} Multiply displacements by @var{scale} before exporting it
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

  [fd, msg] = fopen(filename, "wt");

  if (fd == -1)
    error("failed to open file \"%s\"", filename);
  endif

  unwind_protect
    fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");
    fprintf(fd, "$NodeData\n1\n\"nodal deformation\"\n1\n%e\n3\n%d\n3\n%d\n", t, idx_t - 1, size(sol.def, 1));
    fprintf(fd, "%d %e %e %e\n", [1:size(sol.def, 1); scale * sol.def(:, 1:3, idx_sol).']);
    fputs(fd, "$EndNodeData\n");

    if (isfield(sol, "stress"))
      idxtens = int32([1, 4, 6, 4, 2, 5, 6, 5, 3]);

      stress_type = {"discontinuous stress tensor", "continuous stress tensor", "van Mises stress"};
      stress_field = {"tau", "taum", "vmis"};      
      stress_comp = int32([9, 9, 1]);
      
      for l=1:numel(stress_field)
        if (isfield(sol.stress, stress_field{l}))
          elem_stress = {"iso4", "iso8", "tet10"};

          inumelem = int32(0);
          
          for i=1:numel(elem_stress)
            if (isfield(getfield(sol.stress, stress_field{l}), elem_stress{i}))
              tau = getfield(getfield(sol.stress, stress_field{l}), elem_stress{i});
              inumelem += rows(tau);
            endif
          endfor

          fprintf(fd, "$ElementNodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n%d\n%d\n", stress_type{l}, t, idx_t - 1, stress_comp(l), inumelem);

          inumelem = int32(0);

          for i=1:numel(elem_stress)
            if (isfield(getfield(sol.stress, stress_field{l}), elem_stress{i}))
              switch (elem_stress{i})
                case "iso4"
                  idxnode = int32([1:4]);
                case "iso8"
                  idxnode = int32([5:8, 1:4]);
                case "tet10"
                  idxnode = int32([1:8, 10, 9]);
              endswitch
              
              tau = getfield(getfield(sol.stress, stress_field{l}), elem_stress{i});

              tauout = zeros(2 + numel(idxnode) * stress_comp(l), rows(tau));
              tauout(1, :) = inumelem + (1:rows(tau));
              tauout(2, :) = numel(idxnode);
              
              for k=1:numel(idxnode)
                if (stress_comp(l) == 1)
                  taukl = tau(:, idxnode(k), idx_sol);
                else
                  taukl = tau(:, idxnode(k), idxtens, idx_sol);
                endif
                
                tauout(2 + (k - 1) * stress_comp(l) + (1:stress_comp(l)), :) = reshape(taukl, rows(tau), stress_comp(l)).';
              endfor

              inumelem += rows(tau);
              
              format = ["%d %d", repmat(" %e", 1, numel(idxnode) * stress_comp(l)), "\n"];

              fprintf(fd, format, tauout);
            endif
          endfor

          fputs(fd, "$EndElementNodeData\n");
        endif
      endfor
    endif
  unwind_protect_cleanup
    fclose(fd);
  end_unwind_protect
endfunction
