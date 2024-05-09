## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{sol} = fem_post_cms_sol_merge(@var{mesh}, @var{dof_map}, @var{sol_tot})
##
## Merge the solutions (displacements, stresses) of several modal bodies
##
## @var{mesh} @dots{} the combined mesh returned from fem_post_mesh_merge
##
## @var{dof_map} @dots{} the combined degree of freedom mapping returned from fem_post_mesh_merge
##
## @var{sol_tot} @dots{} nodal solutions of several flexible bodies returned from fem_post_cms_expand
## @seealso{fem_post_cms_sol_import}
## @end deftypefn

function [sol] = fem_post_cms_sol_merge(mesh, dof_map, sol_tot)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  num_nodes = int32(0);
  num_steps = int32(-1);
  num_dof = int32(-1);

  for i=1:numel(sol_tot.bodies)
    num_nodes += size(sol_tot.bodies(i).def, 1);

    if (num_dof == -1)
      num_dof = size(sol_tot.bodies(i).def, 2);
    elseif (num_dof ~= size(sol_tot.bodies(i).def, 2))
      error("invalid number of dof detected");
    endif

    if (num_steps == -1)
      num_steps = size(sol_tot.bodies(i).def, 3);
    elseif (num_steps ~= size(sol_tot.bodies(i).def, 3))
      error("invalid number of load steps detected");
    endif
  endfor

  sol.def = zeros(num_nodes, num_dof, num_steps);

  for i=1:numel(sol_tot.bodies)
    nnodes = size(sol_tot.bodies(i).def, 1);
    sol.def(dof_map.submesh.offset.nodes(i) + (1:nnodes), :, :) = sol_tot.bodies(i).def;
  endfor

  if (isfield(sol_tot.bodies(i), "stress"))
    elem_type = {fem_pre_mesh_elem_type_dim([2, 3]).name};
    stress_type = {"tau", "taum", "vmis"};

    sol.stress = struct();

    for k=1:numel(stress_type)
      for j=1:numel(elem_type)
        num_elem = int32(0);
        num_nodes = int32(-1);
        num_comp = int32(-1);
        num_steps = int32(-1);

        stress_elem_type_m = struct();

        for i=1:numel(sol_tot.bodies)
          if (isfield(sol_tot.bodies(i).stress, stress_type{k}))
            stress_type_k = getfield(sol_tot.bodies(i).stress, stress_type{k});

            if (isfield(stress_type_k, elem_type{j}))
              stress_elem = getfield(stress_type_k, elem_type{j});
              num_elem += size(stress_elem, 1);

              if (num_nodes == -1)
                num_nodes = size(stress_elem, 2);
              elseif (num_nodes ~= size(stress_elem, 2))
                error("invalid number of element nodes detected");
              endif

              if (num_comp == -1)
                num_comp = size(stress_elem, 3);
              elseif (num_comp ~= size(stress_elem, 3))
                error("invalid number of stress components detected");
              endif

              if (num_steps == -1)
                num_steps = size(stress_elem, 4);
              elseif (num_steps ~= size(stress_elem, 4))
                error("invalid number of load steps detected");
              endif
            endif
          endif
        endfor

        if (num_elem > 0)
          stress_elem_m = zeros(num_elem, num_nodes, num_comp, num_steps);
          offset_elem = getfield(dof_map.submesh.offset.elements, elem_type{j});

          for i=1:numel(sol_tot.bodies)
            if (isfield(sol_tot.bodies(i).stress, stress_type{k}))
              stress_type_b = getfield(sol_tot.bodies(i).stress, stress_type{k});

              if (isfield(stress_type_b, elem_type{j}))
                stress_elem_b = getfield(stress_type_b, elem_type{j});
                stress_elem_m(offset_elem(i) + (1:size(stress_elem_b, 1)), :, :, :) = stress_elem_b;
              endif
            endif
          endfor

          stress_elem_type_m = setfield(stress_elem_type_m, elem_type{j}, stress_elem_m);
        endif

        if (numel(fieldnames(stress_elem_type_m)))
          sol.stress = setfield(sol.stress, stress_type{k}, stress_elem_type_m);
        endif
      endfor
    endfor
  endif
endfunction

