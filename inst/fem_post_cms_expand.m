## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{sol_tot} = fem_post_cms_expand(@var{sol_dyn}, @var{cms_data}, @var{idx_t}, @var{options})
##
## Expand a modal solution to a nodal solution, compute stresses and scale deformations
##
## @var{sol_dyn} @dots{} modal solution returned from fem_post_cms_sol_import
##
## @var{cms_data} @dots{} struct array containing mesh, dof_map, and mat_ass fields of several flexible bodies
##
## @var{idx_t} @dots{} index of time steps which should be converted
##
## @var{options}.scale @dots{} scale factor for deformations
##
## @var{options}.scale_type @dots{} options for scaling
##
## Valid values are "modal node", "modal node*", "least square", "least square*", "reference node", "reference node*"
## @seealso{fem_post_cms_expand_body, fem_post_cms_sol_import}
## @end deftypefn

function sol_tot = fem_post_cms_expand(sol_dyn, cms_data, idx_t, options)
  if (nargin ~= 4 || nargout > 1)
    print_usage();
  endif

  sol_tot.bodies = struct("def", cell(1, numel(sol_dyn)));

  if (~isfield(options, "scale"))
    options.scale = 1;
  endif
  
  if (~isfield(options, "scale_type"))
    options.scale_type = "modal node";
  endif

  if (~isfield(options, "output_stress"))
    options.output_stress = int32(-1);
  endif

  for i=1:numel(sol_dyn.bodies)
    if (isfield(cms_data(i).cms_opt, "selected_modes"))
      q = zeros(columns(cms_data(i).mat_ass.Tred), numel(idx_t));
      q(cms_data(i).cms_opt.selected_modes, :) = sol_dyn.bodies(i).q(:, idx_t);
    else
      q = sol_dyn.bodies(i).q(:, idx_t);
    endif
    
    Xn = cms_data(i).mesh.nodes;
    Xm = Xn(cms_data(i).cms_opt.nodes.modal.number, :);
    Un = fem_post_cms_expand_body(cms_data(i).mesh, ...
                                  cms_data(i).dof_map, ...
                                  cms_data(i).mat_ass, ...
                                  q);

    if (options.output_stress ~= -1)
      sol_tot.bodies(i).def = Un;
      sol_tot.bodies(i).stress = fem_ass_matrix(cms_data(i).mesh, ...
                                                cms_data(i).dof_map, ...
                                                options.output_stress, ...
                                                cms_data(i).load_case, ...
                                                sol_tot.bodies(i));
    endif

    switch (options.scale_type)
      case {"modal node", "modal node*"}
        Un *= options.scale;
      case {"least square", "least square*"}
        A = rigid_body_trans_mat(Xn);

        Un3 = zeros(3 * rows(Un), size(Un, 3));

        for j=1:3
          Un3(j:3:end, :) = reshape(Un(:, j, :), rows(Un), size(Un, 3));
        endfor

        U0 = (A.' * A) \ (A.' * Un3);
        Urb3 = A * U0;

        switch (options.scale_type)
          case "least square"
            Un3 = options.scale * (Un3 - Urb3) + Urb3;
          case "least square*"
            Un3 -= Urb3;
        endswitch
        
        for j=1:3
          Un(:, j, :) = Un3(j:3:end, :);
        endfor
      case {"reference node", "reference node*"}
        A = rigid_body_trans_mat(Xn);
      otherwise
        error("unknown value for options.scale_type=%s", options.scale_type);
    endswitch

    switch (options.scale_type)
      case {"modal node*", "least square*", "reference node*"}
        sol_tot.bodies(i).def = Un;
      otherwise
        ln = Un + (Xn - Xm);

        sol_tot.bodies(i).def = zeros(size(Xn, 1), 6, numel(idx_t));

        for j=1:3
          sol_tot.bodies(i).def(:, j, :) = -Xn(:, j) + sol_dyn.bodies(i).X(j, idx_t);

          for k=1:3
            sol_tot.bodies(i).def(:, j, :) += ln(:, k, :) .* sol_dyn.bodies(i).R(j, k, idx_t);
          endfor
        endfor
    endswitch
    
    switch (options.scale_type)
      case "reference node"
        def_rb = zeros(size(Xn, 1), 6, numel(idx_t));
        ln_ref = Xn - Xm;
        
        for j=1:3
          def_rb(:, j, :) = -Xn(:, j) + sol_dyn.bodies(i).X_ref(j, idx_t);

          for k=1:3
            def_rb(:, j, :) += ln_ref(:, k, :) .* sol_dyn.bodies(i).R_ref(j, k, idx_t);
          endfor
        endfor

        sol_tot.bodies(i).def = (sol_tot.bodies(i).def - def_rb) * options.scale + def_rb;
    endswitch
  endfor
endfunction

function A = rigid_body_trans_mat(Xn)
  A = zeros(3 * rows(Xn), 6);

  ## A(i:i+2, :) = [eye(3), -skew(Xn(i, :))]
  for j=1:3
    A(j:3:end, j) = 1;
  endfor

  A(1:3:end, 5) = Xn(:, 3);
  A(2:3:end, 4) = -Xn(:, 3);
  A(1:3:end, 6) = -Xn(:, 2);
  A(3:3:end, 4) = Xn(:, 2);
  A(2:3:end, 6) = Xn(:, 1);
  A(3:3:end, 5) = -Xn(:, 1);
endfunction
