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
## "modal node" @dots{} Deformations of the flexible body are scaled with respect to MBDyn's modal node.
## "least square" @dots{} Subtract the rigid body motion from nodal displacements before scaling.
## "reference node" @dots{} Deformations of the flexible body are scaled with respect to an arbitrary point defined by sol_dyn.bodies(i).X_ref and sol_dyn.bodies(i).R_ref.
## With the options above, rigid body motion of the modal node or reference node will be added.
## All scaling options ending with an asterisk do not add the rigid body motion to nodal displacements.
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

%!test
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!             0,  0.5 * b,  0.5 * c;  ##  2
%!             0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!             0,  0.5 * b, -0.5 * c;  ##  6
%!             0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!             a,  0.5 * b,  0.5 * c;  ##  9
%!             a, -0.5 * b,  0.5 * c;  ## 10
%!             a,  0.5 * b, -0.5 * c;  ## 11
%!             a, -0.5 * b, -0.5 * c,  ## 12
%!         a + d,        0,        0;  ## 13
%!            -e,        0,        0]; ## 14
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! mesh.material_data.E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! cms_data.load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.number_of_threads = 1;
%! cms_opt.algorithm = "eliminate";
%! cms_opt.invariants = true;
%! [cms_data.mesh, ...
%!  cms_data.mat_ass, ...
%!  cms_data.dof_map, ...
%!  cms_data.sol_eig, ...
%!  cms_data.cms_opt] = fem_cms_create(mesh, cms_data.load_case, cms_opt);
%! sol_dyn.t = linspace(0, 1e-3, 100);
%! sol_dyn.bodies.q = zeros(columns(cms_data.mat_ass.Tred), numel(sol_dyn.t));
%! f = 1000 + rand(1, columns(cms_data.mat_ass.Tred)) * 1000;
%! for i=1:rows(sol_dyn.bodies.q)
%!   sol_dyn.bodies.q(i, :) = sin(2 * pi * f(i) * sol_dyn.t);
%! endfor
%! sol_dyn.bodies.X = zeros(3, numel(sol_dyn.t));
%! sol_dyn.bodies.R = zeros(3, 3, numel(sol_dyn.t));
%! for i=1:numel(sol_dyn.t)
%!   sol_dyn.bodies.R(:, :, i) = eye(3);
%! endfor
%! sol_dyn.bodies.X_ref = zeros(3, numel(sol_dyn.t));
%! sol_dyn.bodies.R_ref = sol_dyn.bodies.R;
%! scale_types = {"modal node", "modal node*", "least square", "least square*", "reference node", "reference node*"};
%! sol_tot = struct("bodies", cell(numel(scale_types), 2));
%! for j=1:2
%!   for i=1:numel(scale_types)
%!     opt_exp.scale = 1;
%!     opt_exp.scale_type = scale_types{i};
%!     opt_exp.output_stress = FEM_SCA_STRESS_VMIS;
%!     sol_tot(i, j) = fem_post_cms_expand(sol_dyn, cms_data, 1:10:numel(sol_dyn.t), opt_exp);
%!   endfor
%!   if (j == 1)
%!     cms_data.cms_opt.selected_modes = 1:2:columns(cms_data.mat_ass.Tred);
%!     sol_dyn.bodies.q = sol_dyn.bodies.q(cms_data.cms_opt.selected_modes, :);
%!   endif
%! endfor

%!demo
%! close all;
%! SI_unit_m = 1e-3;
%! SI_unit_kg = 1e3;
%! SI_unit_s = 1e-1;
%! SI_unit_N = SI_unit_kg * SI_unit_m / SI_unit_s^2;
%! SI_unit_Pa = SI_unit_N / SI_unit_m^2;
%! a = 150e-3 / SI_unit_m;
%! b = 20e-3 / SI_unit_m;
%! c = 45e-3 / SI_unit_m;
%! d = 10e-3 / SI_unit_m;
%! e = 10e-3 / SI_unit_m;
%! X = [ 0.5 * a,  0.5 * b,  0.5 * c;  ##  1
%!             0,  0.5 * b,  0.5 * c;  ##  2
%!             0, -0.5 * b,  0.5 * c;  ##  3
%!       0.5 * a, -0.5 * b,  0.5 * c;  ##  4
%!       0.5 * a,  0.5 * b, -0.5 * c;  ##  5
%!             0,  0.5 * b, -0.5 * c;  ##  6
%!             0, -0.5 * b, -0.5 * c;  ##  7
%!       0.5 * a, -0.5 * b, -0.5 * c,  ##  8
%!             a,  0.5 * b,  0.5 * c;  ##  9
%!             a, -0.5 * b,  0.5 * c;  ## 10
%!             a,  0.5 * b, -0.5 * c;  ## 11
%!             a, -0.5 * b, -0.5 * c,  ## 12
%!         a + d,        0,        0;  ## 13
%!            -e,        0,        0]; ## 14
%! mesh.nodes = [X, zeros(rows(X), 3)];
%! mesh.elements.iso8 = int32([1:8; 9, 1, 4, 10, 11, 5, 8, 12]);
%! mesh.materials.iso8 = int32([1; 1]);
%! mesh.elements.rbe3(1).nodes = int32([13, 9, 10, 11, 12]);
%! mesh.elements.rbe3(1).weight = ones(1, 4);
%! mesh.elements.rbe3(2).nodes = int32([14, 2, 3, 6, 7]);
%! mesh.elements.rbe3(2).weight = ones(1, 4);
%! mesh.material_data.E = 210000e6 / (SI_unit_N / SI_unit_m^2);
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850 / (SI_unit_kg / SI_unit_m^3);
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! cms_data.load_case.locked_dof = false(rows(mesh.nodes), 6);
%! cms_opt.verbose = false;
%! cms_opt.modes.number = int32(6);
%! cms_opt.nodes.modal.number = int32(14);
%! cms_opt.nodes.interfaces.number = int32(13);
%! cms_opt.number_of_threads = 1;
%! cms_opt.algorithm = "eliminate";
%! cms_opt.invariants = true;
%! [cms_data.mesh, ...
%!  cms_data.mat_ass, ...
%!  cms_data.dof_map, ...
%!  cms_data.sol_eig, ...
%!  cms_data.cms_opt] = fem_cms_create(mesh, cms_data.load_case, cms_opt);
%! sol_dyn.t = linspace(0, 1e-3, 100);
%! sol_dyn.bodies.q = zeros(columns(cms_data.mat_ass.Tred), numel(sol_dyn.t));
%! f = 1000 + rand(1, columns(cms_data.mat_ass.Tred)) * 1000;
%! for i=1:rows(sol_dyn.bodies.q)
%!   sol_dyn.bodies.q(i, :) = sin(2 * pi * f(i) * sol_dyn.t);
%! endfor
%! sol_dyn.bodies.X = zeros(3, numel(sol_dyn.t));
%! sol_dyn.bodies.R = zeros(3, 3, numel(sol_dyn.t));
%! for i=1:numel(sol_dyn.t)
%!   sol_dyn.bodies.R(:, :, i) = eye(3);
%! endfor
%! sol_dyn.bodies.X_ref = zeros(3, numel(sol_dyn.t));
%! sol_dyn.bodies.R_ref = sol_dyn.bodies.R;
%! opt_exp.scale = 1;
%! opt_exp.scale_type = "least square";
%! opt_exp.output_stress = FEM_SCA_STRESS_VMIS;
%! sol_tot = fem_post_cms_expand(sol_dyn, cms_data, 1:10:numel(sol_dyn.t), opt_exp);
