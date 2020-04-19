## Copyright (C) 2016(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{comp_mat} = fem_ehd_pre_comp_mat_struct(@var{bearing_dimensions}, @var{material}, @var{options})
## Compute a compliance matrix for use with MBDyn's module-hydrodynamic_plain_bearing2 for a simple structured mesh
##
## @var{bearing_dimensions}.bearing_diameter @dots{} Journal plain bearing diameter.
##
## @var{bearing_dimensions}.bearing_width @dots{} Journal plain bearing width.
##
## @var{bearing_dimensions}.bearing_geometry @dots{} Callback function for building mesh and boundary conditions.
##
## @var{material} @dots{} Material data for the Finite Element Mesh.
##
## @var{options}.verbose @dots{} Enable verbose output.
##
## @var{options}.output_matrices @dots{} Return all Finite Element matrices in @var{comp_mat}.
##
## @var{options}.output_solution @dots{} Output the Finite Element solution in @var{comp_mat}.
##
## @var{options}.number_of_modes @dots{} Number of mode shapes used for a mode based EHD method.
##
## @var{options}.reference_pressure @dots{} Nodal pressure values are scaled by @var{options}.reference_pressure.
##
## @seealso{fem_ehd_pre_comp_mat_export, fem_ehd_pre_comp_mat_plot}
## @end deftypefn

function [comp_mat] = fem_ehd_pre_comp_mat_struct(bearing_dimensions, material, options)
  if (nargin ~= 3 || nargout > 1)
    print_usage();
  endif

  if (~isfield(bearing_dimensions, "bearing_geometry"))
    error("missing field bearing_dimensions.bearing_geometry");
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  if (~isfield(options, "output_matrices"))
    options.output_matrices = false;
  endif

  if (~isfield(options, "output_solution"))
    options.output_solution = false;
  endif

  if (~isfield(options, "number_of_modes"))
    options.number_of_modes = -1;
  endif

  if (options.verbose)
    fprintf(stderr, "creating 2D mesh ...\n");
    tic();
  endif

  [geometry, loads, comp_mat.bearing_surf] = feval(bearing_dimensions.bearing_geometry, bearing_dimensions, options);

  if (options.verbose)
    toc();
    fprintf(stderr, "creating 3D mesh ...\n");
    tic();
  endif

  options_mesh.verbose = options.verbose;
  options_mesh.identical_boundary_cond = true;

  [mesh, load_cases] = fem_pre_mesh_struct_create(geometry, loads, material, options_mesh);

  if (options.verbose)
    toc();
    fprintf(stderr, "assembling matrices ...\n");
    tic();
  endif

  dof_map = fem_ass_dof_map(mesh, load_cases(1));
  
  [mat_ass.K, ...
   mat_ass.R, ...
   mat_ass.Rlump] = fem_ass_matrix(mesh, ...
                                   dof_map, ...
                                   [FEM_MAT_STIFFNESS_SYM, ...
                                    FEM_VEC_LOAD_CONSISTENT, ...
                                    FEM_VEC_LOAD_LUMPED], ...
                                   load_cases);

  if (options.verbose)
    toc();
    fprintf(stderr, "building compliance matrix ...\n");
    tic();
  endif

  comp_mat.bearing_dimensions = bearing_dimensions;

  if (options.output_matrices)
    comp_mat.load_cases = load_cases;
    comp_mat.mat_ass = mat_ass;
    comp_mat.loads = loads;
  endif

  comp_mat.reference_pressure = options.reference_pressure;
  comp_mat.mesh = mesh;
  comp_mat.dof_map = dof_map;
  N = length(comp_mat.bearing_surf.s) * length(comp_mat.bearing_surf.t);

  if (options.number_of_modes <= 0)
    sol_stat = fem_sol_static(mesh, dof_map, mat_ass, options);

    if (options.output_solution)
      comp_mat.sol_stat = sol_stat;
    endif

    comp_mat.C = zeros(N, N);

    for i=1:length(comp_mat.bearing_surf.s)
      for j=1:length(comp_mat.bearing_surf.t)
        idxrow = (i - 1) * length(comp_mat.bearing_surf.t) + j;
        inode = mesh.structured.inode_idx(1, i, j + comp_mat.bearing_surf.offset_t);
        for k=1:length(loads)
          idxcol = (loads(k).index.i - 1) * length(comp_mat.bearing_surf.t) + loads(k).index.j;
          U = sol_stat.def(inode, 1:2, k);
          n = mesh.nodes(inode, 1:2).';
          comp_mat.C(idxrow, idxcol) = U * (n / norm(n));
        endfor
      endfor
    endfor
  else
    if (options.verbose)
      fprintf(stderr, "building condensed matrices ...\n");
      tic();
    endif

    master_dof = zeros(1, N * 3, "int32");
    k = 0;

    for i=1:length(comp_mat.bearing_surf.s)
      for j=1:length(comp_mat.bearing_surf.t)
        inode = mesh.structured.inode_idx(1, i, j + comp_mat.bearing_surf.offset_t);
        node_dofs = dof_map.ndof(inode, 1:3);
        master_dof(k + (1:length(node_dofs))) = node_dofs;
        k += length(node_dofs);
      endfor
    endfor

    master_dof = unique(master_dof);

    [Kred, Tred] = fem_cms_red_guyan(mat_ass.K, master_dof);
    Tred = full(Tred);
    Kred = full(Kred);
    Rred = full(Tred.' * mat_ass.R);

    An = zeros(3, rows(mesh.nodes));

    for k=1:(length(load_cases) - length(comp_mat.bearing_surf.t)) ## last t load cases are redundant
      if (isfield(load_cases(k).pressure, "iso4"))
        for l=1:rows(load_cases(k).pressure.iso4.elements)
          for m=1:columns(load_cases(k).pressure.iso4.elements)
            p = load_cases(k).pressure.iso4.p(l, m);
            if (p ~= 0)
              inode = load_cases(k).pressure.iso4.elements(l, m);
              node_dof = dof_map.ndof(inode, 1:3);
              for i=1:length(node_dof)
                An(i, inode) += mat_ass.Rlump(node_dof(i), k) / p;
              endfor
            endif
          endfor
        endfor
      endif
    endfor

    A = zeros(1, length(mat_ass.K));

    for i=1:rows(mesh.nodes)
      node_dof = dof_map.ndof(i, 1:3);
      for j=1:length(node_dof)
        if (node_dof(j) > 0)
          A(node_dof(j)) = norm(An(:, i)); ## Why not using An(j, i) ??
        endif
      endfor
    endfor

    comp_mat.An = An;
    comp_mat.A = A;

    A = diag(A);
    Ared = Tred.' * A * Tred;

    if (options.verbose)
      toc();
      fprintf(stderr, "computing invariants of condensed stiffness matrix ...\n");
      tic();
    endif

    [Phi, lambda] = eig_stiff_mat(Kred, Ared, options.number_of_modes);

    if (options.verbose)
      toc();
      fprintf(stderr, "computing mode shapes ...\n");
      tic();
    endif

    n = zeros(length(comp_mat.bearing_surf.s) * length(comp_mat.bearing_surf.t), length(master_dof));

    for i=1:length(comp_mat.bearing_surf.s)
      for j=1:length(comp_mat.bearing_surf.t)
        inode = mesh.structured.inode_idx(1, i, j + comp_mat.bearing_surf.offset_t);
        node_dofs = dof_map.ndof(inode, 1:3);
        nij = mesh.nodes(inode, 1:2).';
        nij /= norm(nij);
        idxrow = (i - 1) * length(comp_mat.bearing_surf.t) + j;

        for k=1:length(nij)
          idxcol = find(master_dof == node_dofs(k));
          n(idxrow, idxcol) = nij(k);
        endfor
      endfor
    endfor

    kappa = zeros(columns(Phi), 1);

    for i=1:columns(Phi)
      kappa(i) = sqrt(1 / (Phi(:, i).' * Kred * Phi(:, i)));
    endfor

    Phi *= diag(kappa);

    KPhi = Phi.' * Kred * Phi;

    RPhi = Phi.' * Rred;
    Phin = n * Phi;

    if (isfield(options, "modal_error"))
      err_red = zeros(columns(Phi), 1);
      uref = 0;

      pred = zeros(rows(Rred), columns(Rred));
      ured = zeros(rows(Rred), columns(Rred));

      for i=1:columns(Rred)
        pred(:, i) = test_press_for_modal(bearing_dimensions, comp_mat.bearing_surf, loads, Rred, i);
        ured(:, i) = Kred \ pred(:, i);
        uref = max(uref, max(abs(n * ured(:, i))));
      endfor

      for j=1:columns(Phi)
        for i=1:columns(Rred)
          Kred2 = Phi(:, 1:j).' * Kred * Phi(:, 1:j);
          ured2 = Phi(:, 1:j) * (Kred2 \ (Phi(:, 1:j).' * pred(:, i)));
          err_red(j) += max(abs(n * (ured(:, i) - ured2))) / uref;
        endfor
        err_red(j) /= columns(Rred);
      endfor

      comp_mat.err_red = err_red;
    endif

    comp_mat.master_dof = master_dof;
    comp_mat.KPhi = KPhi;
    comp_mat.Phin = Phin;
    comp_mat.RPhi = RPhi;

    if (options.output_solution)
      U = Tred * Phi;

      for i=1:columns(U)
        U(:, i) /= max(abs(U(:, i)));
      endfor

      comp_mat.modal.def = fem_post_def_nodal(mesh, dof_map, U);

      comp_mat.modal_sol.def = fem_post_def_nodal(comp_mat.mesh, ...
                                                       comp_mat.dof_map, ...
                                                       Tred * Phi * (KPhi \ RPhi));
    endif

    comp_mat.lambda = lambda;

    if (options.verbose)
      toc();
    endif
  endif
endfunction

function [x, z] = surf_idx_to_surf_pos(bearing_dimensions, bearing_surf, s, t)
  x = bearing_dimensions.bearing_diameter * pi * (s - bearing_surf.s(1)) ...
      / (bearing_surf.s(end) - bearing_surf.s(1));

  z = bearing_dimensions.bearing_width * ((t - bearing_surf.t(1)) ...
                                         / (bearing_surf.t(end) - bearing_surf.t(1)) - 0.5);
endfunction

function pred = test_press_for_modal(bearing_dimensions, bearing_surf, loads, Rred, i)
  pred = zeros(rows(Rred), 1);

  [xi, zi] = surf_idx_to_surf_pos(bearing_dimensions, ...
                                  bearing_surf, ...
                                  loads(i).position.s, ...
                                  loads(i).position.t);

  for k=1:columns(Rred)
    [xk, zk] = surf_idx_to_surf_pos(bearing_dimensions, ...
                                    bearing_surf, ...
                                    loads(k).position.s, ...
                                    loads(k).position.t);

    if (abs(xk - xi) > 0.5 * pi * bearing_dimensions.bearing_diameter)
      xk -= sign(xk - xi) * bearing_dimensions.bearing_diameter * pi;
    endif

    dr = sqrt((xk - xi)^2 + (zk - zi)^2) / min(bearing_dimensions.bearing_width, ...
                                               pi * bearing_dimensions.bearing_diameter);
    if (dr < 0.25)
      pred += Rred(:, k);
    endif
  endfor
endfunction

function [Phi, lambda] = eig_stiff_mat(K, A, max_modes)
  opts.maxiter = int32(50000);
  opts.disp = int32(0);

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);

    [Phi, lambda, status] = eigs(K, A, max_modes, "sm", opts);

    if (status ~= 0)
      error("eigs failed with status %d", status);
    endif
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  [lambda, idx_lambda] = sort(diag(lambda), "ascend");

  max_modes = min(length(lambda), max_modes);

  lambda = lambda(1:max_modes);
  idx_lambda = idx_lambda(1:max_modes);
  Phi = Phi(:, idx_lambda);
endfunction

