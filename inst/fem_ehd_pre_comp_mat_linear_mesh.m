## Copyright (C) 2025(-2026) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{mat_ass}, @var{dof_map}, @var{cms_opt}, @var{comp_mat}, @var{bearing_surf}, @var{sol_eig}, @var{cond_info}] = fem_ehd_pre_comp_mat_linear_mesh(@var{mesh}, @var{load_case_dof}, @var{cms_opt}, @var{bearing_surf}, @var{options})
## Compute a compliance matrix suitable for MBDyn's module-hydrodynamic_plain_bearing2
##
## @var{mesh} @dots{} Finite Element mesh data structure returned from fem_cms_create.
##
## @var{mat_ass} @dots{} Assembled Finite Element matrices returned from fem_cms_create.
##
## @var{dof_map} @dots{} Degree of freedom mapping returned from fem_cms_create.
##
## @var{cms_opt} @dots{} Component mode synthesis options returned from fem_cms_create.
##
## @var{bearing_surf} @dots{} Struct array describing all hydrodynamic bearings of this mesh.
##
## @var{bearing_surf}.options.bearing_type @dots{} One of "journal", "shell".
##
## @var{bearing_surf}.options.interpolate_interface @dots{} Generate a layer of iso8 elements in order to interpolate between nodes
##
## @var{bearing_surf}.R @dots{} Orientation of the bearing surface. R(:, 3) represents the axis of the journal bearing.
##
## @var{bearing_surf}.X0 @dots{} Location of the bearing center.
##
## @var{bearing_surf}.r @dots{} Radius of the cylindrical bearing surface.
##
## @var{bearing_surf}.w @dots{} Width of the cylindrical bearing surface.
##
## @var{bearing_surf}.nodes @dots{} Node numbers of all Finite Element nodes at the bearing surface.
##
## @var{bearing_surf}.idx_load_case @dots{} Column indices of all load cases in @var{mat_ass}.R related to this bearing surface.
##
## @var{options}.elem_type @dots{} Element type for the surface mesh at the bearing
##
## @seealso{fem_ehd_pre_comp_mat_load_case, fem_ehd_pre_comp_mat_export, fem_ehd_pre_comp_mat_plot}
## @end deftypefn

function [mesh, mat_ass_itf, dof_map_itf, cms_opt, comp_mat, bearing_surf, sol_eig, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, options)
  if (nargin < 4 || nargin > 5 || nargout > 8)
    print_usage();
  endif

  if (~isfield(cms_opt, "solver"))
    cms_opt.solver = "";
  endif

  if (~isfield(cms_opt, "number_of_threads"))
    cms_opt.number_of_threads = int32(1);
  endif

  if (nargin < 5)
    options = struct();
  endif

  if (~isfield(options, "elem_type"))
    options.elem_type = "tria6";
  endif

  if (~isfield(cms_opt, "floating_frame"))
    cms_opt.floating_frame = false;
  endif

  if (~isfield(cms_opt, "svd_threshold"))
    cms_opt.svd_threshold = 1e-10;
  endif

  if (~isfield(cms_opt, "lambda_threshold"))
    cms_opt.lambda_threshold = 1e-10;
  endif

  if (~isfield(cms_opt, "max_cond_D"))
    cms_opt.max_cond_D = 1e8;
  endif

  if (~isfield(cms_opt, "tol_gamma_rel"))
    cms_opt.tol_gamma_rel = 1e-6;
  endif

  if (~isfield(cms_opt, "tol_gamma_abs"))
    cms_opt.tol_gamma_abs = 1e-12;
  endif

  if (~isfield(cms_opt, "verbose"))
    cms_opt.verbose = int32(1);
  endif

  if (~isfield(cms_opt.nodes.interfaces, "include_rigid_body_modes"))
    for i=1:numel(cms_opt.nodes.interfaces)
      cms_opt.nodes.interfaces(i).include_rigid_body_modes = true;
    endfor
  endif

  for i=1:numel(bearing_surf)
    if (~isfield(bearing_surf(i).options, "interpolate_interface"))
      bearing_surf(i).options.interpolate_interface = false;
    endif
  endfor

  cms_opt.solver = fem_sol_select(true, cms_opt.solver);

  options.interpolate_interface = [[bearing_surf.options].interpolate_interface];

  [bearing_surf] = fem_ehd_pre_comp_mat_grid(mesh, bearing_surf, options);

  [mesh, load_case_dof, bearing_surf] = fem_ehd_pre_comp_mat_gen_mesh(mesh, load_case_dof, bearing_surf);

  mesh = fem_ehd_pre_comp_mat_gen_constr(mesh, bearing_surf, options);

  [mesh, load_case_itf] = fem_ehd_pre_comp_mat_gen_loads_cms(mesh, load_case_dof, cms_opt);

  [Phi_p, kappa_p, mat_ass_eig_p, dof_map_eig_p, load_case_eig_p, bearing_surf] = fem_ehd_pre_comp_mat_eig_p(mesh, load_case_dof, bearing_surf, cms_opt, options);

  [Phi_s, Phi_n, lambda_n, mat_ass_itf, dof_map_itf] = fem_ehd_pre_comp_mat_modes_itf(mesh, load_case_dof, load_case_itf, cms_opt);

  num_modes_cb = columns(Phi_s) + columns(Phi_n);

  if (cms_opt.verbose)
    printf("%d Craig Bampton modes in total\n", num_modes_cb);
    printf("%d pressure modes in total\n", columns(Phi_p));
  endif

  [mat_ass_itf, bearing_surf] = fem_ehd_pre_comp_mat_modes_combine(mesh, dof_map_itf, dof_map_eig_p, mat_ass_itf, cms_opt, bearing_surf, Phi_s, Phi_n, Phi_p);

  [comp_mat] = fem_ehd_pre_comp_mat_gen(mesh, dof_map_itf, mat_ass_itf, load_case_itf, cms_opt, bearing_surf);

  [mat_ass_itf, comp_mat, cond_info] = fem_ehd_pre_comp_mat_filter_cond(mat_ass_itf, comp_mat, cms_opt, num_modes_cb);

  [mat_ass_itf, comp_mat, cond_info] = fem_ehd_pre_comp_mat_filter_svd(mat_ass_itf, comp_mat, cms_opt, num_modes_cb, cond_info);

  [mat_ass_itf, comp_mat, cond_info] = fem_ehd_pre_comp_mat_filter_lambda(mat_ass_itf, dof_map_itf, comp_mat, cms_opt, num_modes_cb, cond_info);

  [mat_ass_itf, comp_mat] = fem_ehd_pre_comp_mat_filter_chol(mat_ass_itf, dof_map_itf, comp_mat);

  [mat_ass_itf, sol_eig] = fem_ehd_pre_comp_mat_gen_cms(mesh, dof_map_itf, mat_ass_itf, load_case_itf, lambda_n, kappa_p, cms_opt);

  for i=1:numel(comp_mat)
    comp_mat(i).E = comp_mat(i).D.' * diag(comp_mat(i).A) * comp_mat(i).reference_pressure;
  endfor

  if (nargout >= 8)
    cond_info = fem_ehd_pre_comp_mat_cond(mat_ass_itf, comp_mat, cond_info);
  endif
endfunction

function [mat_ass_itf, comp_mat] = fem_ehd_pre_comp_mat_filter_chol(mat_ass_itf, dof_map_itf, comp_mat);
  Msym = fem_mat_sym(mat_ass_itf.M)(dof_map_itf.idx_node, dof_map_itf.idx_node);

  Mred = fem_cms_matrix_trans(mat_ass_itf.Tred, Msym, "Lower");

  L = chol(Mred, "lower");

  mat_ass_itf.Tred /= L.';

  for i=1:numel(comp_mat)
    comp_mat(i).D /= L.';
  endfor
endfunction

function [mat_ass_itf, comp_mat, cond_info] = fem_ehd_pre_comp_mat_filter_cond(mat_ass_itf, comp_mat, cms_opt, num_modes_cb)
  keep_modes = false(1, columns(mat_ass_itf.Tred));

  keep_modes(1:num_modes_cb) = true;

  num_modes_rejected_gamma = int32(0);
  num_modes_rejected_cond = int32(0);

  for i=1:numel(comp_mat)
    [D, A] = fem_ehd_comp_mat_tot(mat_ass_itf, comp_mat(i));

    selected = false(1, numel(comp_mat(i).mode_idx));
    gamma = zeros(1, numel(comp_mat(i).mode_idx));

    for j=1:numel(comp_mat(i).mode_idx)
      d_j = D(:, comp_mat(i).mode_idx(j));

      if (isempty(selected))
        selected(j) = true;
        gamma(j) = (d_j.' * diag(A) * d_j) / (d_j.' * d_j);
        continue;
      endif

      D_old = D(:, selected);

      ## A-weighted projection
      G   = D_old.' * diag(A) * D_old;
      rhs = D_old.' * diag(A) * d_j;

      coeff = pinv(G) * rhs;

      d_j_orth = d_j - D_old * coeff;

      gamma(j) = (d_j_orth.' * diag(A) * d_j_orth) / (d_j.' * d_j);

      if (gamma(j) > cms_opt.tol_gamma_rel * median(gamma(1:j - 1)) && gamma(j) > cms_opt.tol_gamma_abs)
        selected(j) = true;
      else
        ++num_modes_rejected_gamma;
      endif
    endfor

    selected = find(selected);

    cond_info.comp_mat(i).gamma = gamma;
    cond_info.comp_mat(i).gamma_range = selected;

    cond_info.comp_mat(i).cond_value = zeros(size(selected));

    for j=1:numel(selected)
      Dj = D(:, comp_mat(i).mode_idx(selected(1:j)));

      cond_info.comp_mat(i).cond_value(j) = cond(Dj.' * diag(A) * Dj);

      if (cond_info.comp_mat(i).cond_value(j) > cms_opt.max_cond_D)
        ++num_modes_rejected_cond;
      else
        keep_modes(comp_mat(i).mode_idx(selected(j))) = true;
      endif
    endfor

    cond_info.comp_mat(i).cond_range = find(keep_modes(comp_mat(i).mode_idx(selected)));
  endfor

  if (cms_opt.verbose)
    printf("keeping %d of %d modes (gamma > %.2e)\n", numel(keep_modes) - num_modes_cb - num_modes_rejected_gamma, numel(keep_modes) - num_modes_cb, cms_opt.tol_gamma_rel);
    printf("keeping %d of %d modes (cond < %.2e)\n", numel(keep_modes) - num_modes_cb - num_modes_rejected_gamma - num_modes_rejected_cond, numel(keep_modes) - num_modes_cb - num_modes_rejected_gamma, cms_opt.max_cond_D);
  endif

  mat_ass_itf.Tred = mat_ass_itf.Tred(:, keep_modes);

  mode_idx_new = zeros(size(keep_modes), "int32");
  mode_idx_new(find(keep_modes)) = (1:sum(keep_modes))(:);

  for i=1:numel(comp_mat)
    comp_mat(i).D = comp_mat(i).D(:, keep_modes);
    comp_mat(i).mode_idx = mode_idx_new(comp_mat(i).mode_idx);
    comp_mat(i).mode_idx = comp_mat(i).mode_idx(comp_mat(i).mode_idx > 0);
  endfor
endfunction

function [mat_ass_itf, comp_mat, cond_info] = fem_ehd_pre_comp_mat_filter_lambda(mat_ass_itf, dof_map_itf, comp_mat, cms_opt, num_modes_cb, cond_info)
  [D, A] = fem_ehd_comp_mat_tot(mat_ass_itf, comp_mat);

  G = fem_cms_matrix_trans(D(:, num_modes_cb + 1:end), diag(A), "Lower");

  Msym = fem_mat_sym(mat_ass_itf.M)(dof_map_itf.idx_node, dof_map_itf.idx_node);

  Mred = fem_cms_matrix_trans(mat_ass_itf.Tred(:, num_modes_cb + 1:end), Msym, "Lower");

  [V, lambda] = eig(G, Mred, "vector", "chol");

  idx_hydro = lambda > cms_opt.lambda_threshold * max(lambda);

  cond_info.lambda = lambda;
  cond_info.lambda_range = find(idx_hydro);

  if (cms_opt.verbose)
    fprintf(stderr, "keeping %d of %d modes (lambda > %e)\n", sum(idx_hydro), numel(idx_hydro), cms_opt.lambda_threshold);
  endif

  V = V(:, idx_hydro);

  V *= diag(lambda(idx_hydro).^(-1/2));

  mat_ass_itf.Tred = [mat_ass_itf.Tred(:, 1:num_modes_cb), mat_ass_itf.Tred(:, num_modes_cb + 1:end) * V];

  for i=1:numel(comp_mat)
    comp_mat(i).D = [comp_mat(i).D(:, 1:num_modes_cb), comp_mat(i).D(:, num_modes_cb + 1:end) * V];
  endfor
endfunction

function [mat_ass_itf, bearing_surf] = fem_ehd_pre_comp_mat_modes_combine(mesh, dof_map_itf, dof_map_eig_p, mat_ass_itf, cms_opt, bearing_surf, Phi_s, Phi_n, Phi_p)
  Phi_p = Phi_p(dof_map_eig_p.idx_node, :);
  Phi_n = Phi_n(dof_map_itf.idx_node, :);
  Phi_s = Phi_s(dof_map_itf.idx_node, :);

  Msym = fem_mat_sym(mat_ass_itf.M)(dof_map_itf.idx_node, dof_map_itf.idx_node);

  Phi_sn = [Phi_s, Phi_n];

  Phi_p -= Phi_sn * ((Phi_sn.' * (Msym * Phi_sn)) \ (Phi_sn.' * (Msym * Phi_p)));

  mat_ass_itf.Tred = [Phi_sn, Phi_p];

  for i=1:numel(bearing_surf)
    bearing_surf(i).options.mode_idx += columns(Phi_sn);
  endfor

  if (cms_opt.floating_frame)
    Phi_rb = fem_ehd_pre_comp_mat_gen_rigid_body_mode(mesh, dof_map_itf, cms_opt.nodes.modal.number);
    mat_ass_itf.Tred -= Phi_rb * (Phi_rb \ mat_ass_itf.Tred);
  endif
endfunction

function [D, A] = fem_ehd_comp_mat_tot(mat_ass_itf, comp_mat)
  num_rows_D = int32(0);

  for i=1:numel(comp_mat)
    nz = numel(comp_mat(i).bearing_surf.grid_z);
    num_rows_D += rows(comp_mat(i).D) - nz;
  endfor

  D = zeros(num_rows_D, columns(mat_ass_itf.Tred));
  A = zeros(1, num_rows_D);

  num_rows_D = int32(0);

  for i=1:numel(comp_mat)
    nz = numel(comp_mat(i).bearing_surf.grid_z);
    D(num_rows_D + (1:rows(comp_mat(i).D) - nz), :) = comp_mat(i).D(1:end - nz, :);
    A(num_rows_D + (1:rows(comp_mat(i).D) - nz)) = comp_mat(i).A(1:end - nz);
    num_rows_D += rows(comp_mat(i).D) - nz;
  endfor
endfunction

function cond_info = fem_ehd_pre_comp_mat_cond(mat_ass_itf, comp_mat, cond_info)
  [D, A] = fem_ehd_comp_mat_tot(mat_ass_itf, comp_mat);

  G = D.' * diag(A) * D;

  cond_info.D_rank = rank(D);
  cond_info.D_size = size(D);
  cond_info.D_cond = cond(G);
  cond_info.eta = diag(G);
  cond_info.eta /= max(cond_info.eta);
  cond_info.eta = sort(cond_info.eta, "descend");
endfunction

function [mat_ass_itf, sol_eig] = fem_ehd_pre_comp_mat_gen_cms(mesh, dof_map_itf, mat_ass_itf, load_case_itf, lambda_n, kappa_p, cms_opt)
  Msym = fem_mat_sym(mat_ass_itf.M);

  switch (mat_ass_itf.mat_info.mat_type(3))
    case FEM_MAT_STIFFNESS
      Ksym = mat_ass_itf.K;
    case {FEM_MAT_STIFFNESS_SYM, FEM_MAT_STIFFNESS_SYM_L}
      Ksym = fem_mat_sym(mat_ass_itf.K);
    otherwise
      error("unknown matrix type for mat_ass_itf.K");
  endswitch

  Dsym = fem_mat_sym(mat_ass_itf.D);

  mat_ass_itf.Mred = fem_cms_matrix_trans(mat_ass_itf.Tred, Msym(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");
  mat_ass_itf.Kred = fem_cms_matrix_trans(mat_ass_itf.Tred, Ksym(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");
  mat_ass_itf.Dred = fem_cms_matrix_trans(mat_ass_itf.Tred, Dsym(dof_map_itf.idx_node, dof_map_itf.idx_node), "Lower");

  Phi = zeros(dof_map_itf.totdof, columns(mat_ass_itf.Tred));

  Phi(dof_map_itf.idx_node, :) = mat_ass_itf.Tred;

  sol_eig.def = fem_post_def_nodal(mesh, dof_map_itf, Phi);

  clear Phi;

  ## Needed for FEM_MAT_INERTIA_INV4 and FEM_MAT_INERTIA_INV8
  mesh.nodes -= mesh.nodes(cms_opt.nodes.modal.number, :);

  [mat_ass_itf.Inv3, ...
   mat_ass_itf.Inv4, ...
   mat_ass_itf.Inv5, ...
   mat_ass_itf.Inv8, ...
   mat_ass_itf.Inv9] = fem_ass_matrix(mesh, ...
                                      dof_map_itf, ...
                                      [FEM_MAT_INERTIA_INV3, ...
                                       FEM_MAT_INERTIA_INV4, ...
                                       FEM_MAT_INERTIA_INV5, ...
                                       FEM_MAT_INERTIA_INV8, ...
                                       FEM_MAT_INERTIA_INV9], ...
                                      load_case_itf, ...
                                      sol_eig);
endfunction

function [comp_mat] = fem_ehd_pre_comp_mat_gen(mesh, dof_map_itf, mat_ass_itf, load_case_itf, cms_opt, bearing_surf)
  empty_cell = cell(1, numel(bearing_surf));

  comp_mat = struct("C", empty_cell, ...
                    "D", empty_cell, ...
                    "E", empty_cell, ...
                    "xi", empty_cell, ...
                    "zi", empty_cell, ...
                    "dX", empty_cell, ...
                    "dR", empty_cell, ...
                    "bearing_surf", empty_cell, ...
                    "bearing_dimensions", empty_cell, ...
                    "reference_pressure", empty_cell, ...
                    "A", empty_cell, ...
                    "mode_idx", empty_cell);

  for i=1:numel(bearing_surf)
    switch (bearing_surf(i).options.bearing_type)
      case "shell"
        s = 1;
      case "journal"
        s = -1;
      otherwise
        error("invalid value for parameter bearing_type\"%s\"", bearing_surf(i).options.bearing_type);
    endswitch

    comp_mat(i).bearing_dimensions.bearing_diameter = 2 * bearing_surf(i).r;
    comp_mat(i).bearing_dimensions.bearing_width = bearing_surf(i).w;
    comp_mat(i).reference_pressure = bearing_surf(i).options.reference_pressure;
    comp_mat(i).bearing_surf.grid_x = bearing_surf(i).grid_x;
    comp_mat(i).bearing_surf.grid_z = bearing_surf(i).grid_z;
    comp_mat(i).dX = bearing_surf(i).X0 - mesh.nodes(cms_opt.nodes.modal.number, 1:3).';
    comp_mat(i).dR = bearing_surf(i).R;
    comp_mat(i).mode_idx = bearing_surf(i).options.mode_idx;

    bs_nodes = bearing_surf(i).nodes;

    X = bearing_surf(i).R.' * (mesh.nodes(bs_nodes, 1:3).' - bearing_surf(i).X0);

    n = s * X(1:2, :);
    r = norm(n, "cols");
    n ./= r;

    comp_mat(i).xi = mod(atan2(X(2, :), X(1, :)), 2 * pi) * bearing_surf(i).r;
    comp_mat(i).zi = X(3, :);

    comp_mat(i).D = zeros(numel(bs_nodes), columns(mat_ass_itf.Tred));

    dof_idx = dof_map_itf.ndof(bs_nodes, 1:3);

    Ustatn = zeros(rows(dof_idx), columns(mat_ass_itf.Tred), 3);

    for k=1:columns(dof_idx)
      Ustatn(:, :, k) = mat_ass_itf.Tred(dof_idx(:, k), :);
    endfor

    for k=1:rows(n)
      for l=1:columns(dof_idx)
        comp_mat(i).D += diag(n(k, :)) * Ustatn(:, :, l) * bearing_surf(i).R(l, k);
      endfor
    endfor

    Ai = repmat((bearing_surf(i).grid_x(2) - bearing_surf(i).grid_x(1)) * (bearing_surf(i).grid_z(2) - bearing_surf(i).grid_z(1)), numel(bs_nodes), 1);

    Ai(1:numel(bearing_surf(i).grid_z):end) /= 2;
    Ai(numel(bearing_surf(i).grid_z):numel(bearing_surf(i).grid_z):end) /= 2;
    Ai((numel(bearing_surf(i).grid_x) - 1) * numel(bearing_surf(i).grid_z) + 1:end) = 0;

    comp_mat(i).A = Ai;
  endfor
endfunction

function [mat_ass_itf, comp_mat, cond_info] = fem_ehd_pre_comp_mat_filter_svd(mat_ass_itf, comp_mat, cms_opt, num_modes_cb, cond_info)
  num_rows_D = int32(0);

  for i=1:numel(comp_mat)
    nz = numel(comp_mat(i).bearing_surf.grid_z);
    num_rows_D += rows(comp_mat(i).D) - nz;
  endfor

  Dcb = zeros(num_rows_D, num_modes_cb);
  Dp = zeros(num_rows_D, columns(mat_ass_itf.Tred) - num_modes_cb);

  num_rows_D = int32(0);

  for i=1:numel(comp_mat)
    nz = numel(comp_mat(i).bearing_surf.grid_z);
    Dcb(num_rows_D + (1:rows(comp_mat(i).D) - nz), :) = diag(comp_mat(i).A(1:end - nz).^(1/2)) * comp_mat(i).D(1:end - nz, 1:num_modes_cb);
    Dp(num_rows_D + (1:rows(comp_mat(i).D) - nz), :) = diag(comp_mat(i).A(1:end - nz).^(1/2)) * comp_mat(i).D(1:end - nz, num_modes_cb + 1:end);
    num_rows_D += rows(comp_mat(i).D) - nz;
  endfor

  [Q, ~] = qr(Dcb, "econ");

  Pcb = Q * Q.';

  Dp_tilde = (eye(size(Q, 1)) - Pcb) * Dp;

  [U, S, V] = svd(Dp_tilde, 'econ');

  cond_info.S = diag(S);

  idx_keep = sum(diag(S) > cms_opt.svd_threshold * max(abs(diag(S))));

  cond_info.S_range = find(idx_keep);

  if (cms_opt.verbose)
    fprintf(stderr, "keeping %d of %d modes (S > %e)\n", idx_keep, columns(V), cms_opt.svd_threshold);
  endif

  V = V(:, 1:idx_keep);

  mat_ass_itf.Tred = [mat_ass_itf.Tred(:, 1:num_modes_cb), mat_ass_itf.Tred(:, num_modes_cb + 1:end) * V];

  for i=1:numel(comp_mat)
    comp_mat(i).D = [comp_mat(i).D(:, 1:num_modes_cb), comp_mat(i).D(:, num_modes_cb + 1:end) * V];
  endfor
endfunction

function [Phi, kappa, mat_ass_eig_p, dof_map_eig_p, load_case_eig_p, bearing_surf] = fem_ehd_pre_comp_mat_eig_p(mesh, load_case_dof, bearing_surf, cms_opt, options)
  dof_map_eig_p = fem_ass_dof_map(mesh, load_case_dof);
  dof_map_eig_p.parallel.threads_ass = cms_opt.number_of_threads;

  mat_type_stiffness = FEM_MAT_STIFFNESS_SYM_L;

  switch (cms_opt.solver)
    case {"umfpack", "lu", "mldivide"}
      mat_type_stiffness = FEM_MAT_STIFFNESS; ## enforce full matrix for unsymmetric solvers
  endswitch

  inum_elem_press = int32(0);

  for i=1:numel(bearing_surf)
    inum_elem_press += numel(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements);
  endfor

  elements = zeros(inum_elem_press, columns(getfield(mesh.elements, options.elem_type)), "int32");

  inum_elem_press = int32(0);
  ielem_idx = zeros(numel(bearing_surf), 2, "int32");

  for i=1:numel(bearing_surf)
    elem_idx = getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements;
    idx_slot = inum_elem_press + (1:numel(elem_idx));
    elements(idx_slot, :) = getfield(mesh.elements, options.elem_type)(elem_idx, :);
    ielem_idx(i, :) = idx_slot([1, end]);
    inum_elem_press += numel(elem_idx);
  endfor

  load_case_eig_p = fem_pre_load_case_create_empty(numel(bearing_surf));

  for i=1:numel(bearing_surf)
    press_elem.p = zeros(size(elements));
    press_elem.p(ielem_idx(i, 1):ielem_idx(i, end), :) = 1;
    press_elem.elements = elements;
    load_case_eig_p(i).pressure = setfield(load_case_eig_p(i).pressure, options.elem_type, press_elem);
  endfor

  [mat_ass_eig_p.K, ...
   mat_ass_eig_p.R, ...
   mat_ass_eig_p.mat_info, ...
   mat_ass_eig_p.mesh_info] = fem_ass_matrix(mesh, ...
                                             dof_map_eig_p, ...
                                             [mat_type_stiffness, ...
                                              FEM_VEC_LOAD_CONSISTENT], ...
                                             load_case_eig_p);

  diagA = zeros(rows(mesh.nodes), numel(bearing_surf));

  for i=1:columns(dof_map_eig_p.ndof)
    dof_idx = dof_map_eig_p.ndof(:, i);
    idx_act_dof = find(dof_idx > 0);
    diagA(idx_act_dof, :) += mat_ass_eig_p.R(dof_idx(idx_act_dof), :).^2;
  endfor

  diagA = sqrt(diagA);

  num_modes = int32(0);

  for i=1:numel(bearing_surf)
    num_modes += bearing_surf(i).options.number_of_modes;
  endfor

  Phi = zeros(dof_map_eig_p.totdof, num_modes);
  kappa = zeros(1, num_modes);
  num_modes = int32(0);

  Kfact = fem_sol_factor(mat_ass_eig_p.K, cms_opt);

  for i=1:numel(bearing_surf)
    if (bearing_surf(i).options.number_of_modes)
      Ap = zeros(dof_map_eig_p.totdof, 1);

      for j=1:columns(dof_map_eig_p.ndof)
        dof_idx = dof_map_eig_p.ndof(:, j);
        idx_act_dof = find(dof_idx > 0);
        Ap(dof_idx(idx_act_dof)) = diagA(idx_act_dof, i);
      endfor

      Ap = diag(Ap);

      rndstate = rand("state");

      unwind_protect
        rand("seed", 0);

        oper = cell(1, 2);
        oper{1} = @(x) Ap * x;
        oper{2} = @(x) Kfact \ x;
        sigma = 0;

        opts.disp = 0;
        opts.maxit = 50000;
        opts.tol = 0;

        num_modes_i = bearing_surf(i).options.number_of_modes;

        [Phi_i, kappa_i] = eig_sym(oper, columns(mat_ass_eig_p.K), num_modes_i, sigma, opts);

      unwind_protect_cleanup
        rand("state", rndstate);
      end_unwind_protect

      clear oper Ap;

      kappa_i = diag(kappa_i);
    else
      kappa_i = [];
      Phi_i = [];
    endif

    bearing_surf(i).options.mode_idx = num_modes + (1:bearing_surf(i).options.number_of_modes);
    kappa(bearing_surf(i).options.mode_idx) = kappa_i;
    Phi(:, bearing_surf(i).options.mode_idx) = Phi_i;
    num_modes += bearing_surf(i).options.number_of_modes;
  endfor
endfunction

function [Phi_s, Phi_n, lambda, mat_ass_itf, dof_map_itf] = fem_ehd_pre_comp_mat_modes_itf(mesh, load_case_dof, load_case_itf, cms_opt)
  dof_map_itf = fem_ass_dof_map(mesh, load_case_dof);

  dof_map_itf.parallel.threads_ass = cms_opt.number_of_threads;

  mat_type_stiffness = FEM_MAT_STIFFNESS_SYM_L;

  if (~isfield(cms_opt, "solver"))
    cms_opt.solver = "";
  endif

  cms_opt.solver = fem_sol_select(true, cms_opt.solver);

  switch (cms_opt.solver)
    case {"umfpack", "lu", "mldivide"}
      mat_type_stiffness = FEM_MAT_STIFFNESS; ## enforce full matrix for unsymmetric solvers
  endswitch

  [mat_ass_itf.M, ...
   mat_ass_itf.D, ...
   mat_ass_itf.K, ...
   mat_ass_itf.R, ...
   mat_ass_itf.dm, ...
   mat_ass_itf.S, ...
   mat_ass_itf.J, ...
   mat_ass_itf.mat_info, ...
   mat_ass_itf.mesh_info] = fem_ass_matrix(mesh, ...
                                           dof_map_itf, ...
                                           [FEM_MAT_MASS_SYM_L, ...
                                            FEM_MAT_DAMPING_SYM_L, ...
                                            mat_type_stiffness, ...
                                            FEM_VEC_LOAD_CONSISTENT, ...
                                            FEM_SCA_TOT_MASS, ...
                                            FEM_VEC_INERTIA_M1, ...
                                            FEM_MAT_INERTIA_J], ...
                                           load_case_itf);

  Kfact = fem_sol_factor(mat_ass_itf.K, cms_opt);
  Msym = fem_mat_sym(mat_ass_itf.M);

  Phi_s = Kfact \ mat_ass_itf.R;

  [Phi_n, lambda] = fem_ehd_pre_comp_mat_eigs(Msym, Kfact, cms_opt);
endfunction

function [Phi, lambda] = fem_ehd_pre_comp_mat_eigs(Msym, Kfact, cms_opt)
  if (cms_opt.modes.number == 0)
    Phi = zeros(columns(Msym), 0);
    lambda = [];
    return;
  endif

  opts.disp = 0;
  opts.maxit = 50000;
  opts.tol = 0;

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);

    SIGMA = 0;
    op{1} = @(x) Msym * x;
    op{2} = @(x) Kfact \ x;

    [Phi, mu] = eig_sym(op, columns(Msym), cms_opt.modes.number, SIGMA, opts);
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  lambda = sqrt(-diag(mu)).';
endfunction

function [mesh, load_case] = fem_ehd_pre_comp_mat_gen_loads_cms(mesh, load_case_dof, cms_opt)
  if (~isfield(mesh.elements, "joints"))
    empty_cell = cell(0, 0);
    mesh.elements.joints = struct("C", empty_cell, "nodes", empty_cell);
  endif

  joint_modal.C = eye(6)(~load_case_dof.locked_dof(cms_opt.nodes.modal.number, :), :);
  joint_modal.nodes = cms_opt.nodes.modal.number;

  empty_cell = cell(1, numel(cms_opt.nodes.interfaces));
  joints_itf = struct("C", empty_cell, "nodes", empty_cell);
  num_modes_itf = int32(0);

  for i=1:numel(cms_opt.nodes.interfaces)
    joints_itf(i).nodes = cms_opt.nodes.interfaces(i).number;
    joints_itf(i).C = eye(6)(~load_case_dof.locked_dof(cms_opt.nodes.interfaces(i).number, :), :);

    if (cms_opt.nodes.interfaces(i).include_rigid_body_modes)
      num_modes_itf += rows(joints_itf(i).C);
    endif
  endfor

  if (~isempty(joint_modal.C))
    mesh.elements.joints = [mesh.elements.joints, joint_modal];
  endif

  idx_joint_itf = numel(mesh.elements.joints) + 1;

  mesh.elements.joints = [mesh.elements.joints, joints_itf];

  load_case = fem_pre_load_case_create_empty(num_modes_itf);

  for i=1:numel(load_case)
    load_case(i).joints = struct("U", cellfun(@(C) zeros(rows(C), 1), {mesh.elements.joints.C}, "UniformOutput", false));
  endfor

  idx_mode_itf = int32(0);

  for i=1:numel(cms_opt.nodes.interfaces)
    if (cms_opt.nodes.interfaces(i).include_rigid_body_modes)
      for j=1:rows(mesh.elements.joints(idx_joint_itf + i - 1).C)
        load_case(++idx_mode_itf).joints(idx_joint_itf + i - 1).U(j) = 1;
      endfor
    endif
  endfor
endfunction

function [mesh] = fem_ehd_pre_comp_mat_gen_constr(mesh, bearing_surf, options)
  switch (options.elem_type)
    case "tria6"
      elem_type_sfn = "sfncon6";
    case "tria6h"
      elem_type_sfn = "sfncon6h";
    case "tria10"
      elem_type_sfn = "sfncon10";
    case "quad8"
      elem_type_sfn = "sfncon8";
    case "quad8r"
      elem_type_sfn = "sfncon8r";
    case "quad9"
      elem_type_sfn = "sfncon9";
    case "iso4"
      elem_type_sfn = "sfncon4";
    otherwise
      error("unknown value for options.elem_type=\"%s\"", options.elem_type);
  endswitch

  if (~isfield(mesh.elements, elem_type_sfn))
    empty_cell = cell(0, 0);
    mesh.elements = setfield(mesh.elements, ...
                             elem_type_sfn, ...
                             struct("master", empty_cell, ...
                                    "slave", empty_cell, ...
                                    "maxdist", empty_cell, ...
                                    "constraint", empty_cell));
  endif

  if (~isfield(mesh.elements, "joints"))
    empty_cell = cell(0, 0);
    mesh.elements.joints = struct("C", empty_cell, "nodes", empty_cell);
  endif

  for i=1:numel(bearing_surf)
    sfncon.master = getfield(mesh.elements, options.elem_type)(getfield(mesh.groups, options.elem_type)(bearing_surf(i).group_idx).elements, :);
    sfncon.slave = bearing_surf(i).nodes(:, 1:end - 1)(:);
    sfncon.maxdist = bearing_surf(i).absolute_tolerance + bearing_surf(i).relative_tolerance * bearing_surf(i).r;
    sfncon.constraint = FEM_CT_FIXED;

    if (bearing_surf(i).options.interpolate_interface)
      joints = fem_pre_mesh_constr_surf_to_node(mesh.nodes, struct(elem_type_sfn, sfncon), FEM_DO_STRUCTURAL).joints;
      mesh.elements.joints = [mesh.elements.joints, joints];
    else
      mesh.elements = setfield(mesh.elements, elem_type_sfn, [getfield(mesh.elements, elem_type_sfn), sfncon]);
    endif
  endfor

  if (isempty(getfield(mesh.elements, elem_type_sfn)))
    mesh.elements = rmfield(mesh.elements, elem_type_sfn);
  endif

  if (isempty(mesh.elements.joints))
    mesh.elements = rmfield(mesh.elements, "joints");
  endif
endfunction

function [mesh, load_case_dof, bearing_surf] = fem_ehd_pre_comp_mat_gen_mesh(mesh, load_case_dof, bearing_surf)
  if (~isfield(mesh.elements, "iso4"))
    mesh.elements.iso4 = zeros(0, 4, "int32");
  endif

  if (~isfield(mesh.elements, "iso8"))
    mesh.elements.iso8 = zeros(0, 8, "int32");
  endif

  if (~isfield(mesh.groups, "iso4"))
    empty_cell = cell(0, 0);
    mesh.groups.iso4 = struct("id", empty_cell, "name", empty_cell, "elements", empty_cell, "nodes", empty_cell);
  endif

  if (~isfield(mesh.groups, "iso8"))
    empty_cell = cell(0, 0);
    mesh.groups.iso8 = struct("id", empty_cell, "name", empty_cell, "elements", empty_cell, "nodes", empty_cell);
  endif

  if (~isfield(mesh.materials, "iso8"))
    mesh.materials.iso8 = zeros(0, 1, "int32");
  endif

  for i=1:numel(bearing_surf)
    X0_i = bearing_surf(i).X0;
    R_i = bearing_surf(i).R;

    [x_S, z_S] = meshgrid(bearing_surf(i).grid_x, bearing_surf(i).grid_z);

    Phi_i = x_S / bearing_surf(i).r;
    x_i = cos(Phi_i) * bearing_surf(i).r;
    y_i = sin(Phi_i) * bearing_surf(i).r;
    X_i = zeros(rows(x_S) * columns(x_S), 3);

    for j=1:rows(x_S)
      X_i((j - 1) * columns(x_S) + (1:columns(x_S)), :) = (X0_i + R_i * [x_i(j, :); y_i(j, :); z_S(j, :)]).';
    endfor

    node_idx = zeros(rows(x_S), columns(x_S), "int32");

    for j=1:rows(x_S)
      node_idx(j, :) = (j - 1) * columns(x_S) + [1:columns(x_S) - 1, 1];
    endfor

    elem_iso4 = zeros((columns(x_S) - 1) * (rows(x_S) - 1), 4, "int32");

    for j=1:rows(x_S) - 1
      elem_iso4((j - 1) * (columns(x_S) - 1) + (1:columns(x_S) - 1), :) = [node_idx(j + 1, (1:end - 1) + 1);
                                                                           node_idx(j + 1, (1:end - 1));
                                                                           node_idx(j, (1:end - 1));
                                                                           node_idx(j, (1:end - 1) + 1)].';
    endfor

    switch (bearing_surf(i).options.bearing_type)
      case "shell"
        elem_iso4 = elem_iso4(:, [end:-1:1]);
    endswitch

    elem_iso4 += rows(mesh.nodes);
    node_idx += rows(mesh.nodes);

    grp_iso4.id = bearing_surf(i).group_id_interface;
    grp_iso4.name = bearing_surf(i).name;
    grp_iso4.elements = rows(mesh.elements.iso4) + int32(1:rows(elem_iso4));
    grp_iso4.nodes = unique(elem_iso4).';

    mesh.elements.iso4 = [mesh.elements.iso4;
                          elem_iso4];

    mesh.nodes = [mesh.nodes;
                  [X_i, zeros(rows(X_i), 3)]];

    mesh.groups.iso4 = [mesh.groups.iso4, grp_iso4];

    bearing_surf(i).group_idx_interface = numel(mesh.groups.iso4);
    bearing_surf(i).nodes = node_idx;
  endfor

  for i=1:numel(bearing_surf)
    if (bearing_surf(i).options.interpolate_interface)
      h = mean([x_S(1, 2) - x_S(1, 1), z_S(2, 1) - z_S(1, 1)]);

      [mesh.nodes, elem_iso8] = fem_pre_mesh_extrude_surf(mesh, "iso8", bearing_surf(i).group_id_interface, h);

      grp_iso8.id = bearing_surf(i).group_id_interface;
      grp_iso8.name = bearing_surf(i).name;
      grp_iso8.elements = rows(mesh.elements.iso8) + int32(1:rows(elem_iso8));
      grp_iso8.nodes = unique(elem_iso8).';

      mesh.elements.iso8 = [mesh.elements.iso8;
                            elem_iso8];

      mesh.materials.iso8 = [mesh.materials.iso8; repmat(bearing_surf(i).material_id_interface, rows(elem_iso8), 1)];

      mesh.groups.iso8 = [mesh.groups.iso8, grp_iso8];
    endif
  endfor

  load_case_dof.locked_dof((end + 1):rows(mesh.nodes), :) = false;

  if (isempty(mesh.elements.iso8))
    mesh.elements = rmfield(mesh.elements, "iso8");
    mesh.groups = rmfield(mesh.groups, "iso8");
  endif
endfunction

function Phi_rb = fem_ehd_pre_comp_mat_gen_rigid_body_mode(mesh, dof_map, ref_node_idx)
  Arb = repmat(eye(6), rows(mesh.nodes), 1);

  l = mesh.nodes(:, 1:3) - mesh.nodes(ref_node_idx, 1:3);

  ir = int32([3, 2;
              1, 3;
              2, 1]);

  ic = int32([2, 3;
              3, 1;
              1, 2]);

  dc = [1, -1];

  for j=1:3
    for k=1:2
      Arb(ir(j, k):6:end, ic(j, k) + 3) = -dc(k) * l(:, j);
    endfor
  endfor

  Phi_rb = zeros(dof_map.totdof, columns(Arb));

  for j=1:columns(dof_map.ndof)
    idx_act_dof = find(dof_map.ndof(:, j) > 0);
    Phi_rb(dof_map.ndof(idx_act_dof, j), :) = Arb(6 * (idx_act_dof - 1) + j, :);
  endfor

  Phi_rb = Phi_rb(dof_map.idx_node, :);
endfunction

%!test
%! try
%!   ## TEST 1
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     num_modes = int32(10);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       ri = 8e-3;
%!       ro = 10e-3;
%!       h = 20e-3;
%!       c = 6e-3;
%!       b = h - 2 * c;
%!       scale_def = 5e-3;
%!       mesh_size = 3e-3;
%!       num_modes = 10;
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "ri = %g;\n", ri);
%!       fprintf(fd, "ro = %g;\n", ro);
%!       fprintf(fd, "h = %g;\n", h);
%!       fprintf(fd, "c = %g;\n", c);
%!       fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!       fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!       fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!       fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!       fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!       fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!       fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!       fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,5};\n");
%!       fputs(fd, "Line(5) = {5,6};\n");
%!       fputs(fd, "Line(6) = {6,7};\n");
%!       fputs(fd, "Line(7) = {7,8};\n");
%!       fputs(fd, "Line(8) = {8,1};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!       fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%!       fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%!       fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!     cms_opt.solver = "pardiso";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.name = "node_id_interface1";
%!     cms_opt.element.name = "elem_id_modal";
%!     cms_opt.nodes.modal.name = "node_id_modal";
%!     cms_opt.refine_max_iter = int32(10);
%!     grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%!     grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%!     grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%!     bearing_surf(1).group_idx = grp_id_p1;
%!     bearing_surf(1).name = "journal-surface";
%!     bearing_surf(1).group_id_interface = int32(30);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "journal";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(1).r = ri;
%!     bearing_surf(1).w = b;
%!     bearing_surf(1).X0 = [0; 0; b/2 + c];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(2).group_idx = grp_id_p2;
%!     bearing_surf(2).name = "shell-surface";
%!     bearing_surf(2).group_id_interface = int32(20);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = ro;
%!     bearing_surf(2).w = b;
%!     bearing_surf(2).X0 = [0; 0; b/2 + c];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!     for i=1:numel(bearing_surf)
%!       mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!     endfor
%!     for i=1:numel(bearing_surf)
%!       mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no);
%!     endfor
%!     cms_opt.modes.number = num_modes;
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);
%!     for i=1:numel(comp_mat)
%!       fem_ehd_pre_comp_mat_export(comp_mat(i), bearing_surf(i).options, sprintf("%s_ehd_%d.dat", filename, i));
%!     endfor
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e6);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 2
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 3e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(60);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       scale_def = 5e-3;
%!       mesh_size = 2e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.load_cases = "index";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.solver = "pardiso";
%!     cms_opt.verbose = int32(1);
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + 100;
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.number_of_nodes_x = 60;
%!     bearing_surf(1).options.number_of_nodes_z = 20;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + 100;
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.number_of_nodes_x = 60;
%!     bearing_surf(2).options.number_of_nodes_z = 20;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = true;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = 1;
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end - nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(10);
%!     load_case_post_dof.locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_post(i).loads = zeros(1, 6);
%!       load_case_post(i).loads(i) = 1;
%!     endfor
%!     x1 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     Phi1 = atan2(y1, x1);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(7).pressure.tria6.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%!     load_case_post(8).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(8).pressure.tria6.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!     load_case_post(9).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(9).pressure.tria6.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6);
%!     load_case_post(10).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(10).pressure.tria6.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post_dof);
%!     [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     idx_node = [[mesh.groups.tet10].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 3
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d1 = 8e-3;
%!       D1 = 12e-3;
%!       d2 = 14e-3;
%!       D2 = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       scale_def = 5e-3;
%!       mesh_size = 1.5e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d1 = %g;\n", d1);
%!       fprintf(fd, "D1 = %g;\n", D1);
%!       fprintf(fd, "d2 = %g;\n", d2);
%!       fprintf(fd, "D2 = %g;\n", D2);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {         0.0,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {         0.0,  0.5 * d2, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d2,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {         0.0, -0.5 * d2, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d2,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D2/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D2,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D2/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + int32(100);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(1).r = 0.5 * d1;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d1;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + int32(100);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = 0.5 * d2;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * 0.5 * d2;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!     cms_opt.solver = "umfpack";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     num_modes = int32(30);
%!     num_modes_cms = int32(10);
%!     k1 = 1;
%!     k2 = 0;
%!     err_red = zeros(7, 1);
%!     cms_opt.modes.number = num_modes_cms;
%!     for i=1:numel(bearing_surf)
%!       bearing_surf(i).options.number_of_modes = min(num_modes, floor(numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes) * 3 / 2 - 6));
%!     endfor
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!     [lambda, idx_lambda] = sort(diag(lambda));
%!     qred = qred(:, idx_lambda);
%!     sol_eig_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     for i=1:size(sol_eig_red.def, 3)
%!       sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!     endfor
%!     sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!       switch (i)
%!         case {4, 5}
%!           load_case_itf(i).loads(i) *= w;
%!         case 6
%!           load_case_itf(i).loads(i) *= d1;
%!       endswitch
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     X1 = mesh.nodes(mesh.groups.tria6(bearing_surf(1).group_idx).nodes, 1:3) - bearing_surf(1).X0.';
%!     Phi1 = atan2(X1(:, 2), X1(:, 1));
%!     p1 = k1 * sin(Phi1).^2 * bearing_surf(1).options.reference_pressure;
%!     Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!     z1g = comp_mat(1).bearing_surf.grid_z(:);
%!     p1red = zeros(numel(z1g), numel(Phi1g));
%!     for i=1:columns(p1red)
%!       p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!     endfor
%!     p1red = p1red(:);
%!     X2 = mesh.nodes(mesh.groups.tria6(bearing_surf(2).group_idx).nodes, 1:3) - bearing_surf(2).X0.';
%!     Phi2 = atan2(X2(:, 2), X2(:, 1));
%!     p2 = k2 * cos(Phi2) .* sin(pi * X2(:, 3) / (max(X2(:, 3)) - min(X2(:, 3)))) * bearing_surf(2).options.reference_pressure;
%!     Phi2g = comp_mat(2).bearing_surf.grid_x(:) / (0.5 * comp_mat(2).bearing_dimensions.bearing_diameter);
%!     z2g = comp_mat(2).bearing_surf.grid_z(:);
%!     p2red = zeros(numel(z2g), numel(Phi2g));
%!     for i=1:columns(p2red)
%!       p2red(:, i) = griddata([Phi2 - 2 * pi; Phi2; Phi2 + 2 * pi], repmat(X2(:, 3), 3, 1), repmat(p2, 3, 1), repmat(Phi2g(i), rows(z2g), 1), z2g);
%!     endfor
%!     p2red = p2red(:);
%!     Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf(1).options.reference_pressure + comp_mat(2).E(:, 1:end - nz2) * p2red(1:end -  nz2) / bearing_surf(2).options.reference_pressure;
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     w1red = comp_mat(1).D * qred;
%!     w2red = comp_mat(2).D * qred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(7);
%!     load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!       load_case_post(i).loads = load_case_itf(i).loads;
%!     endfor
%!     load_case_post(7).pressure.tria6.elements = [mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!                                                  mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :)];
%!     pn = zeros(rows(mesh.nodes), 1);
%!     pn(mesh_post.groups.tria6(grp_idx_p1).nodes) = p1;
%!     pn(mesh_post.groups.tria6(grp_idx_p2).nodes) = p2;
%!     load_case_post(7).pressure.tria6.p = [pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :));
%!                                           pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :))];
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!     [mat_ass_post.M, ...
%!      mat_ass_post.K, ...
%!      mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                       dof_map_post, ...
%!                                       [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     node_idx1 = mesh.groups.tria6(bearing_surf(1).group_idx).nodes;
%!     w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!     for i=1:numel(node_idx1)
%!       ni = [mesh.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!       ni /= norm(ni);
%!       w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!     endfor
%!     w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!     for i=1:columns(w1postint)
%!       for m=1:numel(Phi1g)
%!         w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                        repmat(X1(:, 3), 3, 1), ...
%!                                                                        repmat(w1post(:, i), 3, 1), ...
%!                                                                        repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                        z1g);
%!       endfor
%!     endfor
%!     if (do_plot)
%!       for m=1:columns(w1red)
%!         for i=1:numel(z1g)
%!           figure("visible", "off");
%!           hold("on");
%!           plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;r");
%!           plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;k");
%!           ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!           title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!           xlabel("Phi [deg]");
%!           ylabel("w [um]");
%!         endfor
%!       endfor
%!     endif
%!     sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!     for i=1:size(sol_eig_post.def, 3)
%!       sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!     endfor
%!     num_modes_comp = min(10, floor(numel(sol_eig_red.f) / 3));
%!     err_mod = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!     err_w = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     idx_node = [[mesh.groups.tet10].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, 1:3, i) - sol_red.def(idx_node, 1:3, i)))) / max(max(abs(sol_post.def(idx_node, 1:3, i))));
%!     endfor
%!     fprintf(stderr, "flexible interfaces using %d static pressure modes, %d dynamic modes:\n", num_modes, num_modes_cms);
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod);
%!     fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w);
%!     tol_red = 5e-2;
%!     tol_mod = 2e-2;
%!     tol_w = 1e-2;
%!     assert_simple(all(err_red < tol_red));
%!     assert_simple(all(err_mod < tol_mod));
%!     assert_simple(all(err_w < tol_w));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 4
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d1 = 8e-3;
%!       D1 = 12e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_clamp = 3;
%!       scale_def = 5e-3;
%!       mesh_size = 2.5e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d1 = %g;\n", d1);
%!       fprintf(fd, "D1 = %g;\n", D1);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {           0,  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(8) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {          0, 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 1, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 1, 9};\n");
%!       fputs(fd, "Line(7) = {9, 10};\n");
%!       fputs(fd, "Line(8) = {10, 6};\n");
%!       fputs(fd, "Line(9) = {6, 7};\n");
%!       fputs(fd, "Curve Loop(10) = {5, 6, 7, 8, 9};\n");
%!       fputs(fd, "Curve Loop(11) = {1, 2, 3, 4};\n");
%!       fputs(fd, "Plane Surface(12) = {10, 11};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{12}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[7],tmp[8],tmp[9],tmp[10]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"clamp\", %d) = {tmp[5]};\n", grp_id_clamp);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + int32(100);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(1).r = 0.5 * d1;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d1;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_clamp, cms_opt.nodes.modal.number);
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!     cms_opt.inveriants = true;
%!     cms_opt.static_modes = false;
%!     cms_opt.modal_node_constraint = false;
%!     cms_opt.load_cases = "index";
%!     cms_opt.solver = "umfpack";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_modes.shift_A = 0;
%!     opt_modes.refine_max_iter = int32(10);
%!     opt_modes.verbose = int32(0);
%!     opt_modes.solver = cms_opt.solver;
%!     num_modes = int32(15);
%!     num_modes_cms = int32(10);
%!     err_red = zeros(7, 1);
%!     k1 = 1;
%!     cms_opt.modes.number = num_modes_cms;
%!     for i=1:numel(bearing_surf)
%!       bearing_surf(i).options.number_of_modes = min(num_modes, floor(numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes) * 3 / 2));
%!     endfor
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!     [lambda, idx_lambda] = sort(diag(lambda));
%!     qred = qred(:, idx_lambda);
%!     sol_eig_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     for i=1:size(sol_eig_red.def, 3)
%!       sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!     endfor
%!     sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!       switch (i)
%!         case {4, 5}
%!           load_case_itf(i).loads(i) *= w;
%!         case 6
%!           load_case_itf(i).loads(i) *= d1;
%!       endswitch
%!     endfor
%!     [Kitf, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     X1 = mesh.nodes(mesh.groups.tria6(bearing_surf(1).group_idx).nodes, 1:3) - bearing_surf(1).X0.';
%!     Phi1 = atan2(X1(:, 2), X1(:, 1));
%!     p1 = k1 * sin(Phi1).^2 * bearing_surf(1).options.reference_pressure;
%!     Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!     z1g = comp_mat(1).bearing_surf.grid_z(:);
%!     p1red = zeros(numel(z1g), numel(Phi1g));
%!     for i=1:columns(p1red)
%!       p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!     endfor
%!     p1red = p1red(:);
%!     Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf(1).options.reference_pressure;
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     w1red = comp_mat(1).D * qred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(7);
%!     load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!       load_case_post(i).loads = load_case_itf(i).loads;
%!     endfor
%!     load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     pn = zeros(rows(mesh.nodes), 1);
%!     pn(mesh_post.groups.tria6(grp_idx_p1).nodes) = p1;
%!     load_case_post(7).pressure.tria6.p = pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :));
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!     [mat_ass_post.M, ...
%!      mat_ass_post.K, ...
%!      mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                       dof_map_post, ...
%!                                       [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     node_idx1 = mesh.groups.tria6(bearing_surf(1).group_idx).nodes;
%!     w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!     for i=1:numel(node_idx1)
%!       ni = [mesh.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!       ni /= norm(ni);
%!       w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!     endfor
%!     w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!     for i=1:columns(w1postint)
%!       for m=1:numel(Phi1g)
%!         w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                        repmat(X1(:, 3), 3, 1), ...
%!                                                                        repmat(w1post(:, i), 3, 1), ...
%!                                                                        repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                        z1g);
%!       endfor
%!     endfor
%!     if (do_plot)
%!       for m=1:columns(w1red)
%!         for i=1:numel(z1g)
%!           figure("visible", "off");
%!           hold("on");
%!           plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;r");
%!           plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;k");
%!           ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!           title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!           xlabel("Phi [deg]");
%!           ylabel("w [um]");
%!         endfor
%!       endfor
%!     endif
%!     sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!     for i=1:size(sol_eig_post.def, 3)
%!       sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!     endfor
%!     num_modes_comp = min(40, floor(numel(sol_eig_red.f) / 3));
%!     err_mod = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!     err_w = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     idx_node = [[mesh.groups.tet10].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     fprintf(stderr, "flexible interfaces using %d static pressure modes, %d dynamic modes:\n", num_modes, num_modes_cms);
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod);
%!     fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w);
%!     tol_red = 3e-2;
%!     tol_mod = 2e-2;
%!     tol_w = 2e-2;
%!     assert_simple(all(err_red < tol_red));
%!     assert_simple(all(err_mod < tol_mod));
%!     assert_simple(all(err_w < tol_w));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e13);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 5
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d1 = 8e-3;
%!       D1 = 12e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_clamp = 3;
%!       grp_id_itf2 = 4;
%!       scale_def = 5e-3;
%!       mesh_size = 2e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d1 = %g;\n", d1);
%!       fprintf(fd, "D1 = %g;\n", D1);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {           l,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {           l,  0.5 * d1, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d1,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {           l, -0.5 * d1, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d1,       0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {           0,  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {l - Sqrt((D1/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(8) = {l + 0.5 * D1,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9) = {l - Sqrt((D1/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {2 * h, 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {2 * h, 4.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = { h, 4.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = { h, 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {          0, 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 1, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 1, 9};\n");
%!       fputs(fd, "Line(7) = {9, 10};\n");
%!       fputs(fd, "Line(8) = {10, 11};\n");
%!       fputs(fd, "Line(9) = {11, 12};\n");
%!       fputs(fd, "Line(10) = {12, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Line(12) = {14, 6};\n");
%!       fputs(fd, "Line(13) = {6, 7};\n");
%!       fputs(fd, "Curve Loop(14) = {5, 6, 7, 8, 9, 10, 11, 12, 13};\n");
%!       fputs(fd, "Curve Loop(15) = {1, 2, 3, 4};\n");
%!       fputs(fd, "Plane Surface(16) = {14, 15};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{16}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[11],tmp[12],tmp[13],tmp[14]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"clamp\", %d) = {tmp[9]};\n", grp_id_clamp);
%!       fprintf(fd, "Physical Surface(\"itf2\", %d) = {tmp[6]};\n", grp_id_itf2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 3;
%!     cms_opt.nodes.interfaces(1).number = rows(mesh.nodes) + 1;
%!     cms_opt.nodes.interfaces(2).number = rows(mesh.nodes) + 2;
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + int32(100);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).r = 0.5 * d1;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d1;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces(1).number;
%!     bearing_surf(1).options.interpolate_interface = false;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = zeros(1, 3);
%!     mesh.nodes(cms_opt.nodes.interfaces(1).number, 1:3) = bearing_surf(1).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces(2).number, 1:3) = [1.5 * h, 4.5 * h, 0];
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_clamp, cms_opt.nodes.modal.number);
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces(1).number);
%!     mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_itf2, cms_opt.nodes.interfaces(2).number);
%!     cms_opt.solver = "umfpack";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     load_case = fem_pre_load_case_create_empty(6);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     num_modes = int32(10);
%!     num_modes_cms = int32(10);
%!     err_red = zeros(7, 1);
%!     k1 = 1;
%!     cms_opt.modes.number = num_modes_cms;
%!     for i=1:numel(bearing_surf)
%!       bearing_surf(i).options.number_of_modes = min(num_modes, floor(numel(mesh.groups.tria6(bearing_surf(i).group_idx).nodes) * 3 / 2));
%!     endfor
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     [qred, lambda] = eig(mat_ass.Kred, mat_ass.Mred);
%!     [lambda, idx_lambda] = sort(diag(lambda));
%!     qred = qred(:, idx_lambda);
%!     sol_eig_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     for i=1:size(sol_eig_red.def, 3)
%!       sol_eig_red.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_red.def(:, 1:3, i))));
%!     endfor
%!     sol_eig_red.f = imag(sqrt(-lambda)) / (2 * pi);
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces(1).number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = w * d1 * bearing_surf(1).options.reference_pressure;
%!       switch (i)
%!         case {4, 5}
%!           load_case_itf(i).loads(i) *= w;
%!         case 6
%!           load_case_itf(i).loads(i) *= d1;
%!       endswitch
%!     endfor
%!     [Kitf, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     X1 = mesh.nodes(mesh.groups.tria6(bearing_surf(1).group_idx).nodes, 1:3) - bearing_surf(1).X0.';
%!     Phi1 = atan2(X1(:, 2), X1(:, 1));
%!     p1 = k1 * sin(pi/4 + Phi1) * bearing_surf(1).options.reference_pressure;
%!     Phi1g = comp_mat(1).bearing_surf.grid_x(:) / (0.5 * comp_mat(1).bearing_dimensions.bearing_diameter);
%!     z1g = comp_mat(1).bearing_surf.grid_z(:);
%!     p1red = zeros(numel(z1g), numel(Phi1g));
%!     for i=1:columns(p1red)
%!       p1red(:, i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], repmat(X1(:, 3), 3, 1), repmat(p1, 3, 1), repmat(Phi1g(i), rows(z1g), 1), z1g);
%!     endfor
%!     p1red = p1red(:);
%!     Fred = comp_mat(1).E(:, 1:end -  nz1) * p1red(1:end -  nz1) / bearing_surf(1).options.reference_pressure;
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     w1red = comp_mat(1).D * qred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(7);
%!     load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = load_case_itf(i).loaded_nodes;
%!       load_case_post(i).loads = load_case_itf(i).loads;
%!     endfor
%!     load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     pn = zeros(rows(mesh.nodes), 1);
%!     pn(mesh_post.groups.tria6(grp_idx_p1).nodes) = p1;
%!     load_case_post(7).pressure.tria6.p = pn(mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :));
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!     [mat_ass_post.M, ...
%!      mat_ass_post.K, ...
%!      mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                       dof_map_post, ...
%!                                       [FEM_MAT_MASS, FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     node_idx1 = mesh.groups.tria6(bearing_surf(1).group_idx).nodes;
%!     w1post = zeros(numel(node_idx1), size(sol_post.def, 3));
%!     for i=1:numel(node_idx1)
%!       ni = [mesh.nodes(node_idx1(i), 1:2).' - bearing_surf(1).X0(1:2); 0];
%!       ni /= norm(ni);
%!       w1post(i, :) = ni.' * reshape(sol_post.def(node_idx1(i), 1:3, :), 3, size(sol_post.def, 3));
%!     endfor
%!     w1postint = zeros(rows(comp_mat(1).D), columns(w1post));
%!     for i=1:columns(w1postint)
%!       for m=1:numel(Phi1g)
%!         w1postint((m - 1) * numel(z1g) + (1:numel(z1g)), i) = griddata([Phi1 - 2 * pi; Phi1; Phi1 + 2 * pi], ...
%!                                                                        repmat(X1(:, 3), 3, 1), ...
%!                                                                        repmat(w1post(:, i), 3, 1), ...
%!                                                                        repmat(Phi1g(m), numel(z1g), 1), ...
%!                                                                        z1g);
%!       endfor
%!     endfor
%!     if (do_plot)
%!       for m=1:columns(w1red)
%!         for i=1:numel(z1g)
%!           figure("visible", "off");
%!           hold("on");
%!           plot(Phi1g * 180 / pi, 1e6 * w1red(i:numel(z1g):end, m), "-;modal;r");
%!           plot(Phi1g * 180 / pi, 1e6 * w1postint(i:numel(z1g):end, m), "-;nodal;k");
%!           ylim(1e6 * [min(min(w1postint(:, m))), max(max(w1postint(:, m)))]);
%!           title(sprintf("%d modes: i=%d m=%d", num_modes(l), i, m));
%!           xlabel("Phi [deg]");
%!           ylabel("w [um]");
%!         endfor
%!       endfor
%!     endif
%!     sol_eig_post = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, columns(mat_ass.Kred), 0, sqrt(eps), "shift-invert", cms_opt.solver);
%!     for i=1:size(sol_eig_post.def, 3)
%!       sol_eig_post.def(:, :, i) *= 10e-3 / max(max(abs(sol_eig_post.def(:, 1:3, i))));
%!     endfor
%!     num_modes_comp = min(40, floor(numel(sol_eig_red.f) / 3));
%!     err_mod = max(abs(sol_eig_red.f(num_modes_comp) - sol_eig_post.f(num_modes_comp)), [], 2) / max(sol_eig_post.f(num_modes_comp), [], 2);
%!     err_w = max(max(abs(w1red - w1postint))) / max(max(abs(w1postint)));
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     mesh_data(2).mesh.nodes(:, 2) += 40e-3;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     idx_node = [[mesh.groups.tet10].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = norm(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)) / norm(sol_post.def(idx_node, :, i));
%!     endfor
%!     fprintf(stderr, "flexible interfaces using %d static pressure modes, %d dynamic modes:\n",  num_modes, num_modes_cms);
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     fprintf(stderr, "natural frequency: %.1f%%\n", 100 * err_mod);
%!     fprintf(stderr, "bearing radial deformation: %.1f%%\n", 100 * err_w);
%!     tol_red = 1e-2;
%!     tol_mod = 1e-2;
%!     tol_w = 1e-2;
%!     assert_simple(all(err_red < tol_red));
%!     assert_simple(all(err_mod < tol_mod));
%!     assert_simple(all(err_w < tol_w));
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 6
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 4e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(30);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       scale_def = 5e-3;
%!       mesh_size = 1.5e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     cms_opt.algorithm = "shift-invert";
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + int32(100);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.number_of_nodes_x = 60;
%!     bearing_surf(1).options.number_of_nodes_z = 20;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + int32(100);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.number_of_nodes_x = 60;
%!     bearing_surf(2).options.number_of_nodes_z = 20;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number);
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number);
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_modes.refine_max_iter = int32(10);
%!     opt_modes.solver = "pardiso";
%!     opt_modes.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.solver = opt_modes.solver;
%!     cms_opt.number_of_threads = opt_modes.number_of_threads;
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     [qred, lambda_red] = eig(mat_ass.Kred, mat_ass.Mred);
%!     [lambda_red, idx_lambda_red] = sort(diag(lambda_red));
%!     qred = qred(:, idx_lambda_red);
%!     sol_red_modal.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     sol_red_modal.f = sqrt(lambda_red) / (2 * pi);
%!     for i=1:size(sol_red_modal.def, 3)
%!       sol_red_modal.def(:, :, i) *= 10e-3 / max(max(abs(sol_red_modal.def(:, 1:3, i))));
%!     endfor
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = 1;
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(10);
%!     load_case_post(1).locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_post(i).loads = zeros(1, 6);
%!       load_case_post(i).loads(i) = 1;
%!     endfor
%!     x1 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     Phi1 = atan2(y1, x1);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(7).pressure.tria6.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%!     load_case_post(8).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(8).pressure.tria6.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!     load_case_post(9).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(9).pressure.tria6.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6);
%!     load_case_post(10).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(10).pressure.tria6.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post(1));
%!     dof_map_post.parallel.threads_ass = cms_opt.number_of_threads;
%!     [mat_ass_post.M, ...
%!      mat_ass_post.K, ...
%!      mat_ass_post.R] = fem_ass_matrix(mesh_post, ...
%!                                       dof_map_post, ...
%!                                       [FEM_MAT_MASS, ...
%!                                        FEM_MAT_STIFFNESS, ...
%!                                        FEM_VEC_LOAD_CONSISTENT], ...
%!                                       load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, opt_modes);
%!     sol_post_modal = fem_sol_modal(mesh_post, dof_map_post, mat_ass_post, cms_opt.modes.number + 6, 0, sqrt(eps), "shift-invert", opt_modes.solver);
%!     for i=1:size(sol_post_modal.def, 3)
%!       sol_post_modal.def(:, :, i) *= 10e-3 / max(max(abs(sol_post_modal.def(:, 1:3, i))));
%!       if (norm(sol_post_modal.def(:, :, i) + sol_red_modal.def(:, :, i)) < norm(sol_post_modal.def(:, :, i) - sol_red_modal.def(:, :, i)))
%!         sol_post_modal.def(:, :, i) *= -1;
%!       endif
%!     endfor
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     sol_comb_modal.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post_modal.def, 3));
%!     sol_comb_modal.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post_modal.def)), :, :) = sol_post_modal.def;
%!     sol_comb_modal.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red_modal.def)), :, :) = sol_red_modal.def(:, :, 1:size(sol_post_modal.def, 3));
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb_modal.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     idx_node = [[mesh.groups.tet10].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     err_mod_freq = (sol_red_modal.f(1:numel(sol_post_modal.f)) ./ sol_post_modal.f - 1);
%!     MAC = zeros(1, numel(sol_post_modal.f));
%!     for i=1:numel(sol_post_modal.f)
%!       ve = sol_post_modal.def(:, :, i)(:);
%!       vo = sol_red_modal.def(:, :, i)(:);
%!       MAC(i) = (ve.' * vo)^2 / ((ve.' * ve) * (vo.' * vo));
%!     endfor
%!     for i=1:numel(sol_post_modal.f)
%!       fprintf(stderr, ...
%!               "mode %d: reduced: %.1fHz full: %.1fHz difference freq %.1f%% MAC %.4f\n", ...
%!               i, ...
%!               sol_red_modal.f(i), ...
%!               sol_post_modal.f(i), ...
%!               100 * (sol_red_modal.f(i) / sol_post_modal.f(i) - 1), ...
%!               MAC(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 7
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 3e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(30);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       scale_def = 5e-3;
%!       mesh_size = 2e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     opt_msh.elem_type = {"tria6h", "tet10h"};
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6h].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria6h].id] == grp_id_p2);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.load_cases = "index";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.solver = "umfpack";
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + 100;
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.number_of_nodes_x = 60;
%!     bearing_surf(1).options.number_of_nodes_z = 20;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-4 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + 100;
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.number_of_nodes_x = 60;
%!     bearing_surf(2).options.number_of_nodes_z = 20;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = true;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-4 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number, "tria6h");
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number, "tria6h");
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet10h = ones(rows(mesh.elements.tet10h), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_comp_mat.elem_type = "tria6h";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, opt_comp_mat);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = 1;
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(10);
%!     load_case_post_dof.locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_post(i).loads = zeros(1, 6);
%!       load_case_post(i).loads(i) = 1;
%!     endfor
%!     x1 = mesh.nodes(:, 1)(mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.tria6h(mesh.groups.tria6h(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     Phi1 = atan2(y1, x1);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(7).pressure.tria6h.elements = mesh_post.elements.tria6h(mesh_post.groups.tria6h(grp_idx_p1).elements, :);
%!     load_case_post(7).pressure.tria6h.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6h(grp_idx_p1).elements), 6);
%!     load_case_post(8).pressure.tria6h.elements = mesh_post.elements.tria6h(mesh_post.groups.tria6h(grp_idx_p1).elements, :);
%!     load_case_post(8).pressure.tria6h.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!     load_case_post(9).pressure.tria6h.elements = mesh_post.elements.tria6h(mesh_post.groups.tria6h(grp_idx_p2).elements, :);
%!     load_case_post(9).pressure.tria6h.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6h(grp_idx_p2).elements), 6);
%!     load_case_post(10).pressure.tria6h.elements = mesh_post.elements.tria6h(mesh_post.groups.tria6h(grp_idx_p2).elements, :);
%!     load_case_post(10).pressure.tria6h.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post_dof);
%!     [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     idx_node = [[mesh.groups.tet10h].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 8
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     num_modes = int32(10);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       ri = 8e-3;
%!       ro = 10e-3;
%!       h = 20e-3;
%!       c = 6e-3;
%!       b = h - 2 * c;
%!       scale_def = 5e-3;
%!       mesh_size = 3e-3;
%!       num_modes = 10;
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "ri = %g;\n", ri);
%!       fprintf(fd, "ro = %g;\n", ro);
%!       fprintf(fd, "h = %g;\n", h);
%!       fprintf(fd, "c = %g;\n", c);
%!       fprintf(fd, "ms = %g;\n", mesh_size);
%!       fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!       fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!       fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!       fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!       fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!       fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!       fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!       fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,5};\n");
%!       fputs(fd, "Line(5) = {5,6};\n");
%!       fputs(fd, "Line(6) = {6,7};\n");
%!       fputs(fd, "Line(7) = {7,8};\n");
%!       fputs(fd, "Line(8) = {8,1};\n");
%!       fputs(fd, "Transfinite Line(1) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(2) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(3) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(4) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(5) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(6) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(7) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(8) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Surface(1) = {1,2,3,8};\n");
%!       fputs(fd, "Transfinite Surface(2) = {3,4,7,8};\n");
%!       fputs(fd, "Transfinite Surface(3) = {4,5,6,7};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{6}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp1[1]};\n");
%!       fputs(fd, "tmp2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{tmp1[0]}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp1[1],tmp2[1]};\n");
%!       fputs(fd, "Physical Surface(\"clamp\",1) = {tmp1[2],tmp2[2]};\n");
%!       fputs(fd, "Physical Surface(\"load1\",2) = {tmp1[4],tmp2[4]};\n");
%!       fputs(fd, "Physical Surface(\"load2\",3) = {tmp1[8],tmp2[8]};\n");
%!       fputs(fd, "Mesh.ElementOrder = 2;\n");
%!       fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!       fputs(fd, "Mesh.Algorithm = 3;\n");
%!       fputs(fd, "Mesh 3;\n");
%!       fputs(fd, "Coherence Mesh;\n");
%!       fputs(fd, "Mesh.Format = 1;\n");
%!       fprintf(fd, "Save \"%s.msh\";\n", filename);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     ## spawn_wait(spawn("gmsh", {[filename, ".geo"]})); return;
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!     cms_opt.solver = "pardiso";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.name = "node_id_interface1";
%!     cms_opt.element.name = "elem_id_modal";
%!     cms_opt.nodes.modal.name = "node_id_modal";
%!     cms_opt.refine_max_iter = int32(10);
%!     grp_id_clamp = find([[mesh.groups.quad8].id] == 1);
%!     grp_id_p1 = find([[mesh.groups.quad8].id] == 3);
%!     grp_id_p2 = find([[mesh.groups.quad8].id] == 2);
%!     bearing_surf(1).group_idx = grp_id_p1;
%!     bearing_surf(1).name = "journal-surface";
%!     bearing_surf(1).group_id_interface = int32(30);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "journal";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(1).r = ri;
%!     bearing_surf(1).w = b;
%!     bearing_surf(1).X0 = [0; 0; b/2 + c];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(2).group_idx = grp_id_p2;
%!     bearing_surf(2).name = "shell-surface";
%!     bearing_surf(2).group_id_interface = int32(20);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = ro;
%!     bearing_surf(2).w = b;
%!     bearing_surf(2).X0 = [0; 0; b/2 + c];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!     for i=1:numel(bearing_surf)
%!       mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!     endfor
%!     for i=1:numel(bearing_surf)
%!       mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no, "quad8");
%!     endfor
%!     cms_opt.modes.number = num_modes;
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.iso20 = ones(rows(mesh.elements.iso20), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_comp_mat.elem_type = "quad8";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, opt_comp_mat);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);
%!     for i=1:numel(comp_mat)
%!       fem_ehd_pre_comp_mat_export(comp_mat(i), bearing_surf(i).options, sprintf("%s_ehd_%d.dat", filename, i));
%!     endfor
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 9
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     num_modes = int32(10);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       ri = 8e-3;
%!       ro = 10e-3;
%!       h = 20e-3;
%!       c = 6e-3;
%!       b = h - 2 * c;
%!       scale_def = 5e-3;
%!       mesh_size = 3e-3;
%!       num_modes = 10;
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "ri = %g;\n", ri);
%!       fprintf(fd, "ro = %g;\n", ro);
%!       fprintf(fd, "h = %g;\n", h);
%!       fprintf(fd, "c = %g;\n", c);
%!       fprintf(fd, "ms = %g;\n", mesh_size);
%!       fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!       fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!       fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!       fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!       fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!       fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!       fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!       fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,5};\n");
%!       fputs(fd, "Line(5) = {5,6};\n");
%!       fputs(fd, "Line(6) = {6,7};\n");
%!       fputs(fd, "Line(7) = {7,8};\n");
%!       fputs(fd, "Line(8) = {8,1};\n");
%!       fputs(fd, "Transfinite Line(1) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(2) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(3) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(4) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(5) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(6) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(7) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(8) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Surface(1) = {1,2,3,8};\n");
%!       fputs(fd, "Transfinite Surface(2) = {3,4,7,8};\n");
%!       fputs(fd, "Transfinite Surface(3) = {4,5,6,7};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{6}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp1[1]};\n");
%!       fputs(fd, "tmp2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{tmp1[0]}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp1[1],tmp2[1]};\n");
%!       fputs(fd, "Physical Surface(\"clamp\",1) = {tmp1[2],tmp2[2]};\n");
%!       fputs(fd, "Physical Surface(\"load1\",2) = {tmp1[4],tmp2[4]};\n");
%!       fputs(fd, "Physical Surface(\"load2\",3) = {tmp1[8],tmp2[8]};\n");
%!       fputs(fd, "Mesh.ElementOrder = 2;\n");
%!       fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!       fputs(fd, "Mesh.Algorithm = 3;\n");
%!       fputs(fd, "Mesh 3;\n");
%!       fputs(fd, "Coherence Mesh;\n");
%!       fputs(fd, "Mesh.Format = 1;\n");
%!       fprintf(fd, "Save \"%s.msh\";\n", filename);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     ## spawn_wait(spawn("gmsh", {[filename, ".geo"]})); return;
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!     cms_opt.solver = "pardiso";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.name = "node_id_interface1";
%!     cms_opt.element.name = "elem_id_modal";
%!     cms_opt.nodes.modal.name = "node_id_modal";
%!     cms_opt.refine_max_iter = int32(10);
%!     grp_id_clamp = find([[mesh.groups.quad9].id] == 1);
%!     grp_id_p1 = find([[mesh.groups.quad9].id] == 3);
%!     grp_id_p2 = find([[mesh.groups.quad9].id] == 2);
%!     bearing_surf(1).group_idx = grp_id_p1;
%!     bearing_surf(1).name = "journal-surface";
%!     bearing_surf(1).group_id_interface = int32(30);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "journal";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).r = ri;
%!     bearing_surf(1).w = b;
%!     bearing_surf(1).X0 = [0; 0; b/2 + c];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(2).group_idx = grp_id_p2;
%!     bearing_surf(2).name = "shell-surface";
%!     bearing_surf(2).group_id_interface = int32(20);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = ro;
%!     bearing_surf(2).w = b;
%!     bearing_surf(2).X0 = [0; 0; b/2 + c];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!     for i=1:numel(bearing_surf)
%!       mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!     endfor
%!     for i=1:numel(bearing_surf)
%!       mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no, "quad9");
%!     endfor
%!     cms_opt.modes.number = num_modes;
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.iso27 = ones(rows(mesh.elements.iso27), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_comp_mat.elem_type = "quad9";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, opt_comp_mat);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);
%!     for i=1:numel(comp_mat)
%!       fem_ehd_pre_comp_mat_export(comp_mat(i), bearing_surf(i).options, sprintf("%s_ehd_%d.dat", filename, i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 10
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     num_modes = int32(10);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       ri = 8e-3;
%!       ro = 10e-3;
%!       h = 20e-3;
%!       c = 6e-3;
%!       b = h - 2 * c;
%!       scale_def = 5e-3;
%!       mesh_size = 3e-3;
%!       num_modes = 10;
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "ri = %g;\n", ri);
%!       fprintf(fd, "ro = %g;\n", ro);
%!       fprintf(fd, "h = %g;\n", h);
%!       fprintf(fd, "c = %g;\n", c);
%!       fprintf(fd, "ms = %g;\n", mesh_size);
%!       fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!       fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!       fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!       fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!       fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!       fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!       fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!       fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,5};\n");
%!       fputs(fd, "Line(5) = {5,6};\n");
%!       fputs(fd, "Line(6) = {6,7};\n");
%!       fputs(fd, "Line(7) = {7,8};\n");
%!       fputs(fd, "Line(8) = {8,1};\n");
%!       fputs(fd, "Transfinite Line(1) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(2) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(3) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(4) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(5) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(6) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(7) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(8) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Surface(1) = {1,2,3,8};\n");
%!       fputs(fd, "Transfinite Surface(2) = {3,4,7,8};\n");
%!       fputs(fd, "Transfinite Surface(3) = {4,5,6,7};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{6}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp1[1]};\n");
%!       fputs(fd, "tmp2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{tmp1[0]}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp1[1],tmp2[1]};\n");
%!       fputs(fd, "Physical Surface(\"clamp\",1) = {tmp1[2],tmp2[2]};\n");
%!       fputs(fd, "Physical Surface(\"load1\",2) = {tmp1[4],tmp2[4]};\n");
%!       fputs(fd, "Physical Surface(\"load2\",3) = {tmp1[8],tmp2[8]};\n");
%!       fputs(fd, "Mesh.ElementOrder = 1;\n");
%!       fputs(fd, "Mesh.Algorithm = 3;\n");
%!       fputs(fd, "Mesh 3;\n");
%!       fputs(fd, "Coherence Mesh;\n");
%!       fputs(fd, "Mesh.Format = 1;\n");
%!       fprintf(fd, "Save \"%s.msh\";\n", filename);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     ## spawn_wait(spawn("gmsh", {[filename, ".geo"]})); return;
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!     cms_opt.solver = "pardiso";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.name = "node_id_interface1";
%!     cms_opt.element.name = "elem_id_modal";
%!     cms_opt.nodes.modal.name = "node_id_modal";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.verbose = int32(1);
%!     grp_id_clamp = find([[mesh.groups.iso4].id] == 1);
%!     grp_id_p1 = find([[mesh.groups.iso4].id] == 3);
%!     grp_id_p2 = find([[mesh.groups.iso4].id] == 2);
%!     bearing_surf(1).group_idx = grp_id_p1;
%!     bearing_surf(1).name = "journal-surface";
%!     bearing_surf(1).group_id_interface = int32(30);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "journal";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = true;
%!     bearing_surf(1).r = ri;
%!     bearing_surf(1).w = b;
%!     bearing_surf(1).X0 = [0; 0; b/2 + c];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 2 * ro * (1 - cos(mesh_size / (2 * ri)));
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(2).group_idx = grp_id_p2;
%!     bearing_surf(2).name = "shell-surface";
%!     bearing_surf(2).group_id_interface = int32(20);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = true;
%!     bearing_surf(2).r = ro;
%!     bearing_surf(2).w = b;
%!     bearing_surf(2).X0 = [0; 0; b/2 + c];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 2 * ro * (1 - cos(mesh_size / (2 * ri)));
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!     for i=1:numel(bearing_surf)
%!       mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!     endfor
%!     for i=1:numel(bearing_surf)
%!       mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no, "iso4");
%!     endfor
%!     cms_opt.modes.number = num_modes;
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_comp_mat.elem_type = "iso4";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, opt_comp_mat);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);
%!     for i=1:numel(comp_mat)
%!       fem_ehd_pre_comp_mat_export(comp_mat(i), bearing_surf(i).options, sprintf("%s_ehd_%d.dat", filename, i));
%!     endfor
%!     assert(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     for i=1:numel(cond_info)
%!       assert(cond_info(i).D_size(2) == cond_info(i).D_rank);
%!       assert(cond_info(i).D_cond < 1e10);
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 11
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     num_modes = int32(10);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       ri = 8e-3;
%!       ro = 10e-3;
%!       h = 20e-3;
%!       c = 6e-3;
%!       b = h - 2 * c;
%!       scale_def = 5e-3;
%!       mesh_size = 3e-3;
%!       num_modes = 10;
%!       fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "ri = %g;\n", ri);
%!       fprintf(fd, "ro = %g;\n", ro);
%!       fprintf(fd, "h = %g;\n", h);
%!       fprintf(fd, "c = %g;\n", c);
%!       fprintf(fd, "ms = %g;\n", mesh_size);
%!       fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%!       fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%!       fputs(fd, "Point(3) = {ro,0.0,c};\n");
%!       fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%!       fputs(fd, "Point(5) = {ro,0.0,h};\n");
%!       fputs(fd, "Point(6) = {ri,0.0,h};\n");
%!       fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%!       fputs(fd, "Point(8) = {ri,0.0,c};\n");
%!       fputs(fd, "Line(1) = {1,2};\n");
%!       fputs(fd, "Line(2) = {2,3};\n");
%!       fputs(fd, "Line(3) = {3,4};\n");
%!       fputs(fd, "Line(4) = {4,5};\n");
%!       fputs(fd, "Line(5) = {5,6};\n");
%!       fputs(fd, "Line(6) = {6,7};\n");
%!       fputs(fd, "Line(7) = {7,8};\n");
%!       fputs(fd, "Line(8) = {8,1};\n");
%!       fputs(fd, "Transfinite Line(1) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(2) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(3) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(4) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(5) = Round((ro - ri) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(6) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(7) = Round((h - 2 * c) / ms) + 1;\n");
%!       fputs(fd, "Transfinite Line(8) = Round(c / ms) + 1;\n");
%!       fputs(fd, "Transfinite Surface(1) = {1,2,3,8};\n");
%!       fputs(fd, "Transfinite Surface(2) = {3,4,7,8};\n");
%!       fputs(fd, "Transfinite Surface(3) = {4,5,6,7};\n");
%!       fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%!       fputs(fd, "Plane Surface(6) = {5};\n");
%!       fputs(fd, "tmp1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{6}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp1[1]};\n");
%!       fputs(fd, "tmp2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} { Surface{tmp1[0]}; Layers{Round(Pi * ri / ms)}; Recombine;};\n");
%!       fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp2[1]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {tmp1[1],tmp2[1]};\n");
%!       fputs(fd, "Physical Surface(\"clamp\",1) = {tmp1[2],tmp2[2]};\n");
%!       fputs(fd, "Physical Surface(\"load1\",2) = {tmp1[4],tmp2[4]};\n");
%!       fputs(fd, "Physical Surface(\"load2\",3) = {tmp1[8],tmp2[8]};\n");
%!       fputs(fd, "Mesh.ElementOrder = 2;\n");
%!       fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!       fputs(fd, "Mesh.Algorithm = 3;\n");
%!       fputs(fd, "Mesh 3;\n");
%!       fputs(fd, "Coherence Mesh;\n");
%!       fputs(fd, "Mesh.Format = 1;\n");
%!       fprintf(fd, "Save \"%s.msh\";\n", filename);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     ## spawn_wait(spawn("gmsh", {[filename, ".geo"]})); return;
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     opt_msh.elem_type = {"iso20r", "quad8r"};
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 1;
%!     cms_opt.solver = "pardiso";
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.name = "node_id_interface1";
%!     cms_opt.element.name = "elem_id_modal";
%!     cms_opt.nodes.modal.name = "node_id_modal";
%!     cms_opt.refine_max_iter = int32(10);
%!     grp_id_clamp = find([[mesh.groups.quad8r].id] == 1);
%!     grp_id_p1 = find([[mesh.groups.quad8r].id] == 3);
%!     grp_id_p2 = find([[mesh.groups.quad8r].id] == 2);
%!     bearing_surf(1).group_idx = grp_id_p1;
%!     bearing_surf(1).name = "journal-surface";
%!     bearing_surf(1).group_id_interface = int32(30);
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "journal";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).r = ri;
%!     bearing_surf(1).w = b;
%!     bearing_surf(1).X0 = [0; 0; b/2 + c];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(2).group_idx = grp_id_p2;
%!     bearing_surf(2).name = "shell-surface";
%!     bearing_surf(2).group_id_interface = int32(20);
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = ro;
%!     bearing_surf(2).w = b;
%!     bearing_surf(2).X0 = [0; 0; b/2 + c];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * ri;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.modal.number;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.interfaces.number;
%!     for i=1:numel(bearing_surf)
%!       mesh.nodes(bearing_surf(i).master_node_no, 1:3) = bearing_surf(i).X0.';
%!     endfor
%!     for i=1:numel(bearing_surf)
%!       mesh.elements.rbe3(i) = fem_pre_mesh_rbe3_from_surf(mesh, bearing_surf(i).group_idx, bearing_surf(i).master_node_no, "quad8r");
%!     endfor
%!     cms_opt.modes.number = num_modes;
%!     cms_opt.refine_max_iter = int32(10);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_comp_mat.elem_type = "quad8r";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, opt_comp_mat);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     fem_cms_export([filename, "_cms"], mesh, dof_map, mat_ass, cms_opt);
%!     for i=1:numel(comp_mat)
%!       fem_ehd_pre_comp_mat_export(comp_mat(i), bearing_surf(i).options, sprintf("%s_ehd_%d.dat", filename, i));
%!     endfor
%!     assert(cond_info.D_size(2) == cond_info.D_rank);
%!     assert(cond_info.D_cond < 1e10);
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 12
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 3e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(30);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       scale_def = 5e-3;
%!       mesh_size = 4e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "3", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria10].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria10].id] == grp_id_p2);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.load_cases = "index";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.solver = "umfpack";
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + 100;
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + 100;
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = true;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number, "tria10");
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number, "tria10");
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.tet20 = ones(rows(mesh.elements.tet20), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     opt_comp_mat.elem_type = "tria10";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, opt_comp_mat);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = 1;
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(10);
%!     load_case_post_dof.locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_post(i).loads = zeros(1, 6);
%!       load_case_post(i).loads(i) = 1;
%!     endfor
%!     x1 = mesh.nodes(:, 1)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.tria10(mesh.groups.tria10(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     Phi1 = atan2(y1, x1);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(7).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!     load_case_post(7).pressure.tria10.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria10(grp_idx_p1).elements), 10);
%!     load_case_post(8).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p1).elements, :);
%!     load_case_post(8).pressure.tria10.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!     load_case_post(9).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p2).elements, :);
%!     load_case_post(9).pressure.tria10.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria10(grp_idx_p2).elements), 10);
%!     load_case_post(10).pressure.tria10.elements = mesh_post.elements.tria10(mesh_post.groups.tria10(grp_idx_p2).elements, :);
%!     load_case_post(10).pressure.tria10.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post_dof);
%!     [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     idx_node = [[mesh.groups.tet20].nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     assert(cond_info.D_size(2) == cond_info.D_rank);
%!     assert(cond_info.D_cond < 1e10);
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 13
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   function Phi_rb = gen_rigid_body_modes(mesh, dof_map, ref_node_id)
%!     l = mesh.nodes(:, 1:3) - mesh.nodes(ref_node_id, 1:3);
%!     ir = int32([3, 2;
%!                 1, 3;
%!                 2, 1]);
%!     ic = int32([2, 3;
%!                 3, 1;
%!                 1, 2]);
%!     dc = [1, -1];
%!     Arb = repmat(eye(6), rows(mesh.nodes), 1);
%!     for j=1:3
%!       for k=1:2
%!         Arb(ir(j, k):6:end, ic(j, k) + 3) = -dc(k) * l(:, j);
%!       endfor
%!     endfor
%!     Phi_rb = zeros(dof_map.totdof, columns(Arb));
%!     for j=1:columns(dof_map.ndof)
%!       idx_act_dof = find(dof_map.ndof(:, j) > 0);
%!       Phi_rb(dof_map.ndof(idx_act_dof, j), :) = Arb(6 * (idx_act_dof - 1) + j, :);
%!     endfor
%!     Phi_rb = Phi_rb(dof_map.idx_node, :);
%!   endfunction
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 3e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(30);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       scale_def = 5e-3;
%!       mesh_size = 2e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; };\n");
%!       fputs(fd, "Physical Volume(\"volume\", 1) = {tmp[1]};\n");
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.tria6].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.tria6].id] == grp_id_p2);
%!     cms_opt.nodes.interfaces(1).number = rows(mesh.nodes) + 1;
%!     cms_opt.nodes.interfaces(1).include_rigid_body_modes = true;
%!     cms_opt.nodes.interfaces(2).number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces(2).include_rigid_body_modes = false;
%!     cms_opt.floating_frame = true;
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 3;
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.load_cases = "index";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.solver = "umfpack";
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + 100;
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.mesh_size = 1e-3;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces(1).number;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + 100;
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.mesh_size = 1e-3;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-3 * 0.5 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.interfaces(2).number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [l/2, 0, 0];
%!     for i=1:numel(bearing_surf)
%!       mesh.nodes(cms_opt.nodes.interfaces(i).number, 1:3) = bearing_surf(i).X0.';
%!     endfor
%!     mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [grp_id_p1, grp_id_p2], [cms_opt.nodes.interfaces.number]);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     load_case_dof.locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%!     mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     l12 = mesh.nodes(cms_opt.nodes.interfaces(1).number, 1:3) - mesh.nodes(cms_opt.nodes.interfaces(2).number, 1:3);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = [cms_opt.nodes.interfaces.number](:);
%!       load_case_itf(i).loads = zeros(2, 6);
%!       load_case_itf(i).loads(1, i) = 1;
%!       load_case_itf(i).loads(2, :) = -load_case_itf(i).loads(1, :);
%!       load_case_itf(i).loads(2, 4:6) -= cross(l12, load_case_itf(i).loads(1, 1:3));
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(2 * bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end -  nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     Phi_rb = gen_rigid_body_modes(mesh, dof_map, cms_opt.nodes.interfaces(2).number);
%!     Ured = zeros(dof_map.totdof, columns(qred));
%!     Ured(dof_map.idx_node, :) = mat_ass.Tred * qred;
%!     Ured(dof_map.idx_node, :) -= Phi_rb * (Phi_rb \ Ured(dof_map.idx_node, :));
%!     sol_red.def = fem_post_def_nodal(mesh, dof_map, Ured);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(10);
%!     load_case_post_dof.locked_dof = false(size(mesh_post.nodes));
%!     load_case_post_dof.locked_dof(cms_opt.nodes.modal.number, 1:6) = true;
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces(1).number;
%!       load_case_post(i).loads = zeros(1, 6);
%!       load_case_post(i).loads(i) = 1;
%!     endfor
%!     x1 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.tria6(mesh.groups.tria6(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     Phi1 = atan2(y1, x1);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(7).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(7).pressure.tria6.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p1).elements), 6);
%!     load_case_post(8).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p1).elements, :);
%!     load_case_post(8).pressure.tria6.p = sin(2 * Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!     load_case_post(9).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(9).pressure.tria6.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.tria6(grp_idx_p2).elements), 6);
%!     load_case_post(10).pressure.tria6.elements = mesh_post.elements.tria6(mesh_post.groups.tria6(grp_idx_p2).elements, :);
%!     load_case_post(10).pressure.tria6.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.interfaces(2).number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post_dof);
%!     [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!     Phi_rb_post = gen_rigid_body_modes(mesh_post, dof_map_post, cms_opt.nodes.interfaces(2).number);
%!     U_post = fem_sol_linsolve(mat_ass_post.K, mat_ass_post.R, cms_opt);
%!     U_post(dof_map_post.idx_node, :) -= Phi_rb_post * (Phi_rb_post \ U_post(dof_map_post.idx_node, :));
%!     sol_post.def = fem_post_def_nodal(mesh_post, dof_map_post, U_post);
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     opt_merge.group_id = "unique";
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, opt_merge);
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     idx_node = [[mesh.groups.tet10].nodes];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     assert(cond_info.D_size(2) == cond_info.D_rank);
%!     assert(cond_info.D_cond < 1e10);
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! try
%!   ## TEST 14
%!   do_plot = false;
%!   if (do_plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     tol_red = 4e-2;
%!     num_modes_cms = int32(10);
%!     num_modes = int32(60);
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
%!       endif
%!       d = 14e-3;
%!       D = 19.5e-3;
%!       w = 5e-3;
%!       l = 47e-3;
%!       h = 5e-3;
%!       grp_id_volume = 1;
%!       grp_id_p1 = 2;
%!       grp_id_p2 = 3;
%!       p1 = 1;
%!       p2 = 2;
%!       scale_def = 5e-3;
%!       mesh_size = 2e-3;
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "d = %g;\n", d);
%!       fprintf(fd, "D = %g;\n", D);
%!       fprintf(fd, "w = %g;\n", w);
%!       fprintf(fd, "l = %g;\n", l);
%!       fprintf(fd, "h = %g;\n", h);
%!       fprintf(fd, "m = %g;\n", mesh_size);
%!       fputs(fd, "Point(1)  = {          l,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(2)  = {          l,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(3)  = {l + 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(4)  = {          l, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(5)  = {l - 0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(6)  = {        0.0,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(7)  = {        0.0,  0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(8)  = {    0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(9)  = {        0.0, -0.5 * d, -0.5 * w};\n");
%!       fputs(fd, "Point(10) = {   -0.5 * d,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(11) = {l - Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(12) = {l + 0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(13) = {l - Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(14) = {Sqrt((D/2)^2 - (h/2)^2), 0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Point(15) = {   -0.5 * D,      0.0, -0.5 * w};\n");
%!       fputs(fd, "Point(16) = {Sqrt((D/2)^2 - (h/2)^2),  -0.5 * h, -0.5 * w};\n");
%!       fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!       fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!       fputs(fd, "Circle(3) = {4, 1, 5};\n");
%!       fputs(fd, "Circle(4) = {5, 1, 2};\n");
%!       fputs(fd, "Circle(5) = {7, 6, 8};\n");
%!       fputs(fd, "Circle(6) = {8, 6, 9};\n");
%!       fputs(fd, "Circle(7) = {9, 6, 10};\n");
%!       fputs(fd, "Circle(8) = {10, 6, 7};\n");
%!       fputs(fd, "Circle(9) = {11, 1, 12};\n");
%!       fputs(fd, "Circle(10) = {12, 1, 13};\n");
%!       fputs(fd, "Line(11) = {13, 14};\n");
%!       fputs(fd, "Circle(12) = {14, 6, 15};\n");
%!       fputs(fd, "Circle(13) = {15, 6, 16};\n");
%!       fputs(fd, "Line(14) = {16, 11};\n");
%!       fputs(fd, "Curve Loop(15) = {14,13,12,11,10,9};\n");
%!       fputs(fd, "Curve Loop(16) = {4, 3, 2, 1};\n");
%!       fputs(fd, "Curve Loop(17) = {8, 7, 6, 5};\n");
%!       fputs(fd, "Plane Surface(18) = {15, 16, 17};\n");
%!       fputs(fd, "tmp[] = Extrude {0, 0, w} { Surface{18}; Layers{Round(w/m + 1)};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{tmp[0],18};\n");
%!       fprintf(fd, "Physical Volume(\"volume\", %d) = {tmp[1]};\n", grp_id_volume);
%!       fputs(fd, "ReorientMesh Volume{tmp[1]};\n");
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!       fprintf(fd, "Physical Surface(\"small-end\", %d) = {tmp[8],tmp[9],tmp[10],tmp[11]};\n", grp_id_p1);
%!       fprintf(fd, "Physical Surface(\"big-end\", %d) = {tmp[12],tmp[13],tmp[14],tmp[15]};\n", grp_id_p2);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!         fd = -1;
%!       endif
%!     end_unwind_protect
%!     fprintf(stderr, "meshing ...\n");
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "1", "-ho_min", "0.5", "-ho_max", "1.5",  "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       error("gmsh failed with status %d", status);
%!     endif
%!     fprintf(stderr, "loading mesh \"%s\" ...\n", [filename, ".msh"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%!     grp_idx_p1 = find([[mesh.groups.iso4].id] == grp_id_p1);
%!     grp_idx_p2 = find([[mesh.groups.iso4].id] == grp_id_p2);
%!     grp_idx_volume = find([[mesh.groups.iso8.id] == grp_id_volume]);
%!     cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!     cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!     cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!     cms_opt.modes.number = num_modes_cms;
%!     cms_opt.load_cases = "index";
%!     cms_opt.refine_max_iter = int32(10);
%!     cms_opt.solver = "umfpack";
%!     cms_opt.verbose = int32(1);
%!     cms_opt.max_cond_D = 1e5;
%!     cms_opt.lambda_threshold = 1e-3;
%!     cms_opt.tol_gamma_rel = 1e-1;
%!     bearing_surf(1).group_idx = grp_idx_p1;
%!     bearing_surf(1).group_id_interface = grp_id_p1 + 100;
%!     bearing_surf(1).material_id_interface = int32(2);
%!     bearing_surf(1).name = "p1";
%!     bearing_surf(1).options.reference_pressure = 1e9;
%!     bearing_surf(1).options.number_of_nodes_x = 60;
%!     bearing_surf(1).options.number_of_nodes_z = 20;
%!     bearing_surf(1).options.bearing_type = "shell";
%!     bearing_surf(1).options.matrix_type = "modal substruct total";
%!     bearing_surf(1).options.interpolate_interface = false;
%!     bearing_surf(1).r = 0.5 * d;
%!     bearing_surf(1).w = w;
%!     bearing_surf(1).X0 = [l; 0; 0];
%!     bearing_surf(1).R = eye(3);
%!     bearing_surf(1).relative_tolerance = 0;
%!     bearing_surf(1).absolute_tolerance = 1e-2 * 0.5 * d;
%!     bearing_surf(1).options.number_of_modes = num_modes;
%!     bearing_surf(1).master_node_no = cms_opt.nodes.interfaces.number;
%!     bearing_surf(2).group_idx = grp_idx_p2;
%!     bearing_surf(2).group_id_interface = grp_id_p2 + 100;
%!     bearing_surf(2).material_id_interface = int32(2);
%!     bearing_surf(2).name = "p2";
%!     bearing_surf(2).options.reference_pressure = 1e9;
%!     bearing_surf(2).options.number_of_nodes_x = 60;
%!     bearing_surf(2).options.number_of_nodes_z = 20;
%!     bearing_surf(2).options.bearing_type = "shell";
%!     bearing_surf(2).options.matrix_type = "modal substruct total";
%!     bearing_surf(2).options.interpolate_interface = false;
%!     bearing_surf(2).r = 0.5 * d;
%!     bearing_surf(2).w = w;
%!     bearing_surf(2).X0 = [0; 0; 0];
%!     bearing_surf(2).R = eye(3);
%!     bearing_surf(2).relative_tolerance = 0;
%!     bearing_surf(2).absolute_tolerance = 1e-2 * 0.5 * d;
%!     bearing_surf(2).options.number_of_modes = num_modes;
%!     bearing_surf(2).master_node_no = cms_opt.nodes.modal.number;
%!     mesh.nodes(cms_opt.nodes.modal.number, 1:3) = bearing_surf(2).X0.';
%!     mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = bearing_surf(1).X0.';
%!     mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p2, cms_opt.nodes.modal.number, "iso4");
%!     mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, grp_id_p1, cms_opt.nodes.interfaces.number, "iso4");
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     mesh.materials.iso8 = ones(rows(mesh.elements.iso8), 1, "int32");
%!     mesh.material_data(1).E = 210000e6;
%!     mesh.material_data(1).nu = 0.3;
%!     mesh.material_data(1).rho = 7850;
%!     mesh.material_data(2).E = 0.01e6;
%!     mesh.material_data(2).nu = 0.3;
%!     mesh.material_data(2).rho = 0;
%!     options_pre.elem_type = "iso4";
%!     [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig_cms, cond_info] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case_dof, cms_opt, bearing_surf, options_pre);
%!     fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info);
%!     assert_simple(rank(mat_ass.Kred), columns(mat_ass.Kred));
%!     assert_simple(rank(mat_ass.Mred), columns(mat_ass.Mred));
%!     load_case_itf = fem_pre_load_case_create_empty(6);
%!     for i=1:6
%!       load_case_itf(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_itf(i).loads = zeros(1, 6);
%!       load_case_itf(i).loads(i) = 1;
%!     endfor
%!     [~, Ritf] = fem_ass_matrix(mesh, dof_map, [FEM_MAT_STIFFNESS_SYM, FEM_VEC_LOAD_CONSISTENT], load_case_itf);
%!     nx1 = numel(comp_mat(1).bearing_surf.grid_x);
%!     nz1 = numel(comp_mat(1).bearing_surf.grid_z);
%!     nx2 = numel(comp_mat(2).bearing_surf.grid_x);
%!     nz2 = numel(comp_mat(2).bearing_surf.grid_z);
%!     p1red1 = repmat(p1 * bearing_surf(1).options.reference_pressure, (nx1 - 1) * nz1, 1);
%!     p2red1 = repmat(p2 * bearing_surf(2).options.reference_pressure, (nx2 - 1) * nz2, 1);
%!     p1red2 = zeros((nx1 - 1) * nz1, 1);
%!     p2red2 = zeros((nx2 - 1) * nz2, 1);
%!     for i=1:nx1 - 1
%!       p1red2((i - 1) * nz1 + 1:i * nz1) = p1 * sin(bearing_surf(1).grid_x(i) / bearing_surf(1).r) * bearing_surf(1).options.reference_pressure;
%!     endfor
%!     for i=1:nx2 - 1
%!       p2red2((i - 1) * nz2 + 1:i * nz2) = p2 * cos(bearing_surf(2).grid_x(i) / bearing_surf(2).r)^2 * bearing_surf(2).options.reference_pressure;
%!     endfor
%!     Fred = [comp_mat(1).E(:, 1:end - nz1) * [p1red1, p1red2] / bearing_surf(1).options.reference_pressure, ...
%!             comp_mat(2).E(:, 1:end - nz2) * [p2red1, p2red2] / bearing_surf(2).options.reference_pressure];
%!     Fred = [full(mat_ass.Tred.' * Ritf(dof_map.idx_node, :)), Fred];
%!     qred = mat_ass.Kred \ Fred;
%!     sol_red.def = fem_post_cms_expand_body(mesh, dof_map, mat_ass, qred);
%!     mesh_post = mesh;
%!     mesh_post.elements = rmfield(mesh_post.elements, "joints");
%!     load_case_post = fem_pre_load_case_create_empty(10);
%!     load_case_post_dof.locked_dof = false(size(mesh_post.nodes));
%!     for i=1:6
%!       load_case_post(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!       load_case_post(i).loads = zeros(1, 6);
%!       load_case_post(i).loads(i) = 1;
%!     endfor
%!     x1 = mesh.nodes(:, 1)(mesh.elements.iso4(mesh.groups.iso4(grp_idx_p1).elements, :)) - bearing_surf(1).X0(1);
%!     y1 = mesh.nodes(:, 2)(mesh.elements.iso4(mesh.groups.iso4(grp_idx_p1).elements, :)) - bearing_surf(1).X0(2);
%!     x2 = mesh.nodes(:, 1)(mesh.elements.iso4(mesh.groups.iso4(grp_idx_p2).elements, :)) - bearing_surf(2).X0(1);
%!     y2 = mesh.nodes(:, 2)(mesh.elements.iso4(mesh.groups.iso4(grp_idx_p2).elements, :)) - bearing_surf(2).X0(2);
%!     Phi1 = atan2(y1, x1);
%!     Phi2 = atan2(y2, x2);
%!     load_case_post(7).pressure.iso4.elements = mesh_post.elements.iso4(mesh_post.groups.iso4(grp_idx_p1).elements, :);
%!     load_case_post(7).pressure.iso4.p = repmat(p1 * bearing_surf(1).options.reference_pressure, numel(mesh_post.groups.iso4(grp_idx_p1).elements), 4);
%!     load_case_post(8).pressure.iso4.elements = mesh_post.elements.iso4(mesh_post.groups.iso4(grp_idx_p1).elements, :);
%!     load_case_post(8).pressure.iso4.p = sin(Phi1) * p1 * bearing_surf(1).options.reference_pressure;
%!     load_case_post(9).pressure.iso4.elements = mesh_post.elements.iso4(mesh_post.groups.iso4(grp_idx_p2).elements, :);
%!     load_case_post(9).pressure.iso4.p = repmat(p2 * bearing_surf(2).options.reference_pressure, numel(mesh_post.groups.iso4(grp_idx_p2).elements), 4);
%!     load_case_post(10).pressure.iso4.elements = mesh_post.elements.iso4(mesh_post.groups.iso4(grp_idx_p2).elements, :);
%!     load_case_post(10).pressure.iso4.p = cos(Phi2).^2 * p2 * bearing_surf(2).options.reference_pressure;
%!     mesh_post.elements.joints.C = eye(6);
%!     mesh_post.elements.joints.nodes = cms_opt.nodes.modal.number;
%!     dof_map_post = fem_ass_dof_map(mesh_post, load_case_post_dof);
%!     [mat_ass_post.K, mat_ass_post.R] = fem_ass_matrix(mesh_post, dof_map_post, [FEM_MAT_STIFFNESS, FEM_VEC_LOAD_CONSISTENT], load_case_post);
%!     sol_post = fem_sol_static(mesh_post, dof_map_post, mat_ass_post, cms_opt);
%!     mesh_data(1).mesh = mesh;
%!     mesh_data(1).dof_map = dof_map;
%!     mesh_data(1).mesh.nodes(:, 2) += 25e-3;
%!     mesh_data(2).mesh = mesh_post;
%!     mesh_data(2).dof_map = dof_map_post;
%!     [mesh_comb, dof_map_comb] = fem_post_mesh_merge(mesh_data, struct());
%!     sol_comb.def = zeros(rows(mesh_comb.nodes), columns(mesh_comb.nodes), size(sol_post.def, 3));
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(2) + (1:rows(sol_post.def)), :, :) = sol_post.def;
%!     sol_comb.def(dof_map_comb.submesh.offset.nodes(1) + (1:rows(sol_red.def)), :, :) = sol_red.def;
%!     for i=1:size(sol_comb.def, 3)
%!       sol_comb.def(:, :, i) *= 10e-3 / max(max(abs(sol_comb.def(:, 1:3, i))));
%!     endfor
%!     err_red = zeros(1, size(sol_post.def, 3));
%!     idx_node = [mesh.groups.iso8(grp_idx_volume).nodes, cms_opt.nodes.modal.number, cms_opt.nodes.interfaces.number];
%!     for i=1:size(sol_post.def, 3)
%!       err_red(i) = max(max(abs(sol_post.def(idx_node, :, i) - sol_red.def(idx_node, :, i)))) / max(max(abs(sol_post.def(idx_node, :, i))));
%!     endfor
%!     for i=1:size(sol_post.def, 3)
%!       fprintf(stderr, "mode %d: %.1f%%\n", i, 100 * err_red(i));
%!     endfor
%!     assert_simple(all(err_red < tol_red));
%!     assert(cond_info.D_size(2) == cond_info.D_rank);
%!     assert(cond_info.D_cond < 1e10);
%!     if (do_plot)
%!       figure("visible", "off");
%!       semilogy(cond_info.eta);
%!       xlabel("mode #");
%!       ylabel("eta");
%!       grid minor on;
%!       figure("visible", "off");
%!       semilogy(cond_info.S);
%!       xlabel("mode #");
%!       ylabel("S");
%!       grid on;
%!       figure("visible", "off");
%!       semilogy(cond_info.lambda);
%!       xlabel("mode #");
%!       ylabel("lambda");
%!       for i=1:numel(cond_info.comp_mat)
%!        figure("visible", "off");
%!        semilogy(cond_info.comp_mat(i).cond_value);
%!        xlabel("mode #");
%!        ylabel("cond");
%!        grid minor on;
%!        figure("visible", "off");
%!        semilogy(cond_info.comp_mat(i).gamma);
%!        xlabel("mode #");
%!        ylabel("gamma");
%!       endfor
%!     endif
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         status = unlink(fullfile(fn(i).folder, fn(i).name));
%!         if (status ~= 0)
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
