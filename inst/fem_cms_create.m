## Copyright (C) 2019(-2022) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{mat_ass}, @var{dof_map}, @var{sol_eig}, @var{cms_opt}] = fem_cms_create(@var{mesh}, @var{load_case}, @var{cms_opt})
## Build a reduced order model by using the Craig Bampton approach. That model can be exported to MBDyn by means of fem_cms_export.
##
## @var{mesh} @dots{} Finite element mesh data structure containing constraints for the modal node
##
## @var{load_case} @dots{} Struct array of load cases. The first load case will be used for assembly of global finite element matrices. Additional load cases may be defined, and it's load vectors will be returned at the end of @var{mat_ass}.R.
##
## @var{cms_opt}.nodes.modal @dots{} Struct containing node number and node name of the modal node. Appropriate constraints must be defined for that node in @var{mesh} or @var{load_case}
##
## @var{cms_opt}.nodes.interfaces @dots{} Struct array containing node numbers and node names for interface nodes accessible to MBDyn. Static mode shapes will be generated for those nodes.
##
## @var{cms_opt}.verbose @dots{} Enable verbose output
##
## @var{cms_opt}.algorithm @dots{} Algorithm used for eigenanalysis
##
## @var{cms_opt}.number_of_threads @dots{} Number of threads used for linear solver and assembly
##
## @var{cms_opt}.threshold_elem @dots{} Minimum number of elements for multithreaded assembly
##
## @var{cms_opt}.refine_max_iter @dots{} Maximum number of refinement iterations for the linear solver
##
## @var{cms_opt}.scaling @dots{} Various options for scaling mode shapes
##
## @var{cms_opt}.load_cases @dots{} Define if static mode shapes will be generated for all load cases, or for interface nodes only.
##
## @var{cms_opt}.invariants @dots{} Define if in-variants will be assembled directly, or if the are computed by MBDyn from the diagonal mass matrix.
## @seealso{fem_cms_export}
## @end deftypefn

function [mesh, mat_ass, dof_map, sol_eig, cms_opt] = fem_cms_create(mesh, load_case, cms_opt)
  if (nargin ~= 3 || nargout > 5)
    print_usage();
  endif

  if (~isstruct(load_case))
    error("load_case must be a scalar struct");
  endif

  if (~isfield(cms_opt, "verbose"))
    cms_opt.verbose = int32(0);
  endif

  if (~isfield(cms_opt, "algorithm"))
    cms_opt.algorithm = "shift-invert";
  endif

  if (~isfield(cms_opt, "number_of_threads"))
    cms_opt.number_of_threads = int32(1);
  endif

  if (~isfield(cms_opt, "threshold_elem"))
    cms_opt.threshold_elem = int32(10000);
  endif

  if (~isfield(cms_opt, "refine_max_iter"))
    cms_opt.refine_max_iter = int32(10);
  endif

  if (~isfield(cms_opt, "pre_scaling"))
    cms_opt.pre_scaling = false;
  endif

  if (~isfield(cms_opt, "solver"))
    cms_opt.solver = fem_sol_select(true);
  endif

  if (~isfield(cms_opt, "scaling"))
    cms_opt.scaling = "diag M";
  endif

  if (~isfield(cms_opt, "load_cases"))
    cms_opt.load_cases = "interface";
  endif

  if (~isfield(cms_opt, "invariants"))
    cms_opt.invariants = true;
  endif

  if (~isfield(cms_opt, "modes"))
    cms_opt.modes = struct();
  endif

  if (~isfield(cms_opt.modes, "number"))
    cms_opt.modes.number = int32(0);
  endif

  if (~isfield(cms_opt, "tol"))
    cms_opt.tol = 1e-4;
  endif

  if (~isfield(cms_opt, "static_modes"))
    cms_opt.static_modes = true;
  endif

  if (~isfield(cms_opt, "modal_node_constraint"))
    cms_opt.modal_node_constraint = true;
  endif

  if (~isfield(cms_opt, "floating_frame"))
    cms_opt.floating_frame = false;
  endif

  if (cms_opt.modal_node_constraint)
    node_idx_itf = int32([cms_opt.nodes.modal.number]);
    istart_idx_joint = int32(2);
  else
    node_idx_itf = [];
    istart_idx_joint = int32(1);
  endif

  if (cms_opt.static_modes)
    node_idx_itf = int32([node_idx_itf, cms_opt.nodes.interfaces.number]);
  endif

  dof_in_use = fem_cms_dof_active(mesh);

  if (isfield(mesh.elements, "joints"))
    cms_opt.first_joint_idx_itf = numel(mesh.elements.joints);
  else
    cms_opt.first_joint_idx_itf = int32(0);
  endif

  cms_opt.num_modes_itf = int32(0);

  master_node_joint_idx_curr = cms_opt.first_joint_idx_itf;
  master_node_joint_idx = zeros(1, numel(node_idx_itf), "int32");

  for i=1:numel(node_idx_itf)
    dof_idx_master_node = find(dof_in_use(node_idx_itf(i), :) & ~load_case(1).locked_dof(node_idx_itf(i), :));

    if (~numel(dof_idx_master_node))
      continue;
    endif

    master_node_joint_idx(i) = ++master_node_joint_idx_curr;
    mesh.elements.joints(master_node_joint_idx(i)).nodes = node_idx_itf(i);
    mesh.elements.joints(master_node_joint_idx(i)).C = eye(6)(dof_idx_master_node, :);

    for j=1:numel(load_case)
      load_case(j).joints(master_node_joint_idx(i)).U = zeros(rows(mesh.elements.joints(master_node_joint_idx(i)).C), 1);
    endfor

    if (i > 1)
      cms_opt.num_modes_itf += rows(mesh.elements.joints(master_node_joint_idx(i)).C);
    endif
  endfor

  for i=1:numel(load_case)
    for j=1:numel(load_case(i).joints)
      if (isempty(load_case(i).joints(j).U))
        load_case(i).joints(j).U = zeros(rows(mesh.elements.joints(j).C), 1);
      endif
    endfor
  endfor

  if (cms_opt.num_modes_itf > 0)
    ## Reduce memory consumption: do not duplicate load_case.locked_dof!
    load_case_itf = repmat(setfield(load_case(1), "locked_dof", []), 1, cms_opt.num_modes_itf - 1);

    load_case_cms = fem_pre_load_case_merge(load_case(1), load_case_itf, load_case(2:end));

    clear load_case_itf;
  else
    load_case_cms = load_case;
  endif

  idx_load_case = int32(0);

  for i=istart_idx_joint:numel(node_idx_itf)
    if (master_node_joint_idx(i))
      for j=1:rows(mesh.elements.joints(master_node_joint_idx(i)).C)
        load_case_cms(++idx_load_case).joints(master_node_joint_idx(i)).U(j) = 1;
      endfor
    endif
  endfor

  dof_map = fem_ass_dof_map(mesh, load_case(1));

  dof_map.parallel.threads_ass = cms_opt.number_of_threads;
  dof_map.parallel.threshold_elem = cms_opt.threshold_elem;
  mat_type_stiffness = FEM_MAT_STIFFNESS_SYM_L;
  mat_type_mass = FEM_MAT_MASS_SYM_L;
  mat_type_damping = FEM_MAT_DAMPING_SYM_L;

  cms_opt.solver = fem_sol_select(true, cms_opt.solver);
  cms_opt.symmetric = true;

  switch (cms_opt.solver)
    case {"umfpack", "lu", "mldivide"}
      mat_type_stiffness = FEM_MAT_STIFFNESS;
      cms_opt.symmetric = false;
  endswitch

  switch (cms_opt.algorithm)
    case "eliminate"
      mat_type_stiffness = FEM_MAT_STIFFNESS;
      mat_type_mass = FEM_MAT_MASS;
  endswitch

  if (cms_opt.invariants)
    [mat_ass.M, ...
     mat_ass.K, ...
     mat_ass.D, ...
     mat_ass.R, ...
     mat_ass.dm, ...
     mat_ass.S, ...
     mat_ass.J, ...
     mat_ass.mat_info] = fem_ass_matrix(mesh, ...
                                        dof_map, ...
                                        [mat_type_mass, ...
                                         mat_type_stiffness, ...
                                         mat_type_damping, ...
                                         FEM_VEC_LOAD_CONSISTENT, ...
                                         FEM_SCA_TOT_MASS, ...
                                         FEM_VEC_INERTIA_M1, ...
                                         FEM_MAT_INERTIA_J], ...
                                        load_case_cms);
  else
    [mat_ass.M, ...
     mat_ass.K, ...
     mat_ass.D, ...
     mat_ass.Mlumped, ...
     mat_ass.R, ...
     mat_ass.mat_info] = fem_ass_matrix(mesh, ...
                                        dof_map, ...
                                        [mat_type_mass, ...
                                         mat_type_stiffness, ...
                                         mat_type_damping, ...
                                         FEM_MAT_MASS_LUMPED, ...
                                         FEM_VEC_LOAD_CONSISTENT], ...
                                        load_case_cms);
  endif

  if (cms_opt.verbose)
    fprintf(stderr, "K: %d DOF\n", dof_map.totdof);
    fprintf(stderr, "K: %d nonzeros %.1fMB\n", nnz(mat_ass.K), sizeof(mat_ass.K) / 2^20);
  endif

  switch (cms_opt.load_cases)
    case "interface"
      R_itf = mat_ass.R(:, 1:cms_opt.num_modes_itf);
    case "all"
      R_itf = mat_ass.R;
    case "index"
      R_itf = mat_ass.R(:, cms_opt.load_cases_index);
    otherwise
      error("invalid option cms_opt.load_cases=\"%s\"", cms_opt.load_cases);
  endswitch

  switch (mat_type_stiffness)
    case {FEM_MAT_STIFFNESS_SYM, FEM_MAT_STIFFNESS_SYM_L}
      Ksym = fem_mat_sym(mat_ass.K);
    otherwise
      Ksym = mat_ass.K;
  endswitch

  switch (mat_type_mass)
    case {FEM_MAT_MASS_SYM, FEM_MAT_MASS_SYM_L}
      Msym = fem_mat_sym(mat_ass.M);
    otherwise
      Msym = mat_ass.M;
  endswitch

  Dsym = fem_mat_sym(mat_ass.D);

  mat_ass.Tred = zeros(numel(dof_map.idx_node), columns(R_itf) + cms_opt.modes.number);

  opt_sol.number_of_threads = cms_opt.number_of_threads;
  opt_sol.refine_max_iter = cms_opt.refine_max_iter;
  opt_sol.pre_scaling = cms_opt.pre_scaling;
  opt_sol.solver = cms_opt.solver;
  opt_sol.verbose = cms_opt.verbose;

  if (mat_type_stiffness == FEM_MAT_STIFFNESS_SYM_L || mat_type_mass == FEM_MAT_MASS_SYM_L)
    opt_sol.pre_scaling = false;
  endif

  Kfact = fem_sol_factor(mat_ass.K, opt_sol);

  if (columns(R_itf))
    mat_ass.Tred(:, cms_opt.modes.number + (1:columns(R_itf))) = (Kfact \ R_itf)(dof_map.idx_node, :);
  endif

  if (cms_opt.verbose)
    fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
    whos();
  endif

  if (cms_opt.modes.number > 0)
    if (cms_opt.verbose)
      fprintf(stderr, "solving for normal modes ...\n");
      tic();
    endif

    switch (cms_opt.algorithm)
      case {"unsymmetric", "shift-invert", "diag-shift-invert"}
        opts.disp = 0;
        opts.maxit = 50000;
        opts.tol = 0;

        rndstate = rand("state");

        unwind_protect
          rand("seed", 0);

          switch (cms_opt.algorithm)
            case {"shift-invert", "diag-shift-invert"}
              SIGMA = 0;
              op{1} = @(x) Msym * x;
              op{2} = @(x) Kfact \ x;

              [PHI_d, mu] = eig_sym(op, columns(Msym), cms_opt.modes.number, SIGMA, opts);

              clear op;
            case "unsymmetric"
              SIGMA = "LM";
              iter_eig = int32(0);
              max_iter_eig = int32(100);
              opts.issym = eigs_sym(Kfact);
              opts.isreal = true;

              eigs_callback = @(x) eigs_func(Kfact, Msym, x);

              while (true)
                [PHI_d, mu, status] = eigs(eigs_callback, columns(Msym), cms_opt.modes.number, SIGMA, opts);

                if (status ~= 0)
                  error("eigs failed with status %d", status);
                endif

                if (isreal(PHI_d))
                  break;
                endif

                ## If PHI_d is not real, call eigs again with a new real starting vector

                if (++iter_eig > max_iter_eig)
                  error("eigs returned complex eigenvectors");
                endif

                opts.v0 = real(PHI_d(:, 1));
              endwhile

              PHI_d = eigs_post(Kfact, PHI_d);

              clear opts eigs_callback;
          endswitch
        unwind_protect_cleanup
          rand("state", rndstate);
        end_unwind_protect

        switch (cms_opt.algorithm)
          case "unsymmetric"
            [mu, idx_mu] = sort(1 ./ diag(mu));

            PHI_d = PHI_d(:, idx_mu);

            lambda = sqrt(-mu).';
          case {"shift-invert", "diag-shift-invert"}
            lambda = sqrt(-diag(mu)).';
        endswitch

        PHI_d *= diag(1 ./ max(abs(PHI_d), [], 1)); ## Unit scale results in better condition numbers of Mred and Kred

        err = zeros(columns(PHI_d), 1);

        for i=1:columns(PHI_d)
          v1 = Ksym * PHI_d(:, i);
          v2 = Msym * PHI_d(:, i) * lambda(i)^2;
          err(i) = norm(v1 + v2) / max(norm(v1), norm(v2));
        endfor

        if (cms_opt.verbose)
          for i=1:columns(PHI_d)
            fprintf(stderr, "mode %d: f=%.1f error = %g\n", i, imag(lambda(i)) / (2 * pi), err(i));
          endfor
        endif

        if (any(err > cms_opt.tol))
          warning("eigs algorithm %s failed to converge max(err)=%g", cms_opt.algorithm, max(err));
        endif

        mat_ass.Tred(:, 1:cms_opt.modes.number) = PHI_d(dof_map.idx_node, :);

        if (cms_opt.verbose)
          fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
          whos();
        endif

        clear PHI_d;
      case "eliminate"
        if (cms_opt.verbose)
          fprintf(stderr, "eliminating multipliers of Lagrange ...\n");
        endif

        [Tc, Kc, Mc] = fem_cms_constr_elim(mesh, dof_map, mat_ass);

        if (cms_opt.verbose)
          fprintf(stderr, "solving for normal modes ...\n");
        endif

        [Phi_d, lambda] = fem_sol_eigs(Kc, Mc, cms_opt.modes.number);

        mat_ass.Tred(:, 1:cms_opt.modes.number) = Tc * Phi_d;

        clear Phi_d Tc Kc Mc;
      otherwise
        error("unknown algorithm \"%s\"", cms_opt.algorithm);
    endswitch

    if (cms_opt.verbose)
      toc();
    endif
  else
    Phi_d = zeros(numel(dof_map.idx_node), 0);
    lambda = [];
  endif

  clear Kfact;

  if (cms_opt.verbose)
    fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
    whos();
    fprintf(stderr, "building reduced matrices ...\n");
  endif

  if (cms_opt.floating_frame)
    A = zeros(6 * rows(mesh.nodes), 6);

    for j=1:6
      A(j:6:end, j) = 1;
    endfor

    dX = mesh.nodes(:, 1:3) - mesh.nodes(cms_opt.nodes.modal.number, 1:3);

    A(1:6:end, 5) = dX(:, 3);
    A(2:6:end, 4) = -dX(:, 3);
    A(1:6:end, 6) = -dX(:, 2);
    A(3:6:end, 4) = dX(:, 2);
    A(2:6:end, 6) = dX(:, 1);
    A(3:6:end, 5) = -dX(:, 1);

    idxtrans = zeros(3 * rows(mesh.nodes), 1, "int32");

    for j=1:3
      idxtrans(j:3:end) = j:6:6 * rows(mesh.nodes);
    endfor

    idxmodal = (cms_opt.nodes.modal.number - 1) * 6 + (1:6);

    A(idxmodal, :) = 0;
    Adisp = A(idxtrans, :);

    for j=1:columns(mat_ass.Tred)
      Uj = zeros(6 * rows(mesh.nodes), 1);

      for k=1:6
        ndofk = dof_map.ndof(:, k);
        idxndofk = find(ndofk > 0);
        idxnodek = k:6:rows(Uj);
        Uj(idxnodek(idxndofk)) = mat_ass.Tred(ndofk(idxndofk), j);
      endfor

      Uj -= A * (Adisp \ Uj(idxtrans));

      for k=1:6
        ndofk = dof_map.ndof(:, k);
        idxndofk = find(ndofk > 0);
        idxnodek = k:6:rows(Uj);
        mat_ass.Tred(ndofk(idxndofk), j) = Uj(idxnodek(idxndofk));
      endfor
    endfor
  endif

  mat_ass.Mred = fem_cms_matrix_trans(mat_ass.Tred, Msym(dof_map.idx_node, dof_map.idx_node), "Lower");
  mat_ass.Kred = fem_cms_matrix_trans(mat_ass.Tred, Ksym(dof_map.idx_node, dof_map.idx_node), "Lower");
  mat_ass.Dred = fem_cms_matrix_trans(mat_ass.Tred, Dsym(dof_map.idx_node, dof_map.idx_node), "Lower");

  switch (cms_opt.algorithm)
    case "diag-shift-invert"
      [PHI_diag, lambda_diag] = eig(mat_ass.Kred, mat_ass.Mred, "chol", "vector");
      mat_ass.Mred = PHI_diag.' * mat_ass.Mred * PHI_diag;
      mat_ass.Kred = PHI_diag.' * mat_ass.Kred * PHI_diag;
      mat_ass.Dred = PHI_diag.' * mat_ass.Dred * PHI_diag;
      mat_ass.Tred *= PHI_diag;

      clear PHI_diag lambda_diag;
  endswitch

  switch (cms_opt.scaling)
    case "none"
    otherwise
      switch (cms_opt.scaling)
        case "max K"
          s = max(abs(mat_ass.Kred), [], 2);
        case "max M"
          s = max(abs(mat_ass.Mred), [], 2);
        case "max K,M"
          s = max(max(abs(mat_ass.Kred), [], 2),  max(abs(mat_ass.Mred), [], 2));
        case "norm K"
          s = norm(mat_ass.Kred, "cols");
        case "norm M"
          s = norm(mat_ass.Mred, "cols");
        case "norm K,M"
          s = max(norm(mat_ass.Kred, "cols"), norm(mat_ass.Mred, "cols"));
        case "diag K"
          s = abs(diag(mat_ass.Kred));
        case "diag M"
          s = abs(diag(mat_ass.Mred));
        case "lambda"
          s = abs(diag(mat_ass.Kred)) ./ abs(diag(mat_ass.Mred));
        case "Tred"
          s = max(abs(mat_ass.Tred), [], 1).^2;
        case "mean M,K"
          s = diag(mat_ass.Kred) * mean(abs(diag(mat_ass.Mred)) ./ abs(diag(mat_ass.Kred)));
        case "mean K,M"
          s = diag(mat_ass.Mred) * mean(abs(diag(mat_ass.Kred)) ./ abs(diag(mat_ass.Mred)));
        otherwise
          error("invalid option cms_opt.scaling=\"%s\"", cms_opt.scaling);
      endswitch

      s = 1 ./ sqrt(s);

      if (cms_opt.verbose)
        fprintf(stderr, "condition number before scaling ...\n");
        fprintf(stderr, "cond(Kred)=%.1e\n", cond(mat_ass.Kred));
        fprintf(stderr, "cond(Mred)=%.1e\n", cond(mat_ass.Mred));
      endif

      mat_ass.Kred = diag(s) * mat_ass.Kred * diag(s);
      mat_ass.Mred = diag(s) * mat_ass.Mred * diag(s);
      mat_ass.Dred = diag(s) * mat_ass.Dred * diag(s);
      mat_ass.Tred *= diag(s);

      clear s;

      if (cms_opt.verbose)
        fprintf(stderr, "condition number after scaling ...\n");
      endif
  endswitch

  warning("error", "Octave:singular-matrix", "local");
  warning("error", "Octave:nearly-singular-matrix", "local");

  if (cms_opt.verbose)
    fprintf(stderr, "cond(Kred)=%.1e\n", cond(mat_ass.Kred));
    fprintf(stderr, "cond(Mred)=%.1e\n", cond(mat_ass.Mred));
  endif

  if (~cms_opt.invariants)
    diagMlumped = full(diag(mat_ass.Mlumped));
    mat_ass.diagM = zeros(6 * rows(mesh.nodes), 1);
  endif

  if (~cms_opt.invariants)
    ridx_node(dof_map.idx_node) = 1:numel(dof_map.idx_node);
    for i=1:columns(dof_map.ndof)
      idx_act_dof = find(dof_map.ndof(:, i) > 0);
      idx_glob_dof = dof_map.ndof(idx_act_dof, i);
      mat_ass.diagM((idx_act_dof - 1) * 6 + i) = diagMlumped(idx_glob_dof);
    endfor
  endif

  if (nargout >= 4 || cms_opt.invariants)
    sol_eig.def = zeros(rows(mesh.nodes), columns(mesh.nodes), columns(mat_ass.Tred));
    PHI = zeros(dof_map.totdof, 1);

    for i=1:columns(mat_ass.Tred)
      PHI(dof_map.idx_node) = mat_ass.Tred(:, i);
      sol_eig.def(:, :, i) = fem_post_def_nodal(mesh, dof_map, PHI);
    endfor

    clear PHI;
    sol_eig.f = [imag(lambda) / (2 * pi), repmat(-inf, 1, columns(R_itf))];
  endif

  if (cms_opt.invariants)
    if (cms_opt.verbose)
      fprintf(stderr, "total number of modes: %d\n", size(sol_eig.def, 3));
      fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
      whos();
      fprintf(stderr, "building invariants ...\n");
      tic();
    endif

    [mat_ass.Inv3, ...
     mat_ass.Inv4, ...
     mat_ass.Inv5, ...
     mat_ass.Inv8, ...
     mat_ass.Inv9] = fem_ass_matrix(setfield(mesh, "nodes", mesh.nodes - mesh.nodes(cms_opt.nodes.modal.number, :)), ...
                                    dof_map, ...
                                    [FEM_MAT_INERTIA_INV3, ...
                                     FEM_MAT_INERTIA_INV4, ...
                                     FEM_MAT_INERTIA_INV5, ...
                                     FEM_MAT_INERTIA_INV8, ...
                                     FEM_MAT_INERTIA_INV9], ...
                                    load_case_cms, ...
                                    sol_eig);
    if (cms_opt.verbose)
      toc();
      fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
      whos();
    endif
  endif
endfunction

