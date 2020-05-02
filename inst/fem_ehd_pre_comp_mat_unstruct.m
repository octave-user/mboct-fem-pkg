## Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{comp_mat} = fem_ehd_pre_comp_mat_unstruct(@var{mesh}, @var{mat_ass}, @var{dof_map}, @var{cms_opt}, @var{bearing_surf})
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
## @var{bearing_surf}.R @dots{} Orientation of the bearing surface. R(:, 3) represents the axis of the journal bearing.
##
## @var{bearing_surf}.X0 @dots{} Location of the bearing centre.
##
## @var{bearing_surf}.r @dots{} Radius of the cylindrical bearing surface.
##
## @var{bearing_surf}.w @dots{} Width of the cylindrical bearing surface.
##
## @var{bearing_surf}.nodes @dots{} Node numbers of all Finite Element nodes at the bearing surface.
##
## @var{bearing_surf}.matrix_type @dots{} One of "nodal", "nodal substruct", "nodal substruct total".
##
## @var{bearing_surf}.master_node_no @dots{} Node number of the master node connected to the RBE3 element of the bearing surface.
##
## @var{bearing_surf}.idx_load_case @dots{} Column indices of all load cases in @var{mat_ass}.R related to this bearing surface.
##
## @seealso{fem_ehd_pre_comp_mat_load_case, fem_ehd_pre_comp_mat_export, fem_ehd_pre_comp_mat_plot}
## @end deftypefn

function [comp_mat] = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf)
  if (nargin ~= 5 || nargout > 1)
    print_usage();
  endif

  warning("error", "Octave:singular-matrix", "local");
  warning("error", "Octave:nearly-singular-matrix", "local");

  empty_cell = cell(1, numel(bearing_surf));

  comp_mat = struct("C", empty_cell, ...
                    "condest_C", empty_cell, ...
                    "D", empty_cell, ...
                    "E", empty_cell, ...
                    "xi", empty_cell, ...
                    "zi", empty_cell, ...
                    "bearing_surf", empty_cell, ...
                    "bearing_dimensions", empty_cell, ...
                    "reference_pressure", empty_cell);

  have_PHI = false;
  tol_modal_res = sqrt(eps);
  chunk_size = int32(200);

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
    
    Nxz = [numel(comp_mat(i).bearing_surf.grid_x), numel(comp_mat(i).bearing_surf.grid_z)];

    xi_S = [comp_mat(i).bearing_surf.grid_x(end - 1) - 2 * pi * bearing_surf(i).r, ...
            comp_mat(i).bearing_surf.grid_x, ...
            comp_mat(i).bearing_surf.grid_x(2) + 2 * pi * bearing_surf(i).r];

    zi_S = [2 * comp_mat(i).bearing_surf.grid_z(1) - comp_mat(i).bearing_surf.grid_z(2), ...
            comp_mat(i).bearing_surf.grid_z, ...
            2 * comp_mat(i).bearing_surf.grid_z(end) - comp_mat(i).bearing_surf.grid_z(end - 1)];

    [xx_S, zz_S] = meshgrid(comp_mat(i).bearing_surf.grid_x, comp_mat(i).bearing_surf.grid_z);

    X = bearing_surf(i).R.' * (mesh.nodes(bearing_surf(i).nodes, 1:3).' - bearing_surf(i).X0);

    comp_mat(i).xi = mod(atan2(X(2, :), X(1, :)), 2 * pi) * bearing_surf(i).r;
    comp_mat(i).zi = X(3, :);
    
    xi_U = [comp_mat(i).xi - 2 * pi * bearing_surf(i).r, comp_mat(i).xi, comp_mat(i).xi + 2 * pi * bearing_surf(i).r];
    zi_U = repmat(comp_mat(i).zi, 1, 3);

    ## Triangulate data in advance
    tri_U = delaunay(xi_U, zi_U);
    
    n = s * X(1:2, :);
    r = norm(n, "cols");
    n ./= r;

    switch (bearing_surf(i).options.matrix_type)
      case {"nodal substruct", "nodal substruct total"}
        use_modal_contrib = true;
      case "nodal"
        use_modal_contrib = false;
      otherwise
        error("compliance model \"%s\" not implemented", bearing_surf(i).options.matrix_type);
    endswitch

    switch (bearing_surf(i).options.matrix_type)
      case "nodal substruct total"
        use_total_substruct = true;
      otherwise
        use_total_substruct = false;
    endswitch

    if (~use_total_substruct)
      C_U = zeros(numel(bearing_surf(i).nodes), chunk_size);
      C_S = zeros(prod(Nxz), numel(bearing_surf(i).idx_load_case));
    endif

    dof_idx = dof_map.ndof(bearing_surf(i).nodes, 1:3);

    if (any(any(find(dof_idx <= 0))))
      error("nodes at bearing surface must not be constraint");
    endif

    joint_idx = -1;

    for j=1:numel(mesh.elements.joints)
      if (numel(mesh.elements.joints(j).nodes) == 1 && ...
          mesh.elements.joints(j).nodes == bearing_surf(i).master_node_no)

        if (joint_idx ~= -1)
          error("too many joints for bearing surface %d with master_node_no %d", i, bearing_surf(i).master_node_no);
        endif

        joint_idx = j;
      endif
    endfor

    if (joint_idx < 0)
      error("cannot find joint for bearing surface %d with master_node_no %d", i, bearing_surf(i).master_node_no);
    endif

    idx_lambda = dof_map.edof.joints(joint_idx, :);

    idx_lambda = idx_lambda(find(idx_lambda > 0));

    idx_dof_master = dof_map.ndof(bearing_surf(i).master_node_no, :);
    idx_ndof_master_act = find(idx_dof_master > 0);
    idx_dof_master_act = idx_dof_master(idx_ndof_master_act);

    if (use_total_substruct)
      idx_dof_ref_node = dof_map.ndof(cms_opt.nodes.modal.number, :);
    else
      idx_dof_ref_node = idx_dof_master;
    endif

    idx_ndof_ref_node_act = find(idx_dof_ref_node > 0);
    idx_dof_ref_node_act = idx_dof_ref_node(idx_ndof_ref_node_act);

    have_pastix = fem_sol_check_func("pastix");
    have_mumps = fem_sol_check_func("mumps");
    have_umfpack = fem_sol_check_func("umfpack");

    if (have_pastix)
      opt_pastix.matrix_type = PASTIX_API_SYM_YES;
      opt_pastix.factorization = PASTIX_API_FACT_LDLT;
      opt_pastix.refine_max_iter = cms_opt.refine_max_iter;
      opt_pastix.bind_thread_mode = PASTIX_API_BIND_NO;
      opt_pastix.number_of_threads = cms_opt.number_of_threads;
      opt_pastix.verbose = cms_opt.verbose;
    elseif (have_mumps)
      opt_mumps.matrix_type = MUMPS_MAT_SYM;
      opt_mumps.refine_max_iter = cms_opt.refine_max_iter;

      if (cms_opt.verbose)
        opt_mumps.verbose = MUMPS_VER_WARN;
      else
        opt_mumps.verbose = MUMPS_VER_ERR;
      endif

      opt_mumps.workspace_inc = int32(50);
    elseif (have_umfpack)
      opt_umfpack.refine_max_iter = cms_opt.refine_max_iter;
    endif

    if (~use_total_substruct)
      K = mat_ass.K; ## Make a copy of the stiffness matrix only but don't copy mat_ass itself
      ## Save stiffness matrix
      CC_T = K(idx_dof_master_act, idx_lambda);
      CC = K(idx_lambda, idx_dof_master_act);
      LAMBDA = K(idx_lambda, idx_lambda);
      ## Disable those joint which clamps the master node of the current bearing
      ## Rigid body motion has to be subtracted later
      ## If we don't do this, our compliance matrix will be distorted
      K(idx_dof_master_act, idx_lambda) = 0;
      K(idx_lambda, idx_dof_master_act) = 0;
      K(idx_lambda, idx_lambda) = eye(numel(idx_lambda));

      if (have_pastix)
        Kfact = fem_fact_pastix(K, opt_pastix);
      elseif (have_mumps)
        Kfact = fem_fact_mumps(K, opt_mumps);
      elseif (have_umfpack)
        Kfact = fem_fact_umfpack(K, opt_umfpack);
      else
        Kfact = fem_fact_lu(K);
      endif

      ## Restore stiffness matrix
      K(idx_dof_master_act, idx_lambda) = CC_T;
      K(idx_lambda, idx_dof_master_act) = CC;
      K(idx_lambda, idx_lambda) = LAMBDA;
      clear CC CC_T LAMBDA;
    endif
    
    A = zeros(rows(dof_idx) * 3, 6);

    if (~use_total_substruct)
      A_X0 = bearing_surf(i).X0;
    else
      A_X0 = mesh.nodes(cms_opt.nodes.modal.number, 1:3).';
    endif

    A_dX = mesh.nodes(bearing_surf(i).nodes, 1:3).' - A_X0;

    for j=1:rows(dof_idx)
      A((j - 1) * 3 + (1:3), :) = [eye(3), -skew(A_dX(:, j))];
    endfor

    clear A_dX A_X0;

    fres = 0;

    if (use_modal_contrib)
      if (~have_PHI)
        ## In order to build a proper 'a' matrix our set of mode shapes has to be extended
        ## by rigid body motions of the modal node.
        ## That's because the modal node is assumed to be constraint when 'Tred' was built.
        phi_rb = repmat(eye(6), rows(mesh.nodes), 1);
        dx_m = mesh.nodes(:, 1:3) - mesh.nodes(cms_opt.nodes.modal.number, 1:3);

        phi_rb(1:6:end, 5) =  dx_m(:, 3);
        phi_rb(1:6:end, 6) = -dx_m(:, 2);
        phi_rb(2:6:end, 4) = -dx_m(:, 3);
        phi_rb(2:6:end, 6) =  dx_m(:, 1);
        phi_rb(3:6:end, 4) =  dx_m(:, 2);
        phi_rb(3:6:end, 5) = -dx_m(:, 1);

        PHI = zeros(columns(mat_ass.K), columns(phi_rb) + columns(mat_ass.Tred));
        PHI(dof_map.idx_node, columns(phi_rb) + 1:end) = mat_ass.Tred;

        for k=1:columns(phi_rb)
          idx_ndof_act = find(dof_map.ndof(:, k) > 0);
          PHI(dof_map.ndof(idx_ndof_act, k), 1:columns(phi_rb)) = phi_rb(columns(phi_rb) * (idx_ndof_act - 1) + k, :);
        endfor

        have_PHI = true;
      endif

      Umaster = zeros(6, columns(PHI));

      Umaster(idx_ndof_ref_node_act, :) = PHI(idx_dof_ref_node_act, :);

      Unoderb = A * Umaster;

      Ustatn = zeros(rows(dof_idx), columns(PHI), 3);

      for k=1:columns(dof_idx)
        Ustatn(:, :, k) = PHI(dof_idx(:, k), :) - Unoderb(k:3:end, :);
      endfor

      D_U = zeros(numel(bearing_surf(i).nodes), columns(PHI));
      E_S = zeros(columns(PHI), numel(bearing_surf(i).idx_load_case));

      for k=1:rows(n)
        for l=1:columns(dof_idx)
          D_U += diag(n(k, :)) * Ustatn(:, :, l) * bearing_surf(i).R(l, k);
        endfor
      endfor

      a = zeros(columns(PHI), numel(bearing_surf(i).idx_load_case));
    endif

    for j=1:chunk_size:numel(bearing_surf(i).idx_load_case)
      idx_rhs = j:min(numel(bearing_surf(i).idx_load_case), (j - 1 + chunk_size));

      Rstat = mat_ass.R(:, cms_opt.num_modes_itf + bearing_surf(i).idx_load_case(idx_rhs));

      if (~use_total_substruct)
        Ustat = Kfact \ Rstat;

        ## Always consider Umaster, also for the total substruct model, because joints at the modal node could be inactive!
        Umaster = zeros(6, columns(Ustat));
        Umaster(idx_ndof_ref_node_act, :) = Ustat(idx_dof_ref_node_act, :);
        Unoderb = A * Umaster;

        Ustatn = zeros(rows(dof_idx), columns(Ustat), 3);

        for k=1:columns(dof_idx)
          Ustatn(:, :, k) = Ustat(dof_idx(:, k), :) - Unoderb(k:3:end, :);
        endfor

        C_U(:, 1:numel(idx_rhs)) = 0;
        
        for k=1:rows(n)
          for l=1:columns(dof_idx)
            C_U(:, 1:numel(idx_rhs)) += diag(n(k, :)) * Ustatn(:, :, l) * bearing_surf(i).R(l, k);
          endfor
        endfor

        clear Ustatn;
      endif

      if (use_modal_contrib)
        if (~use_total_substruct)
          persistent idx_M = [2, 1, 1;
                              3, 3, 2];

          persistent sign_M = [1, -1, 1;
                               -1, 1, -1];

          persistent idx_X = [3, 3, 2;
                              2, 1, 1];

          X_master = mesh.nodes(bearing_surf(i).master_node_no, 1:3).';
          F_master = zeros(3, numel(idx_rhs));
          M_master = zeros(3, numel(idx_rhs));
          dX_k = mesh.nodes(bearing_surf(i).nodes, 1:3) - X_master.';

          for k=1:3
            dof_idx_k = dof_map.ndof(bearing_surf(i).nodes, k);
            idx_act_k = find(dof_idx_k > 0);
            F_k = zeros(numel(dof_idx_k), numel(idx_rhs));
            F_k(idx_act_k, :) = Rstat(dof_idx_k(idx_act_k), :);
            F_master(k, :) += sum(F_k, 1);

            for l=1:2
              M_master(idx_M(l, k), :) += sign_M(l, k) * sum(F_k .* dX_k(:, idx_X(l, k)), 1);
            endfor
          endfor

          R_master = sparse([], [], [], dof_map.totdof, numel(idx_rhs));

          idx_master_dof = dof_map.ndof(bearing_surf(i).master_node_no, 1:6);
          idx_act_master_dof = find(idx_master_dof > 0);
          R_master(idx_master_dof(idx_act_master_dof), :) = [F_master;
                                                             M_master](idx_act_master_dof, :);

          U_master = Kfact \ R_master; ## Interface modes of the current master node are not constraint in Kfact
          a(:, idx_rhs) = PHI(dof_map.idx_node, :) \ U_master(dof_map.idx_node, :);
          E_S(:, idx_rhs) = PHI.' * (Rstat - R_master);
          fres = max(fres, max(norm(U_master(dof_map.idx_node, :) - PHI(dof_map.idx_node, :) * a(:, idx_rhs), "columns") ./ norm(U_master(dof_map.idx_node, :), "columns")));
        else
          E_S(:, idx_rhs) = PHI.' * Rstat;
        endif
      endif

      if (~use_total_substruct)
        for l=1:numel(idx_rhs)
          w_U = repmat(C_U(:, l), 3, 1);

          w_S = griddata_prepared(xi_U, ...
                                  zi_U, ...
                                  w_U, ...
                                  xx_S, ...
                                  zz_S, ...
                                  tri_U, ...
                                  "linear/nearest");

          if (~all(all(isfinite(w_S))))
            error("interpolation of unstructured grid failed");
          endif

          for k=1:Nxz(1)
            C_S((k - 1) * Nxz(2) + (1:Nxz(2)), idx_rhs(l)) = w_S(:, k);
          endfor
        endfor
      endif
      
      clear Unoderb Umaster Ustatn Rstat Ustat R_master U_master F_master M_master F_k idx_act_dof dof_idx_k;
    endfor

    if (use_modal_contrib && ~use_total_substruct && cms_opt.verbose)
      fprintf(stderr, "%d: modal residual error: %e\n", i, fres);
    endif

    if (fres > tol_modal_res)
      warning("modal residual error %g for bearing %d out of tolerance %g", fres, i, tol_modal_res);
    endif

    clear A Kfact fres;

    if (use_modal_contrib)
      D_S = zeros(prod(Nxz), columns(D_U));

      for j=1:columns(D_U)
        w_U = repmat(D_U(:, j), 3, 1);

        w_S = griddata_prepared(xi_U, ...
                                zi_U, ...
                                w_U, ...
                                xx_S, ...
                                zz_S, ...
                                tri_U, ...
                                "linear/nearest");

        if (~all(all(isfinite(w_S))))
          error("interpolation of unstructured grid failed");
        endif

        for k=1:Nxz(1)
          D_S((k - 1) * Nxz(2) + (1:Nxz(2)), j) = w_S(:, k);
        endfor
      endfor

      clear D_U;
    endif
    
    if (use_modal_contrib)
      D_S = D_S(:, columns(phi_rb) + 1:end); ## Remove rigid body modes not present in our .fem file
      E_S = E_S(columns(phi_rb) + 1:end, :);

      if (~use_total_substruct)
        a = a(columns(phi_rb) + 1:end, :);
        C_S -= D_S * (a + mat_ass.Kred \ E_S);
        clear a;
      endif
    endif

    if (use_modal_contrib)
      comp_mat(i).D = D_S;
      comp_mat(i).E = [E_S, E_S(:, 1:Nxz(2))];
    endif
    
    if (~use_total_substruct)
      comp_mat(i).condest_C = condest(C_S(1:end - Nxz(2), :));
      comp_mat(i).C = [C_S, C_S(:, 1:Nxz(2))];
    endif
  endfor
endfunction

%!demo
%! close all;
%! fd = -1;
%! filename = "";
%! unwind_protect
%! filename = tempname();
%! if (ispc())
%!   filename(filename == "\\") = "/";
%! endif
%! unwind_protect
%! [fd, msg] = fopen([filename, ".geo"], "wt");
%! if (fd == -1)
%!   error("failed to open file \"%s.geo\"", filename);
%! endif
%! ri = 8e-3;
%! ro = 10e-3;
%! h = 12e-3;
%! c = 2e-3;
%! b = h - 2 * c;
%! mesh_size = 1.5e-3;
%! fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%! fprintf(fd, "ri = %g;\n", ri);
%! fprintf(fd, "ro = %g;\n", ro);
%! fprintf(fd, "h = %g;\n", h);
%! fprintf(fd, "c = %g;\n", c);
%! fputs(fd, "Point(1) = {ri,0.0,0.0};\n");
%! fputs(fd, "Point(2) = {ro,0.0,0.0};\n");
%! fputs(fd, "Point(3) = {ro,0.0,c};\n");
%! fputs(fd, "Point(4) = {ro,0.0,h - c};\n");
%! fputs(fd, "Point(5) = {ro,0.0,h};\n");
%! fputs(fd, "Point(6) = {ri,0.0,h};\n");
%! fputs(fd, "Point(7) = {ri,0.0,h - c};\n");
%! fputs(fd, "Point(8) = {ri,0.0,c};\n");
%! fputs(fd, "Line(1) = {1,2};\n");
%! fputs(fd, "Line(2) = {2,3};\n");
%! fputs(fd, "Line(3) = {3,4};\n");
%! fputs(fd, "Line(4) = {4,5};\n");
%! fputs(fd, "Line(5) = {5,6};\n");
%! fputs(fd, "Line(6) = {6,7};\n");
%! fputs(fd, "Line(7) = {7,8};\n");
%! fputs(fd, "Line(8) = {8,1};\n");
%! fputs(fd, "Line Loop(5) = {1,2,3,4,5,6,7,8};\n");
%! fputs(fd, "Plane Surface(6) = {5};\n");
%! fputs(fd, "tmp[] = Extrude {{0, 0, 1}, {0, 0, 0}, 2*Pi} { Surface{6}; };\n");
%! fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%! fputs(fd, "Physical Surface(\"clamp\",1) = {tmp[2]};\n");
%! fputs(fd, "Physical Surface(\"load1\",2) = {tmp[4]};\n");
%! fputs(fd, "Physical Surface(\"load2\",3) = {tmp[8]};\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! pid = spawn("gmsh", {"-format", "msh2", "-3", "-order", "2", "-optimize_ho", "-clmin", sprintf("%g", 0.75 * mesh_size), "-clmax", sprintf("%g", 1.25 *mesh_size), [filename, ".geo"]});
%! status = spawn_wait(pid);
%! if (status ~= 0)
%!  error("gmsh failed with status %d", status);
%! endif

%! mesh = fem_pre_mesh_import([filename, ".msh"], "gmsh");
%! fprintf(stderr, "%d nodes\n", rows(mesh.nodes));
%! grp_id_clamp = find([[mesh.groups.tria6].id] == 1);
%! grp_id_p1 = find([[mesh.groups.tria6].id] == 3);
%! grp_id_p2 = find([[mesh.groups.tria6].id] == 2);
%! node_idx_modal = rows(mesh.nodes) + 1;
%! node_idx_itf1 = rows(mesh.nodes) + 2;
%! node_idx_itf2 = rows(mesh.nodes) + 3;
%! mesh.nodes(node_idx_itf2, 1:3) = [0, 0, c + 0.5 * b];
%! mesh.nodes(node_idx_itf1, 1:3) = [0, 0, c + 0.5 * b];
%! mesh.nodes(node_idx_modal, :) = 0;
%! mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [grp_id_p1, grp_id_p2], [node_idx_itf1, node_idx_itf2], "tria6");
%! bearing_surf(1).group_idx = grp_id_p1;
%! bearing_surf(1).options.reference_pressure = 1e9;
%! bearing_surf(1).options.mesh_size = 1.5 * mesh_size;
%! bearing_surf(1).options.bearing_type = "shell";
%! bearing_surf(1).options.matrix_type = "nodal";
%! bearing_surf(1).master_node_no = node_idx_itf1;
%! bearing_surf(1).nodes = mesh.groups.tria6(grp_id_p1).nodes;
%! bearing_surf(1).r = ri;
%! bearing_surf(1).w = b;
%! bearing_surf(1).X0 = [0; 0; b/2 + c];
%! bearing_surf(1).R = eye(3);
%! bearing_surf(1).relative_tolerance = 0;
%! bearing_surf(1).absolute_tolerance = sqrt(eps) * ri;
%! bearing_surf(2).group_idx = grp_id_p2;
%! bearing_surf(2).options.reference_pressure = 1e9;
%! bearing_surf(2).options.mesh_size = 1.5 * mesh_size;
%! bearing_surf(2).options.bearing_type = "journal";
%! bearing_surf(2).options.matrix_type = "nodal";
%! bearing_surf(2).master_node_no = node_idx_itf2;
%! bearing_surf(2).nodes = mesh.groups.tria6(grp_id_p2).nodes;
%! bearing_surf(2).r = ro;
%! bearing_surf(2).w = b;
%! bearing_surf(2).X0 = [0; 0; b/2 + c];
%! bearing_surf(2).R = eye(3);
%! bearing_surf(2).relative_tolerance = 0;
%! bearing_surf(2).absolute_tolerance = sqrt(eps) * ri;
%! [load_case_press, bearing_surf] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);
%! load_case.locked_dof = false(rows(mesh.nodes), 6);
%! load_case.locked_dof(mesh.groups.tria6(grp_id_clamp).nodes, :) = true;
%! load_case = fem_pre_load_case_merge(load_case, load_case_press);
%! mesh.materials.tet10 = ones(rows(mesh.elements.tet10), 1, "int32");
%! mesh.material_data.E = 210000e6;
%! mesh.material_data.nu = 0.3;
%! mesh.material_data.rho = 7850;
%! mesh.material_data.C = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%! cms_opt.nodes.modal.number = node_idx_modal;
%! cms_opt.nodes.modal.name = "node_id_modal";
%! cms_opt.nodes.interfaces(1).number = node_idx_itf1;
%! cms_opt.nodes.interfaces(1).name = "node_id_itf1";
%! cms_opt.nodes.interfaces(2).number = node_idx_itf2;
%! cms_opt.nodes.interfaces(2).name = "node_id_itf2";
%! cms_opt.number_of_threads = int32(4);
%! [mesh, mat_ass, dof_map, sol_eig, cms_opt] = fem_cms_create(mesh, load_case, cms_opt);
%! comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, mat_ass, dof_map, cms_opt, bearing_surf);
%! opt_plot.plot_nodal = true;
%! opt_plot.contour_levels = int32(15);
%! opt_plot.plot_step.x = int32(2);
%! opt_plot.plot_step.z = int32(2);
%! fem_ehd_pre_comp_mat_plot(comp_mat, opt_plot);
%! figure_list();
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
