function sol = fem_sol_static_incr(mesh, dof_map, load_case, t, opt_sol)
  R_prev = zeros(dof_map.totdof, 1);
  sol.def = zeros(rows(mesh.nodes), columns(mesh.nodes), numel(t));
  sol_prev.def = zeros(size(mesh.nodes));

  ## mat_ass.KTAU0 = sparse([],[],[], dof_map.totdof, dof_map.totdof);

  for i=1:numel(t)
    mesh_i = setfield(mesh, "nodes", mesh.nodes + sol_prev.def);
    [mat_ass.K, ...
     mat_ass.R] = fem_ass_matrix(mesh_i, ...
                                 dof_map, ...
                                 [FEM_MAT_STIFFNESS, ...
                                  FEM_VEC_LOAD_CONSISTENT], ...
                                 load_case);


    R_i = mat_ass.R * t(i);
    mat_ass.R = R_i - R_prev;
    R_prev = R_i;

    #mat_ass.K += mat_ass.KTAU0;

    sol_i = fem_sol_static(mesh_i, dof_map, mat_ass, opt_sol);

    sol_prev.def += sol_i.def;
    sol.def(:, :, i) = sol_prev.def;

    fprintf(stderr, "t=%.2f U=%.3e\n", t(i), max(abs(sol_prev.def(:))));
    
    ## sol_i.stress = fem_ass_matrix(mesh_i, ...
    ##                               dof_map, ...
    ##                               [FEM_VEC_STRESS_CAUCH], ...
    ##                               load_case, ...
    ##                               sol_prev);

    ## load_case.tau0 = sol_i.stress.tau;

    ## mat_ass.KTAU0 = fem_ass_matrix(mesh_i, ...
    ##                                dof_map, ...
    ##                                [FEM_MAT_STIFFNESS_TAU0], ...
    ##                                load_case);
  endfor
endfunction
