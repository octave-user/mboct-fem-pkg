function sol = fem_sol_static_incr(mesh, dof_map, load_case, t, opt_sol)
  R_prev = zeros(dof_map.totdof, 1);
  sol.def = zeros(rows(mesh.nodes), columns(mesh.nodes), numel(t));
  sol_prev.def = zeros(size(mesh.nodes));

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

    sol_i = fem_sol_static(mesh_i, dof_map, mat_ass, opt_sol);

    sol_prev.def += sol_i.def;
    sol.def(:, :, i) = sol_prev.def;

    fprintf(stderr, "t=%.2f U=%.3e\n", t(i), max(abs(sol_prev.def(:))));
  endfor
endfunction
