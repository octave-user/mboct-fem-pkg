## fem_pre_load_case_merge.m:01
%!test
%! load_case1.locked_dof = false(10, 6);
%! load_case2.pressure.tria6 = zeros(5, 6, "int32");
%! load_case3(1).joints.U = zeros(3, 1);
%! load_case3(2).joints.U = ones(3, 1);
%! load_case = fem_pre_load_case_merge(load_case1, load_case2, load_case3);
%! assert_simple(numel(load_case), 4);
%! assert_simple(isfield(load_case, "locked_dof"));
%! assert_simple(isfield(load_case, "pressure"));
%! assert_simple(isfield(load_case, "joints"));
%! assert_simple(all(all(load_case(1).locked_dof == load_case1.locked_dof)));
%! assert_simple(all(all(load_case(2).pressure.tria6 == load_case2.pressure.tria6)));
%! assert_simple(all(all(load_case(3).joints.U == load_case3(1).joints.U)));
%! assert_simple(all(all(load_case(4).joints.U == load_case3(2).joints.U)));
