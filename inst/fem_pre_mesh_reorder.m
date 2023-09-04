## Copyright (C) 2019(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}] = fem_pre_mesh_reorder(@var{mesh})
## @deftypefnx {} [@dots{}] = fem_pre_mesh_reorder(@dots{}, @var{options})
## Apply a fill-in reducing ordering based on METIS
##
## @var{mesh} @dots{} Finite Element Mesh
##
## @var{options}.algorithm @dots{} one of "metis", "symrcm", "csymamd" or "symamd" (default is "metis")
##
## @end deftypefn

function [mesh, PERM, IPERM] = fem_pre_mesh_reorder(mesh, options)
  if (nargin < 1 || nargin > 2 || nargout > 3)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "algorithm"))
    options.algorithm = "ndmetis";
  endif

  options.elem_type = {"tet10", "tria6", "tet10h", "tria6h", "tet20", "tria10", "tet4", "penta6", "tria3", "iso8", "iso4", "iso20", "iso27", "quad8", "iso20r", "quad8r", "penta15", "pyra5"};

  algorithm = {"ndmetis", "csymamd", "symamd", "amd", "symrcm"};

  algorithm_found = false;

  for i=1:numel(algorithm)
    if (fem_sol_check_func(options.algorithm))
      algorithm_found = true;
      break;
    endif
    options.algorithm = algorithm{i};
  endfor

  if (~algorithm_found)
    error("fem_pre_mesh_reorder failed because no algorithm is available");
  endif

  num_nz = int32(0);
  num_el = int32(0);

  for i=1:numel(options.elem_type)
    if (~isfield(mesh.elements, options.elem_type{i}))
      continue;
    endif

    elnodes = getfield(mesh.elements, options.elem_type{i});

    num_nz += numel(elnodes);
    num_el += rows(elnodes);
  endfor

  EPTR = zeros(num_el, 1, "int32");
  EIND = zeros(num_nz, 1, "int32");

  num_nz = int32(0);
  num_el = int32(0);

  for i=1:numel(options.elem_type)
    if (~isfield(mesh.elements, options.elem_type{i}))
      continue;
    endif

    elnodes = getfield(mesh.elements, options.elem_type{i});

    EPTR(num_el + (1:rows(elnodes))) = num_nz + (1:columns(elnodes):numel(elnodes));
    EIND(num_nz + (1:numel(elnodes))) = (elnodes.')(:);

    num_nz += numel(elnodes);
    num_el += rows(elnodes);
  endfor

  switch (options.algorithm)
    case "ndmetis"
      [PERM, IPERM] = feval(options.algorithm, rows(mesh.nodes), EPTR, EIND);
    case {"symrcm", "csymamd", "symamd", "amd"}
      num_nz = sum((EPTR(2:end) - EPTR(1:end-1)).^2);

      RIDX = CIDX = zeros(num_nz, 1, "int32");

      num_nz = int32(0);

      for i=2:rows(EPTR)
        num_nodes = EPTR(i) - EPTR(i - 1);
        for j=1:num_nodes
          for k=1:num_nodes
            ++num_nz;
            RIDX(num_nz) = EIND(EPTR(i - 1) + j);
            CIDX(num_nz) = EIND(EPTR(i - 1) + k);
          endfor
        endfor
      endfor

      A = sparse(RIDX, CIDX, ones(size(RIDX)), rows(mesh.nodes), rows(mesh.nodes));
      PERM = int32(feval(options.algorithm, A)(:));
      IPERM(PERM) = int32(1:rows(PERM));
      IPERM = IPERM(:);
      clear A RIDX CIDX;
    otherwise
      error("unknown value for options.algorithm = \"%s\"", options.algorithm);
  endswitch

  mesh.nodes = mesh.nodes(PERM, :);

  for i=1:numel(options.elem_type)
    if (~isfield(mesh.elements, options.elem_type{i}))
      continue;
    endif

    elnodes = getfield(mesh.elements, options.elem_type{i});

    mesh.elements = setfield(mesh.elements, options.elem_type{i}, reshape(IPERM(elnodes), size(elnodes)));
  endfor

  if (isfield(mesh, "groups"))
    for i=1:numel(options.elem_type)
      if (~isfield(mesh.groups, options.elem_type{i}))
        continue;
      endif

      grp = getfield(mesh.groups, options.elem_type{i});

      for j=1:numel(grp)
        grp(j).nodes = reshape(IPERM(grp(j).nodes), size(grp(j).nodes));
      endfor

      mesh.groups = setfield(mesh.groups, options.elem_type{i}, grp);
    endfor
  endif
endfunction

%!test
%! mesh.nodes = [0, 0, 0;
%!               0, 1, 0;
%!               0, 2, 0;
%!               1, 0, 0;
%!               1, 1, 0;
%!               1, 2, 0;
%!               2, 0, 0;
%!               2, 1, 0;
%!               2, 2, 0];
%!
%! mesh.elements.iso4 = int32([1,2,5,4;
%!                             2,3,6,5;
%!                             4,5,8,7;
%!                             5,6,9,8]);
%! mesh.groups.iso4(1).elements = int32([1, 2]);
%! mesh.groups.iso4(2).elements = int32([3, 4]);
%! mesh.groups.iso4(1).id = 1;
%! mesh.groups.iso4(2).id = 2;
%! mesh.groups.iso4(1).name = "group1";
%! mesh.groups.iso4(2).name = "group2";
%!
%! for i=1:numel(mesh.groups.iso4)
%!   mesh.groups.iso4(i).nodes = sort(unique(mesh.elements.iso4(mesh.groups.iso4(i).elements, :)(:)));
%! endfor
%!
%! x = rand(rows(mesh.nodes), 10);
%! s = zeros(rows(mesh.elements.iso4), 1);
%!
%! for i=1:rows(mesh.elements.iso4)
%!   s(i) = sum(x(mesh.elements.iso4(i, :)));
%! endfor
%!
%! g = zeros(numel(mesh.groups.iso4), 1);
%!
%! for i=1:numel(mesh.groups.iso4)
%!   g(i) += sum(x(mesh.groups.iso4(i).nodes));
%! endfor
%!
%! algorithm = {"metis", "amd", "symrcm", "csymamd", "symamd"};
%!
%! for j=1:numel(algorithm)
%!   [mesh2, perm, iperm] = fem_pre_mesh_reorder(mesh, struct("algorithm", algorithm{j}));
%!
%!   x2 = x(perm, :);
%!   s2 = zeros(rows(mesh2.elements.iso4), 1);
%!
%!   for i=1:rows(mesh2.elements.iso4)
%!     s2(i) = sum(x2(mesh2.elements.iso4(i, :)));
%!   endfor
%!
%!   A = zeros(rows(mesh.nodes), rows(mesh.nodes));
%!
%!   for i=1:rows(mesh.elements.iso4)
%!     A(mesh.elements.iso4(i, :), mesh.elements.iso4(i, :)) += x(mesh.elements.iso4(i, :));
%!   endfor
%!
%!   A2 = zeros(rows(mesh2.nodes), rows(mesh2.nodes));
%!
%!   for i=1:rows(mesh2.elements.iso4)
%!     A2(mesh2.elements.iso4(i, :), mesh2.elements.iso4(i, :)) += x2(mesh2.elements.iso4(i, :));
%!   endfor
%!
%!   assert(s2, s);
%!   assert(mesh2.nodes(iperm,:), mesh.nodes)
%!   assert(A2(iperm, iperm), A);
%!   assert(A2, A(perm, perm));
%!   for i=1:numel(mesh.groups.iso4)
%!     assert(mesh2.groups.iso4(i).elements, mesh.groups.iso4(i).elements);
%!     assert(mesh2.groups.iso4(i).name, mesh.groups.iso4(i).name);
%!     assert(perm(mesh2.groups.iso4(i).nodes), mesh.groups.iso4(i).nodes);
%!   endfor
%!   g2 = zeros(numel(mesh2.groups.iso4), 1);
%!   for i=1:numel(mesh2.groups.iso4)
%!    g2(i) += sum(x2(mesh2.groups.iso4(i).nodes));
%!   endfor
%! endfor
