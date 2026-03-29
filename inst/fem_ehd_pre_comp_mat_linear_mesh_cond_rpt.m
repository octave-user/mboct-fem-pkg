## Copyright (C) 2026(-2026) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(@var{cond_info})
## Print the condition information returned from fem_ehd_pre_comp_mat_linear_mesh
##
## @var{mesh} @dots{} Finite Element mesh data structure returned from fem_cms_create.
## @seealso{fem_ehd_pre_comp_mat_linear_mesh}
## @end deftypefn

function fem_ehd_pre_comp_mat_linear_mesh_cond_rpt(cond_info)
  for i=1:numel(cond_info)
    printf("size(D)=[%d,%d]\n", cond_info(i).D_size(1), cond_info(i).D_size(2));
    printf("rank(D)=%d\n", cond_info(i).D_rank);
    printf("cond(D'*diag(A)*D)=%.2e\n", cond_info(i).D_cond);

    if (cond_info(i).D_size(2) == cond_info(i).D_rank)
      printf("info: D matrix has full rank.\n");
    else
      printf("warning: D matrix does not have full rank!\n");
    endif

    if (cond_info(i).D_cond < 1e6)
      printf("info: The condition number of D'*diag(A)*D is good!\n");
    endif

    if (cond_info(i).D_cond > 1e10)
      printf("warning: The condition number of D'*diag(A)*D is bad!\n");
    endif

    printf("min(eta)=%.2e\n", min(cond_info(i).eta));
  endfor
endfunction
