## Copyright (C) 2021(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {} @var{fact} = fem_fact_pastix(@var{A}, @var{opts})
## Create a factor object which uses PaStiX solve @var{A} * x = b via @var{x} = @var{Afact} \ @var{b}
## @seealso{pastix}
## @end deftypefn

function fact = fem_fact_pastix_ref(A, opts)
  if (~(nargin == 2 && ismatrix(A) && issquare(A) && isstruct(opts)))
    print_usage();
  endif

  if (~isfield(opts, "refine_max_iter"))
    opts.refine_max_iter = int32(0);
  endif

  if (~isfield(opts, "epsilon_refine"))
    opts.epsilon_refine = eps^0.8;
  endif

  if (~isfield(opts, "symmetric"))
    opts.symmetric = false;
  endif
  
  fact.opts = opts;

  opts.refine_max_iter = int32(0);

  if (opts.symmetric)
    fact.A = fem_mat_sym(A);
  else
    fact.A = A;
  endif
  
  fact.pasobj = pastix(A, opts);

  fact = class(fact, "fem_fact_pastix_ref");
endfunction
