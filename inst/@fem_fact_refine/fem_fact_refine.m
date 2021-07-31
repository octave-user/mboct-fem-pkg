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
## Create a factor object with iterative refinement
## @seealso{pastix}
## @end deftypefn

function fact = fem_fact_refine(A, opts)
  if (~(nargin == 2 && ismatrix(A) && issquare(A) && isstruct(opts)))
    print_usage();
  endif

  if (~isfield(opts, "refine_max_iter"))
    opts.refine_max_iter = int32(3);
  endif

  if (~isfield(opts, "epsilon_refinement"))
    opts.epsilon_refinement = -1;
  endif

  if (opts.epsilon_refinement <= 0)
    opts.epsilon_refinement = eps^0.8;
  endif

  if (~isfield(opts, "symmetric"))
    opts.symmetric = true;
  endif
  
  if (~isfield(opts, "verbose"))
    opts.verbose = false;
  endif

  fact.opts = opts; ## Do not overwrite fact.opts

  opts.refine_max_iter = int32(0);

  if (opts.symmetric && ~issymmetric(A))
    fact.A = fem_mat_sym(A);
  else
    fact.A = A;
  endif
  
  fact.Afact = fem_sol_factor(A, opts);

  fact = class(fact, "fem_fact_refine");
endfunction
