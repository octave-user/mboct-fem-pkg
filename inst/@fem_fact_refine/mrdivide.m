## Copyright (C) 2021(-2024) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{X} = mldivide(@var{Afact}, @var{B})
## Solve @var{Afact} * @var{X} = @var{B} by using the factor object @var{Afact}.
## @end deftypefn

function X = mrdivide(B, fact)
  narginchk(2, 2);

  i = int32(0);
  X = zeros(rows(B), rows(fact.A));
  R = B;

  fconverged = false;

  do
    X += R / fact.Afact;
    AX = X * fact.A;
    R = B - AX;
    f = max(norm(R, "rows") ./ max(1, norm(B + AX, "rows"))); ## Avoid division by zero

    if (fact.opts.verbose)
      fprintf(stderr, "iteration %d: %e\n", i, f);
    endif

    if (f <= fact.opts.epsilon_refinement)
      fconverged = true;
      break;
    endif
  until (++i > fact.opts.refine_max_iter)

  if (~fconverged)
    warning("iterative refinement did not converge after %d iterations f=%.3e tol=%.3e", i, f, fact.opts.epsilon_refinement);
  endif
endfunction
