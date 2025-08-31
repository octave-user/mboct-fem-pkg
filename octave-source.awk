## Copyright (C) 2018(-2025) Reinhard <octave-user@a1.net>

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

BEGIN {
    RS = "[<>\"]";
    FPAT = "([^-.]+)";
    majorver = -1;
    minorver = -1;
    subver = -1;
}

/^octave-[0-9]+\.[0-9]+\.[0-9]+\.tar\.gz$/ {
    if ($2 > majorver) {
        majorver = $2;
        minorver = $3;
        subver = $4;
    } else if ($2 == majorver && $3 > minorver) {
        minorver = $3;
        subver = $4;
    } else if ($2 == majorver && $3 == minorver && $4 > subver) {
        subver = $4;
    }
}

END {
    if (majorver >= 0 && minorver >= 0 && subver >= 0) {
        printf("octave-%d.%d.%d.tar.gz\n", majorver, minorver, subver);
    }
}
