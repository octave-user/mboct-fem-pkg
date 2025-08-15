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
