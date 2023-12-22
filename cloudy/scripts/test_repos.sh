#!/bin/sh
\rm -f fix_properties
scriptdir=$(dirname $0)
find . -exec $scriptdir/test_properties.pl '{}' ';'
if [ -s fix_properties ]; then
    echo "---------------------------"
    echo ""
    echo "Please review the proposed changes above, keeping in mind"
    echo "that file types may not have been detected correctly..."
    echo "Pay particular attention to (perl) scripts as they will not"
    echo "be detected as such if the correct header line is missing!"
    echo "(for perl this is '#\!/usr/bin/perl' as the first line)."
    echo ""
    echo "If you agree with the proposed changes, execute this command:"
    echo ""
    echo "  source fix_properties"
    echo ""
    echo "and rerun this script to double-check the result."
    echo ""
    echo "If files have been misidentified, edit them to enable correct"
    echo "identification. This can be checked by typing:"
    echo "  file <filename>"
    echo "Then rerun this script."
fi
if [ -f fix_properties -a ! -s fix_properties ]; then
    \rm -f fix_properties
fi
