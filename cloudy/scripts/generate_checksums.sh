#!/bin/sh
echo -n "generating checksums.dat... "
rm -f checksums.dat
find . -type f -exec ../source/vh128sum.exe '{}' '+' > checksums.tmp6R4sQl
sed 's$  ./$  $' checksums.tmp6R4sQl | grep -v 'checksums.tmp6R4sQl' | LC_COLLATE=C sort -k 2 > checksums.dat
rm -f checksums.tmp6R4sQl
echo "done."
