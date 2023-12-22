#!/bin/bash
#
# Automate the generation of a bibliography for Stout,
# and of a table showing the provenance of the species
# known to Cloudy.
#
# Chatzikos, June 14, 2018
# Chatzikos, April 14, 2022
#	Use docs/SpeciesLabels.txt instead of tsuite/auto/func_lines.spclab
#	Move PDFs to docs/ instead of copying them
# Chatzikos, Aug 17, 2023
#	Add command-line option, "all, "to delete _all_ files.  The net effect
#	is to create two PDF files in ../docs.
#

echo "Gathering references..."
./db-ref-bib2json.pl -ni
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Inspecting error files..."
[ -f empty-files.txt ] && echo "> empty-files.txt" && cat empty-files.txt && echo
[ -f broken-bibtex.txt ] && echo "> broken-bibtex.txt" && cat broken-bibtex.txt && echo
[ -f unresolved-refs.txt ] && echo "> unresolved-refs.txt" && cat unresolved-refs.txt && echo
echo
[ ! -f empty-files.txt ] && [ -f broken-bibtex.txt ] && [ -f unresolved-refs.txt ] && echo "None found!"
echo "Done."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Preparing LaTeX table for Stout..."
./db-ref-json2tex.pl -e -f=s -nr=40
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Compiling Stout table PDF..."
source mktable-stout-refs-list.sh
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

if [ ! -f ../docs/SpeciesLabels.txt ];
then
	if [ ! -x ../source/cloudy.exe ];
	then
		if [[ $OSTYPE =~ "linux" ]];
		then
			nthreads=`nproc --all`
			nht=$( expr `grep '^flags\b' /proc/cpuinfo | tail -1 | grep ht | wc -l` + 1 )
			nthreads=`echo $nthreads / $nht | bc`
		elif [[ $OSTYPE =~ "darwin" ]];
		then
			nthreads=`sysctl -n hw.physicalcpu`
		fi

		echo "Compiling Cloudy with $nthreads threads..."
		cd ../source
		make -j $nthreads > /dev/null
		[ $? == 0 ] && echo "Done."
		cd -
	fi
	
	echo "Running docs/SpeciesLabels.in..."
	cd ../docs
	../source/cloudy.exe -r LineLabels
	[ $? != 0 ] && exit "Something went wrong.  Exiting..."
	cd -
	echo; echo "Hit enter to continue..."
	read < /dev/stdin
fi

echo "Preparing table of species' provenance..."
./db-species-tex.pl -e -np=4 -nr=35 -f=s ../docs/SpeciesLabels.txt
[ $? != 0 ] && exit "Something went wrong.  Exiting..."
echo; echo "Hit enter to continue..."
read < /dev/stdin

echo "Compiling species' provenance table PDF..."
source mktable-species-db-list.sh
[ $? != 0 ] && exit "Something went wrong.  Exiting..."

echo ; echo "Stout PDF bibliography completed successfully"

echo ; echo "Cleaning up..."
if [ "$1" == "all" ]; then
	source cleanup-stout-refs-list.sh all
	source cleanup-species-db-list.sh all
	rm empty-files.txt stout-refs.tex species-db.tex

else
	source cleanup-stout-refs-list.sh
	source cleanup-species-db-list.sh
fi
echo "Done!"

echo ; echo "Moving PDFs to ../docs ..."
mv *.pdf ../docs
command="ls -trl ../docs | tail -n 2"
echo "> $command"
eval $command
echo "Done!"

echo ; echo "All done!"
