#!/usr/bin/env bash
#
# Label the version of Cloudy to be shown on its output, using git metadata.
# The branch name will be used, along with the most recent tag, if found, or the
# SHA1 string of the most recent commit.  If the code has been modified, a
# string indicating so is also appended.  The label is dash-delimited.
#
#
# Created: Dec 05, 2020
# Author: M. Chatzikos
#
# Updated: Dec 11, 2020
# Author: M. Chatzikos
# Comment: Added support for git tags.
#
# Updated: Oct 18, 2022
# Author: M. Chatzikos
# Comment: Added support for compilation sans git (e.g., release tarballs).
#
is_repo=`git rev-parse --is-inside-work-tree 2>&1 | grep true`
if [[ $is_repo != 'true' ]]; then
	echo
	exit 1
fi
branch=`git branch | grep '^\*' | awk '{ print $2 }'`
branch=`echo $branch | sed 's/(no//'`
tag=`git describe --tags --abbrev=0 2> /dev/null`
if [ -z "$tag" ]; then
	tag=`git log --oneline | head -n 1 | awk '{print $1}'`
fi
[[ -z "`git status -s -uno`" ]] && modified="" || modified="-modified"
if [ -z "$branch" ];
then
	echo $tag$modified
else
	echo $branch-$tag$modified
fi
