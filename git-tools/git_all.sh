#!/usr/bin/env bash

# ------------------------------------------------------------------------------
#
# Date: 170313
# Author: Gabriele Girelli
# Description:	Push/Pull all repositories in dir.
#
# ------------------------------------------------------------------------------



# PARAM ========================================================================

# Help string
helps="
 usage: ./git_all.sh [-h][-i dir]

 Description:
  Pull all repositories in dir.

 Mandatory arguments:
 -m mode	Either 'push' or 'pull'.

 Optional arguments:
  -h	Show this help page.
  -i dir	Directory containing repositories. Default: '.'
"

# Default values
dir='.'

# Parse options
while getopts hi:m: opt; do
	case $opt in
		h)
			echo -e "$helps"
			exit 0
		;;
		i)
			if [ -d $OPTARG ]; then
				dir=$OPTARG
			else
				msg="Invalid option -i, folder not found.\nFolder: $OPTARG"
				echo -e "$helps\n!!! $msg"
				exit 1
			fi
		;;
		m)
			if [ "pull" == "$OPTARG" -o "push" == "$OPTARG" ]; then
				mode=$OPTARG
			else
				msg="Invalid option -m. Possible values: pull, push."
				echo -e "$helps\n!!! $msg"
				exit 1
			fi
		;;
	esac
done

# Check mandatory arguments
if [ -z "$mode" ]; then
	msg="Missing mandatory -m option."
	echo -e "$helps\n!!! $msg"
	exit 1
fi

# RUN ==========================================================================

# Set mode
case $mode in
	'push')
		function do_action() {
			echo -e "> Pushing '$file'..."
			git push
		}
	;;
	'pull')
		function do_action() {
			echo -e "> Pulling '$file'..."
			git pull
		}
	;;
esac

# Run on repos
for	file in `ls`; do
	if [ -d "$file/.git" ]; then
		cd $file
		do_action
		cd ..
		echo ""
	fi
done

echo -e "Done."

# END ==========================================================================

################################################################################
