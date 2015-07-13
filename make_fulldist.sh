#!/bin/sh
#
# This scripts make a full .tar.gz that contains dependencies for
# release.  This is only needed to create an easy distribution if
# installing networkx and cloning pcd is too much.

VER=`git describe HEAD --always`
echo "Creating distribution version $VER"

DIST="dynbench-fulldist-$VER"
mkdir $DIST
DIST_FULLPATH=$PWD/$DIST

cp --parents `git ls-files` $DIST


# In subshell
(
    mkdir $DIST_FULLPATH/pcd/
    cd /home/darstr1/proj/pcd/ ;
    cp --parents `git --git-dir=/home/darstr1/proj/pcd/.git ls-files` $DIST_FULLPATH/pcd/
)

cp -r ~/modules/install/networkx-1.9.1/networkx/ $DIST_FULLPATH/networkx/

tar czf $DIST.tar.bz2 $DIST/
