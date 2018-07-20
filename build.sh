#!/bin/bash

# R refuses to build packages that mark themselves as Priority: Recommended
mv DESCRIPTION DESCRIPTION.old
grep -v '^Priority: ' DESCRIPTION.old > DESCRIPTION

# force anaconda not to use system libraries
echo ".libPaths(.Library)" > ".Rprofile"
$R CMD INSTALL --build .

