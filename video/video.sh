#!/usr/bin/env bash
#
# Jim Teresco, CSIS 335, Siena College, Fall 2021
#
# Converts each of the output files from an output level 3 run of
# the Jacobi quadtree solver to a PNG file using gnuplot, then
# pastes them together into am MP4 movie with ffmpeg
#
# First command-line parameter is the number of refinement iterations,
# second command-line parameter is the maximum number of Jacobi
# iterations for any refinement iteration.
#
if [ $# -ne 2 ]; then
    echo "Usage: video.sh num_ref max_jac"
    exit 1;
fi
cp /dev/null images.txt
for r in `seq 1 $1`; do
    for i in `seq 1 $2`; do
	if [ -f "solution.$r.$i.dat" ]; then
	    echo "processing solution.$r.$i.dat"
	    cp solution.$r.$i.dat solution.dat
	    gnuplot solutionplot.gp
	    mv solution.png solution.$r.$i.png
	    echo "file 'solution.$r.$i.png'" >> images.txt
	fi
    done
done
ffmpeg -y -r 30 -f concat -safe 0 -i "images.txt" -c:v libx264 -vf "fps=25,format=yuv420p" "solution.mp4"
