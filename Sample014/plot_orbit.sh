#!/bin/sh

GNUPLOT_CMD="
	set term png size 10000,10000;
	set output 'orbit.png';
	set size ratio -1;
	plot[-2:2][-2:2] 'data.dat' u 1:2 w l lw 2
"

for i in `seq 2 10`; do
	xi=`echo "${i}*2 - 1" | bc`
	yi=`echo "${i}*2" | bc`
	GNUPLOT_CMD="
		${GNUPLOT_CMD}
		, '' u ${xi}:${yi} w l
	"
done

echo ${GNUPLOT_CMD}

gnuplot -e "${GNUPLOT_CMD}"
