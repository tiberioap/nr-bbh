#!/usr/bin/gnuplot

reset

out_scalar = 0.006

key = 1
Nf = 100
n = (Nf-1)*key

set term png

if(out_scalar != 0.0){
	system("echo '\nplotting scalar fields...\n'")
}

set size ratio 1

unset xtics
unset ytics
unset key
set colorbox
unset border
unset grid
unset xl
unset yl


# chi:_______________________________________________________________________________________________
if(out_scalar == 0.1 || out_scalar == 0.02 || out_scalar == 0.03 || out_scalar == 0.006){

	load '../palette/inferno.pal'

	do for[t=0:n]{

		Name = sprintf("../database/chi%d.dat", t)

		set output sprintf("../plots/chi%d.png", t)

		plot Name u 1:2:($3==NaN? 0: $3) matrix w image notitle
	}
}

# alpha:_____________________________________________________________________________________________
if(out_scalar == 0.2 || out_scalar == 0.02 || out_scalar == 0.06 || out_scalar == 0.006){

	load '../palette/light.pal'

	do for[t=0:n]{

		Name = sprintf("../database/alpha%d.dat", t)

		set output sprintf("../plots/alpha%d.png", t)

		plot Name u 1:2:3 matrix w image notitle
	}
}

# K:_________________________________________________________________________________________________
if(out_scalar == 0.3 || out_scalar == 0.03 || out_scalar == 0.06 || out_scalar == 0.006){

	load '../palette/matlab.pal'

	set cbr[-0.1:0.6]

	do for[t=0:n]{

		Name = sprintf("../database/K%d.dat", t)

		set output sprintf("../plots/K%d.png", t)

		plot Name u 1:2:3 matrix w image notitle
	}
}
