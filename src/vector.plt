#!/usr/bin/gnuplot

reset

out_vector = 1.2

key = 1
Nf = 100
n = (Nf-1)*key

set term png

if(out_vector != 0.0){
	system("echo '\nplotting vector fields...\n'")
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

# beta_I:________________________________________________________________________________
if(out_vector == 1.1 || out_vector == 1.32 || out_vector == 1.43 || out_vector == 1.716){

	load '../palette/inferno.pal'

	do for[t=0:n]{

		Name = sprintf("../database/beta_I%d.dat", t)

		set output sprintf("../plots/beta_I%d.png", t)

		plot Name u 1:2:(50*$3):(50*$4) w vec notitle
	}
}

# GammaTil_I:____________________________________________________________________________
if(out_vector == 1.2 || out_vector == 1.32 || out_vector == 1.56 || out_vector == 1.716){

	load '../palette/inferno.pal'

	do for[t=0:n]{

		Name = sprintf("../database/GammaTil_I%d.dat", t)

		set output sprintf("../plots/GammaTil_I%d.png", t)

		plot Name u 1:2:(10*$3):(10*$4) w vec notitle
	}
}

# x_I:___________________________________________________________________________________
if(out_vector == 1.3 || out_vector == 1.43 || out_vector == 1.56 || out_vector == 1.716){

	set xr[0:200]
	set yr[0:200]

	Name = sprintf("../database/x_I.dat")

	set output sprintf("../plots/x_I.png")

	plot Name u 2:1 ls 7 notitle, Name u 4:3 ls 7 notitle
}

