#!/usr/bin/gnuplot

reset

out_tensor = 0.0

key = 1
Nf = 100
n = (Nf-1)*key

set term png

if(out_tensor != 0.0){
	system("echo '\nplotting tensor fields...\n'")
}

unset xtics
unset ytics
unset ztics
unset key
set colorbox
unset border
unset xl
unset yl
unset zl

unset surface 
set view 130,10
set size ratio 1
set pm3d implicit at s
set pm3d scansbackward

# gammaTil_ij:_____________________________________________________
if(out_tensor == 2.1 || out_tensor == 4.62){

	load '../palette/inferno.pal'

	set view 115,130

	set xr[20:180] 	
	set yr[20:180] 	
	set zr[-20:120]
	set cbr[0:20]

	do for[t=0:n]{

		Name = sprintf("../database/gammaTil_ij%d.dat", t)

		set output sprintf("../plots/gammaTil_ij%d.png", t)

		splot Name u 1:2:3 notitle
	}
}

# ATil_ij:_________________________________________________________
if(out_tensor == 2.2 || out_tensor == 4.62){

	load '../palette/inferno.pal'

	set xr[20:180] 	
	set yr[20:180] 	
	set zr[-25:2]	
	set cbr[-2:1]

	set view 115,40

	do for[t=0:n]{

		Name = sprintf("../database/ATil_ij%d.dat", t)

		set output sprintf("../plots/ATil_ij%d.png", t)

		unset cbr
		plot Name u 1:2:3 ls 7 lc palette notitle
	}
}

