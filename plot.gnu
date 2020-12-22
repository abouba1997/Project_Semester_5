set view 40, 30, 0.85, 2.
set pm3d
set isosamples 51
t=0
do for [i=0:40] {
	tittel = sprintf('{/File:Italic t} = %4.2f', i)
	set title tittel
	file_name(i) = sprintf(".u.dat.%03d", i)
	splot file_name(i) using 1:2:3 with lines
	pause 0.2
}