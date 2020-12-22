set termoption dashed
t=0
do for [i=0:20] {
	file_name(n) = sprintf(".u.dat.%03d", n)
	k = i/20
	tittel = sprintf('Time t = %.3f', k);
	splot file_name(i) title tittel
	pause 0.5
}