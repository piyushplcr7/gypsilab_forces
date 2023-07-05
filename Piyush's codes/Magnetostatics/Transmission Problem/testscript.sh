 
for k in {1..1}
do
	export KK=$k
	for a in {1..1}
	do
		export AA=$a
		for b in {1..1}
		do
			export BB=$b
			for c in {1..1}
			do
				export CC=$c
				echo "Shape derivative computations for a = $a, b = $b, c = $c, k = $k"
				outname="gen_shape_derivative_${a}_${b}_${c}_${k}"
				echo $outname
				#sbatch --ntasks=5 --time=3-1 --output=$outname --open-mode=append --mem-per-cpu=10g --wrap="matlab -nodisplay -singleCompThread -r TP_general_shape_derivative_comparison"
				echo "checking KK= ${KK}"
				matlab -nodisplay -singleCompThread -r TP_general_shape_derivative_comparison
			done
		
		done
	
	done

done
