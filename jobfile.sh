
NPRODS=150
NCPUS=6

prod=1
for ((i=1;i<=10;i++))
do
	j=1
	while ((j <= $NCPUS)) && ((prod <= $NPRODS))
	do
		echo $j $i $prod
		((b=$prod+4))
		printf -v x "%03d" $prod
		printf -v y "%03d" $b
		mkdir $x.$y.Covar
		cd $x.$y.Covar
		sed -e s/AAA/$prod/g -e s/BBB/$b/g -e s/XXX/$x/g -e s/YYY/$y/g < ../sample.config > $x.$y.Covar.config
		time ../Covar_calc.py $x.$y.Covar.config > Covar.output &
		cd ../
		((j=$j+1))
		((prod=$prod+5))
	done
	wait
done

