#!/bin/bash

# Run louvain for different K's

for k in {10..50..2}
	do
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/BMMC_benchmark.csv "lou" $k --f "BMMC"
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/AML_benchmark.csv "lou" $k 
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/PANORAMA "lou" $k --f "samples"
done
# Pheno
for k in {10..50..2}
	do
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/BMMC_benchmark.csv "pheno" $k --f "BMMC"
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/AML_benchmark.csv "pheno" $k 
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/PANORAMA "pheno" $k --f "samples"
done
# MCL
for k in {10..50..2}
	do
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/BMMC_benchmark.csv "mcl" $k --f "BMMC" --i 1.1
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/AML_benchmark.csv "mcl" $k --i 1.1
	bash qsubber.sh /exports/lkeb-hpc/pderaadt/PANORAMA "mcl" $k --f "samples" --i 1.1
done
