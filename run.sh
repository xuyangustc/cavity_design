#!/bin/bash

tar -xvzf energy_data.tar.gz

#first put high quality loop lib in folder named hqll as PDBformat
ls hqll > list_hqll

mkdir temp
Initial_Chain_Gen pars/par_initial_chain_gen

mkdir data/

# generate 10 structure for example, user can generate to modify ${num_generate}, parameter named {nout} in pars/par_initial_chain_gen should also be changed
num_generate=10
for (( i=0; i<$num_generate; i++ ))
do
	mkdir data/$i
	mv temp/$i data/$i/pdb_0
	cp pars/par_far data/$i/
	nres=$(tail -n 1 data/$i/pdb_0 | awk '{print $6}')
	sed -i "8c 0 0-${nres} CA" data/$i/par_far
	sed -i "53c 0 0-${nres} CA" data/$i/par_far
	sed -i "62c 0 0-${nres} CA" data/$i/par_far
	cd data/$i
	SDrun_far ../../pars/par_relax par_far
	SDrun_far ../../pars/par_sampling par_far
	SDrun_far ../../pars/par_relax_2 par_far
	SDrun_far ../../pars/par_sampling_2 par_far
	SDrun_far ../../pars/par_optimization par_far
	if [ ! -f "pdb_5" ]
	then
		cd ../../
		continue
	fi
	bash ../../SSinfo.sh pdb_5 ss_5
	Loop_replace pdb_5 ss_5 ../../list_hqll ../../hqll/ pdb_6 resfile_6 list_6
	if [ ! -f "list_6" ]
	then
		continue
	fi
	while read p6
	do
		DesignSeq -in pdb_6_$p6 -out seq6 -log log6 -n 5 -resFile resfile_6_$p6
		seq_min=$(paste log6 lb | sort -n -k 2 | head -n 1 | awk '{print $3}')
		mv seq_min pdb_7_$p6
		cp pdb_7_$p6 pdb_7
		SDrun ../../pars/par_relax_3
		SDrun ../../pars/par_sampling_3
		if [ ! -f "pdb_9" ]
		then
			continue
		fi
		mv pdb_9 pdb_9_$p6
		bash ../../SSinfo.sh pdb_9_$p6 ss_9_$p6
		Loop_replace_2 pdb_9_$p6 ss_9_$p6 ../../list_hqll ../../hqll/ out
		if [ ! -f "out" ]
		then
			rm out
			continue
		fi
		DesignSeq -in pdb_9_$p6 -out seq9 -log log9 -n 5
		seq_min=$(paste log9 lb | sort -n -k 2 | head -n 1 | awk '{print $3}')
		# final_* is the final results
		mv seq_min final_$p6
	done < list_6
	cd ../../
done

