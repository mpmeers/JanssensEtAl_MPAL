#!/usr/bin/env bash

set -ue


if [ $# -lt 2 ]
then
	echo "
		Need at least two inputs: info file and output path!
	"
	exit 1
fi

names=0
while IFS=$'\t' read -r a b c d e f
do
	name="${d}_${e}_${f}"	
	name2="${d}_${e}_${f}_IgG"	
	output="${2}/${d}_${e}_${f}"
	output2="${2}/${d}_${e}_${f}_IgG"
	sbatch -n 1 -J $name  -Q --wrap="bash SEACR_1.3.sh $a $b norm relaxed $output"
	sbatch -n 1 -J $name2  -Q --wrap="bash SEACR_1.3.sh $b 0.01 norm stringent $output2"
	names="${names},${name},${name2}"
	echo "$name" >> $2/DAlChrE_temp.txt
	echo -e "${c}\t${e}_${f}\t${d}\t${name}" >> $2/DAlChrE_beds.txt
done < $1

finished="FALSE"
while [[ $finished == "FALSE" ]]
do
	running=`squeue -n $names | wc -l`
	if [[ $running -gt 1 ]]
	then
		let count=$running-1
		echo "Waiting for $count jobs to finish at $(date +%T)..."
		sleep 10
	else
		finished="TRUE"
		fileout=`ls $2 | grep -f $2/DAlChrE_temp.txt | grep -v "IgG" | wc -l`
		fileout2=`wc -l $2/DAlChrE_temp.txt | awk '{print $1}'`
		if [[ $fileout -lt $fileout2 ]]
		then
			let failed=$fileout2-$fileout
			echo "$failed jobs failed--would you still like to continue? (y/n)"
#			echo "All jobs finished, would you like to continue? (y/n)"
			read response
			while [[ $response != "y" ]] && [[ $response != "n" ]]
			do
				echo "Please provide a (y/n) response"
				read response
			done
			if [[ $response == "y" ]]
			then
				echo "Proceeding!"
			else [[ $response == "n" ]] 
				echo "Aborting!"
				exit 1
			fi
		else
			echo "All jobs succeeded--proceeding!"
		fi
	fi
done

#module load bedops
echo "Generating control mask file at $(date +%T)"
bedops -m $2/*_IgG.stringent.bed > $2/IgG_comb.bed
for i in `cat $2/DAlChrE_temp.txt`
do
	target="${2}/${i}.relaxed.bed"
	IgG="${2}/${i}_IgG.stringent.bed"
	output3="${2}/${i}_IgGfilter.relaxed.bed"
	thislen=`wc -l $target | awk '{print $1}'`
	if [[ $thislen -gt 0 ]]
	then
#		bedtools intersect -v -a $target -b $2/IgG_comb.bed > $output3
		cat $target > $output3
		thatlen=`wc -l $output3 | awk '{print $1}'`
		if [[ $thatlen -lt 1000 ]]
		then
			rm $output3
#			awk -v var=$i '$4 != var {print $0}' $2/DAlChrE_beds.txt > $2/DAlChrE_beds.tmp ## These lines can be omitted to preserve mapping from files for which there are few peaks
#			mv $2/DAlChrE_beds.tmp $2/DAlChrE_beds.txt                                     ## These lines can be omitted to preserve mapping from files for which there are few peaks
		fi
	fi
	rm $target
	rm $IgG
done

echo "Combining replicates at $(date +%T)"
#thesefiles=`ls $2/*_IgGfilter.relaxed.bed`
#for o in $thesefiles
#do	
#	thatlen=`wc -l $o | awk '{print $1}'`
#	if [[ $thatlen -lt 1000 ]]
#	then
#		rm $o
#	fi
#done
targetcond=`basename -a $2* | grep "_IgGfilter.relaxed.bed" | awk 'BEGIN{FS="_"}; {print $1"_"$2}' | sort -u`
for n in $targetcond
do	
	output8="${2}/${n}_IgGfilter_comb.relaxed.bed"
	echo "Writing $output8 at $(date +%T)"
	num=`ls $2/$n*_IgGfilter.relaxed.bed | wc -l`
	if [[ num -gt 1 ]]
	then
#		bedops -i $2/$n*_IgGfilter.relaxed.bed > $output8
#		bedops -m $2/$n*_IgGfilter.relaxed.bed > $output8
		bedtools multiinter -i $2/$n*_IgGfilter.relaxed.bed | awk 'BEGIN{s=1}; {if(s==1){num=$4; chr=$1; start=$2; stop=$3; s++}else{if($2==stop){stop=$3; if($4>num){num=$4}}else{print chr"\t"start"\t"stop"\t"num; num=$4; chr=$1; start=$2; stop=$3}}}' | awk -v var=$num '$4==var {print $0}' > $output8
		rm $2/$n*_IgGfilter.relaxed.bed
	else
	thatfile=`ls $2/$n*_IgGfilter.relaxed.bed`
		mv $thatfile $output8
	fi
done

factors=`awk 'BEGIN{FS="_"}; {print $1}' $2/DAlChrE_temp.txt | sort -u`

echo "Generating factor-specific count files at $(date +%T)"
for j in $factors
do
#	list2=`ls $2/* | grep -f $2/DAlChrE_temp.txt | grep "IgGfilter" | grep $j`
	list2=`ls $2/* | grep "_IgGfilter_comb.relaxed.bed" | grep $j`
#	avg=`bedops -m $list2 | awk '{s+=$3-$2; t++} END {print s/(10*t)}'`
	output4="${2}/${j}.regions.bed"
#	bedops -m $list2 > $output4
#	bedops -m $list2 | awk -v var=$avg 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; end=$3; s++}else{if(chr==$1 && $2 < end+var){end=$3}else{print chr"\t"start"\t"end; chr=$1; start=$2; end=$3}}}' > $output4
#	bedops -m $list2 | awk '{add=($3-$2)/10; print $1"\t"$2-add"\t"$3+add"\t"$2"\t"$3}' | awk 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; end=$3; realstart=$4; realend=$5; s++}else{if(chr==$1 && $2 < end){end=$3; realend=$5}else{print chr"\t"realstart"\t"realend; chr=$1; start=$2; end=$3; realstart=$4; realend=$5}}}' > $output4
	awk -v var=$j '$3==var {print $0}' $2/DAlChrE_beds.txt > $2/DAlChrE_beds.txt.$j


	count=1
	first=1
	
	while [[ $count -gt 0 ]]
	do
		out2="${2}${j}.remove.txt"
		if [[ $first -eq 1 ]]
		then
			bedtools multiinter -i $2$j*relaxed.bed | cut -f 1,2,3,4 | awk 'BEGIN{s=0}; {if(s==0){chr=$1; start=$2; end=$3; regions=$4; sum=$4*($3-$2); s++}else{if(chr==$1 && end==$2){end=$3; sum=sum+($4*($3-$2)); if($4>regions){regions=$4}}else{print chr"\t"start"\t"end"\t"regions"\t"sum/(regions*(end-start))"\t"s; chr=$1; start=$2; end=$3; regions=$4; sum=$4*($3-$2); s++}}} END {print chr"\t"start"\t"end"\t"regions"\t"sum/(regions*(end-start))"\t"s}' | awk '$5 < 0.5 {print $6}' > $out2 
		else	
			bedtools multiinter -i $2$j*temp.bed | cut -f 1,2,3,4 | awk 'BEGIN{s=0}; {if(s==0){chr=$1; start=$2; end=$3; regions=$4; sum=$4*($3-$2); s++}else{if(chr==$1 && end==$2){end=$3; sum=sum+($4*($3-$2)); if($4>regions){regions=$4}}else{print chr"\t"start"\t"end"\t"regions"\t"sum/(regions*(end-start))"\t"s; chr=$1; start=$2; end=$3; regions=$4; sum=$4*($3-$2); s++}}} END {print chr"\t"start"\t"end"\t"regions"\t"sum/(regions*(end-start))"\t"s}' | awk '$5 < 0.5 {print $6}' > $out2 
		fi
		count=`wc -l $out2 | awk '{print $1}'`
		echo "Resolving $count remaining $j regions at $(date +%T)..."
		if [[ $count -gt 0 ]]
		then
			out3="${2}${j}.keep1.bed"
			out4="${2}${j}.keep2.bed"
			out5="${2}${j}.remove.bed"
			if [[ $first -eq 1 ]] 
			then
				cat $2$j*relaxed.bed | awk 'BEGIN{s=0}; {if(s==0){s++; chr=$1; print $1"\t"$2"\t"$3"\t"$3-$2"\t"s}else{if(chr!=$1 && $1=="chr1"){s++}; chr=$1; print $1"\t"$2"\t"$3"\t"$3-$2"\t"s}}' | sort -k1,1 -k2,2n -k3,3n | bedtools cluster -i - | sort -k6,6n -k4,4nr | awk 'BEGIN{s=0}; {if(s==0){s++; print $0"\t"s; cluster=$6}else{if($6==cluster){s++; print $0"\t"s; cluster=$6}else{s=1; print $0"\t"s; cluster=$6}}}' | awk -v outkeep1=$out3 'NR==FNR {pats[$1]=0; next} {s=0; for(p in pats) if($6==p){s=1; break}; if(s==0){print > outkeep1}else{print $0}}' $out2 - | awk -v outkeep2=$out4 -v outremove=$out5 'BEGIN{s=1}; {if(s==1){print >> outremove; cluster=$6; s++}else{if(cluster!=$6){print >> outremove}else{print > outkeep2}; cluster=$6}}'
			first=2
			else	
				cat $2$j*temp.bed | awk 'BEGIN{s=0}; {if(s==0){s++; chr=$1; print $1"\t"$2"\t"$3"\t"$3-$2"\t"s}else{if(chr!=$1 && $1=="chr1"){s++}; chr=$1; print $1"\t"$2"\t"$3"\t"$3-$2"\t"s}}' | sort -k1,1 -k2,2n -k3,3n | bedtools cluster -i - | sort -k6,6n -k4,4nr | awk 'BEGIN{s=0}; {if(s==0){s++; print $0"\t"s; cluster=$6}else{if($6==cluster){s++; print $0"\t"s; cluster=$6}else{s=1; print $0"\t"s; cluster=$6}}}' | awk -v outkeep1=$out3 'NR==FNR {pats[$1]=0; next} {s=0; for(p in pats) if($6==p){s=1; break}; if(s==0){print > outkeep1}else{print $0}}' $out2 - | awk -v outkeep2=$out4 -v outremove=$out5 'BEGIN{s=1}; {if(s==1){print >> outremove; cluster=$6; s++}else{if(cluster!=$6){print >> outremove}else{print > outkeep2}; cluster=$6}}'
			fi
	#		wc -l $out3 | awk '{print $1}'
	#		wc -l $out4 | awk '{print $1}'
	#		wc -l $out5 | awk '{print $1}'
	#		cat $out3 $out4 $out5 | wc -l
			cat $out3 $out4 | sort -k1,1 -k2,2n -k3,3n | awk -v path=$2 -v target=$j '{print > (path target "." $5 ".temp.bed")}'
		
		fi
	done
	
	out6="${2}${j}_merge.bed"
	bedops -m $2$j*.temp.bed | awk '{add=($3-$2)/10; print $1"\t"$2-add"\t"$3+add"\t"$2"\t"$3}' | awk 'BEGIN{s=1}; {if(s==1){chr=$1; start=$2; end=$3; realstart=$4; realend=$5; s++}else{if(chr==$1 && $2 < end){end=$3; realend=$5}else{print chr"\t"realstart"\t"realend; chr=$1; start=$2; end=$3; realstart=$4; realend=$5}}}' > $out6
	
	out7="${2}${j}.remove.sorted.bed"
	sort -k1,1 -k2,2n -k3,3n $out5 > $out7
	
	bedtools multiinter -i $out6 $out7 | cut -f 1,2,3 > $output4
	rm $2$j.remove.txt
	rm $2$j.remove.bed
	rm $2$j.remove.sorted.bed
	rm $2$j.keep1.bed
	rm $2$j.keep2.bed
	rm $2$j*temp.bed
	rm $out6
	echo "Done generating $j regions at $(date +%T)"

	names=0
#	while IFS=$'\t' read -r bedfile condition factor
	while IFS=$'\t' read -r bedfile condition factor fullname
	do
		output5="${output4}.${condition}.bed"
		echo "Writing $output5 at $(date +%T)"
#		bedtools intersect -c -a $output4 -b $bedfile | cut -f 4 | echo -e "$condition\n`cat -`" > $output5
		sbatch -n 1 -J $output5 -Q --wrap="bedtools intersect -c -a $output4 -b $bedfile | cut -f 4 > $output5"
		names="${names},${output5}"
	done < $2/DAlChrE_beds.txt.$j
	rm $2/DAlChrE_beds.txt.$j ## Add this back

	finished="FALSE"
	while [[ $finished == "FALSE" ]]
	do
		running=`squeue -n $names | wc -l`
		if [[ $running -gt 1 ]]
		then
			let count=$running-1
			echo "Waiting for $count jobs for $j to finish at $(date +%T)..."
			sleep 10
		else
			finished="TRUE"
			echo "All jobs for $j finished at $(date +%T)!"
		fi
	done
	echo "Merging count files for $j at $(date +%T)"
	header=`echo -e "chr\tstart\tend"`
#	ls $output4.*
	list=`ls $output4.* | awk 'BEGIN{FS="."}; {print $(NF-1)}'`
	for m in $list
	do
#		top=`awk 'BEGIN{FS="."}; {print $(NF-1)}' $m`
		thisfile=`ls $output4.* | grep "$m\."`
		echo -e "$m\n`cat $thisfile`" > $thisfile.header
		rm $thisfile
	done
	output6="${output4}_counts.bed"
	echo -e "$header\n`cat $output4`" | paste - $output4.* > $output6
	rm $output4.*  ## Add this back
done
rm $2/DAlChrE_*
echo "Performing PCA and t-SNE cluster analysis at $(date +%T)"
module load R
names=0
list2=`ls $2/*_counts.bed`
for l in $list2
do
#	output7=`head -c -12 $l`
	sbatch -n 1 -J $l --wrap="Rscript scripts/DAlChrE_SEACR_1.1.R --input=$l --output=$l"
	names="${names},${l}"
done
numtargets=`ls $2/*_counts.bed | wc -l`

finished="FALSE"
while [[ $finished == "FALSE" ]]
do
	numout=`ls $2/*.normalized.bed | wc -l`
	if [[ $numtargets == $numout ]]
	then
		finished="TRUE"
		echo "Generating combined PCA at $(date +%T)"
		sbatch -n 1 -J combined --wrap="Rscript scripts/DAlChrE_SEACR_comb_1.1.R --output=$2" ## FINISH THIS COMMAND
		names="${names},combined"
	else
		sleep 10
	fi
done

finished="FALSE"
while [[ $finished == "FALSE" ]]
do
	running=`squeue -n $names | wc -l`
	if [[ $running -gt 1 ]]
	then
		let count=$running-1
		echo "Waiting for $count jobs to finish at $(date +%T)..."
		sleep 30
	else
		finished="TRUE"
		echo "All jobs finished at $(date +%T)..."
	fi
done
rm $2/*IgGfilter_comb.relaxed.bed
echo "Done at $(date +%T)!"

