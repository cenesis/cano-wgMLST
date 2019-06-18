#!/bin/bash

if [ $# != 3 ]; then
	echo -e "\nUsage:   ${0} input_contigs_dir output_dir number_of_threads\n"
	echo -e "Example: ${0} ./sample_run/contigFiles ./sample_run 12\n"
	exit 1
else
	echo "The input contigs dir is '${1}', the output dir is '${2}' and using ${3} threaders"
fi

#echo -e "0\t"`date +%s` > ${2}/running_time

echo -e "\nRunning 1.Contig-Annotator ...\n"
1.Contig-Annotator.pl -i ${1} -o ${2}/1.contigAnn -p ${3}
#echo -e "1\t"`date +%s` >> ${2}/running_time

echo -e "\nRunning 2.PGAdb-Builder ...\n"
2.PGAdb-Builder.pl -i ${2}/1.contigAnn -o ${2}/2.PGAdb -p ${3}
#echo -e "2\t"`date +%s` >> ${2}/running_time

echo -e "\nRunning 3.wgMLST-Profiler ...\n"
3.wgMLST-Profiler.pl -i ${1} -d ${2}/2.PGAdb -o ${2}/3.wgProfiles -p ${3}
#echo -e "3\t"`date +%s` >> ${2}/running_time

echo -e "\nRunning 4.Dendro-Plotter ...\n"
4.Dendro-Plotter.pl -i ${2}/3.wgProfiles -o ${2}/4.DendroPlot
#echo -e "4\t"`date +%s` >> ${2}/running_time

echo -e "\nRunning 5.Loci-Extractor ...\n"
5.Loci-Extractor.pl -i ${2}/4.DendroPlot -o ${2}/5.ExtractedLoci
#echo -e "5\t"`date +%s` >> ${2}/running_time

