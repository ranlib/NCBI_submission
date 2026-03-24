#!/bin/bash
#
# run tostadas pipeline
# to submit fastq files and their metadata to NCBI
#

# get meta data
ln -s /home/dbest/data/Covid19/CovidSeq/WWPilot_2026-03-04_AW/SampleSheet.csv
TOSTADAS-Metadata-Parser.py SampleSheet.csv

# add fastq information
# note assumption:
# 1) ncbi-spuid-sra exactly matches Sample_ID
# 2) NoLaneSplitting = true
# 3) Paired-end sequencing
# 4) One FASTQ per read per sample
FASTQ_DIR=$HOME/Analysis/Covid19/CovidSeq/WWPilot_2026-03-04_AW/untrimmed
add_fastq_information.py -i Metadata-Output/reportable_WWPilot_2026-03-04_AW.csv \
			 -o reportable_WWPilot_2026-03-04_AW.xlsx \
			 -y config_list.yaml\
			 -s SampleSheet.csv \
			 -d $FASTQ_DIR

exit 0

# doesn't work
#nextflow run main.nf -profile test,docker --species virus

# that works
#nextflow run $HOME/Software/tostadas/main.nf -profile mpox,test,docker --workflow biosample_and_sra
##nextflow run $HOME/Software/tostadas/main.nf -profile nwss,test,docker --workflow biosample_and_sra

# nextflow run $HOME/Software/tostadas/main.nf \
# 	 -profile docker\
# 	 --workflow biosample_and_sra \
# 	 --species virus \
# 	 --submission \
# 	 --annotation \
# 	 --outdir ./dieter_test \
# 	 --meta_path ./dieter_test/metadata_file.xlsx
# 	 --submission_config ./dieter_test/submission_config.yaml

#	 --species virus \
nextflow run $HOME/Software/tostadas/main.nf \
	 -profile test,docker\
	 --workflow biosample_and_sra \
	 --annotation \
	 --outdir ./tests \
	 --meta_path ./tests/wastewater_biosample_template.xlsx

exit  0
