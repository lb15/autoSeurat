#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o /wynton/group/reiter/lauren/log_autoS                        #-- output directory (fill in)
#$ -e /wynton/group/reiter/lauren/log_autoS                         #-- error directory (fill in)
#$ -cwd                            #-- tell the job that it should start in your working directory
#$ -r y                            #-- tell the system that if a job crashes, it should be restarted
#$ -j y                            #-- tell the system that the STDERR and STDOUT should be joined
#$ -l mem_free=50G
#$ -l scratch=50G
#$ -l h_rt=72:00:00
#$ -m ea                           #--email when done
#$ -M Lauren.Byrnes@ucsf.edu        #--email

module load CBI
module load r/4.0.3

if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

cd "$TMPDIR"


basename=$1
version=$2
para="$1"_"$2"_parameters.csv
do_marks=$3
script=/wynton/group/reiter/lauren/autoSeurat/seurat_processing.R

mkdir $basename

cp -r /wynton/group/reiter/lauren/$basename/outs/filtered_feature_bc_matrix $basename/
cp /wynton/group/reiter/lauren/$basename/$para $basename/

Rscript $script $basename $version $para $do_marks

rsync -azvrP --update $basename/$version /wynton/group/reiter/lauren/$basename/testseurat

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
