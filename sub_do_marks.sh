#!/bin/bash                         #-- what is the language of this shell
#                                  #-- Any line that starts with #$ is an instruction to SGE
#$ -S /bin/bash                     #-- the shell for the job
#$ -o log                        #-- output directory (fill in)
#$ -e log                        #-- error directory (fill in)
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
dir=$3

script=/wynton/group/reiter/lauren/autoSeurat/do_marks_seurat.R

if [[ $dir == "" ]];then dir=/wynton/group/reiter/lauren; fi

mkdir $basename

rsync -av -f"+ */" -f"- *" $dir/$basename/$version $basename/
cp -r $dir/$basename/$version/**/*.rds $basename/$version/
cp $dir/$basename/$para $basename/

Rscript $script $basename $version $para

rsync --exclude *.rds -azvrP --update $basename/$version $dir/$basename/

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
