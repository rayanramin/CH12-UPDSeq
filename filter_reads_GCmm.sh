#!/bin/bash
#### function modified on 06/22/2022 
#### R.S.
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "This function filters a bam file for sequencing reads without at least one mismatch at any G:C position"
   echo "! requires samtools !"
   echo 
   echo "Syntax: filter_reads_GCmm.sh [h] -i input -o output" 
   echo "options:"
   echo "h     Print this Help."
   echo "i     input.bam"
   echo "o     GCmmfiltered.bam; use "-" for stdout (sam format))"
   echo
   echo "example: filter_reads_GCmm.sh -i input.bam -o GCmmfiltered.bam"
}

############################################################
while getopts ":i:o:h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) # input file
         input=$OPTARG;;
      o) # output file
         output=$OPTARG;;
     \?) # incorrect option
         echo "Error: Invalid option"
         exit;;
   esac
done
############################################################
## get a random number to use as a temporary file name
rand=$RANDOM
############################################################
# Main Function                                            #
############################################################
if [[ $# -eq 0 ]] ; then
    echo 'Error: No argument supplied'
    echo 'use -h option for help'
    exit 1
fi

## test if input file exists
if [ -f $input ]; then
   echo "Filtering reads with at least one mismatch at any G:C position from $input"
else
    echo "Error: Input file not found"
    exit 1
fi

### test if $output argument is missing
if [ -z "$output" ]; then
    echo "Error: Output file not specified"
    exit 1
fi
## test is $output is equal to "-"
if [ $output == "-" ]; then
    # echo "Error: Output file not specified"
samtools view -H $input > Xheader_$rand.txt
samtools view -h $input | grep -v "NM:i:0" | grep -E "MD:Z:.*([0-9]{1,3}[CG][0-9]{1,3})" | cat Xheader_$rand.txt - |  samtools view -S 
rm Xheader_$rand.txt
fi


if [ -f $output ]; then
   echo "Output $output already exists"
   echo "Will overwrite $output"
else
   echo "filtered bam file will be written to $output"
fi

samtools view -H $input > Xheader_$rand.txt
samtools view -h $input | grep -v "NM:i:0" | grep -E "MD:Z:.*([0-9]{1,3}[CG][0-9]{1,3})" | cat Xheader_$rand.txt - |  samtools view -Sb  > $output
rm Xheader_$rand.txt

#### test if output file is created
if [ -f $output ]; then
    echo "Completed; G:C mismatch filtered bam file is written to $output"
else
    echo "Error: G:C mismatch filtered file not created"
    exit 1
fi
############################################################
