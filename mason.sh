#!/bin/bash
set -e
#set -o pipefail
################################################################################
#                              MASON command line                              #
#                                                                              #
# MASON can be used to design or analyse peptide nucleic acid or other ASOs    #
# for off-target effects and other sequence-specific attributes.               #
#                                                                              #
################################################################################
################################################################################
################################################################################
#                                                                              #
#  Created by Jakob Jung                                                       #
#                                                                              #
################################################################################
################################################################################
################################################################################

################################################################################
# Help                                                                         #
################################################################################
Help()
{
   # Display Help
   echo "MASON can be used to design or analyse peptide nucleic acid or other ASOs"
   echo
   echo "Version:     v1.0.0"
   echo "About:       Developed by Jakob Jung @Barquistlab, HIRI"
   echo "Docs & Code: https://github.com/BarquistLab/mason_commandline"
   echo "Mail:        jakobjung@tutanota.com"
   echo
   echo "Usage:         sh mason.sh [OPTIONS] -f <fasta> -g <gff/gff3/gtf>  -m <nr_mismatches> "
   echo
   echo "Options:"
   echo "Required:      "
   echo "               -f  FASTA file of target organism"
   echo "                   - (PATH)"
   echo
   echo "               -g  GFF annotation file of target organism"
   echo "                   - (PATH)"
   echo
   echo "               -m  Number of allowed mismatches for off-target"
   echo "                   screening"
   echo "                   - (INTEGER)"
   echo
   echo "               -i  ID of result; determines name of directory where"
   echo "                   results are stored"
   echo "                   - (STRING)"
   echo
   echo "Required if you want to design ASO sequences for a scpecific gene:"
   echo "                -t target gene; Use the locus tag, as annotated in the"
   echo "                   GFF, column 9, e.g. [...];locus_tag=b1253;[...]"
   echo "                   - (STRING)"
   echo
   echo "                -l length of ASO sequence; we usually recommend"
   echo "                   10-13 mers"
   echo "                   - (INTEGER)"
   echo
   echo "Optional if you want to design ASO sequences for a scpecific gene:"
   echo "                 -b bases before the start codon, in which Sequences will be designed"
   echo "                    - (INTEGER)"
   echo
   echo "Required if you already have ASO sequences in FASTA format:"
   echo "                -p PNA/ASO input, FASTA file of sequences to be screened"
   echo
   echo "Help:"
   echo "                -h print this help menu"
   echo
   echo "Version:"
   echo "                -V print version"
   echo
}

################################################################################
################################################################################
# Main program                                                                 #
################################################################################
################################################################################
################################################################################
# Process the input options. Add options as needed.                            #
################################################################################
# I start with assigning the flags (user inputs):

result_id=result
version="v1.0.0"

while getopts f:g:t:l:m:p:i:b:hV flag
do
    case "${flag}" in
	f) fasta=${OPTARG};;
	g) gff=${OPTARG};;
	t) target=${OPTARG};;
	l) length=${OPTARG};;
	m) mismatches=${OPTARG};;
	i) result_id=${OPTARG};;
	p) pna_input=${OPTARG};;
	b) bases_before=${OPTARG};;
	h) Help
	   exit;;
	V) echo "$version"
	   exit;;
    esac
done


# I print them out to be sure it worked out:
echo "fasta: $fasta";
echo "gff: $gff";

if  [ -z "$fasta"] || [-z "$gff" ] || [-z "$mismatches"]; then
        echo 'ERROR: Missing required arguments for fasta (-f) or gff (-g) or mismatches (-m)' >&2
        exit 1
fi




mkdir -p "./data/$result_id"
RES="./data/$result_id"
REF="$RES/reference_sequences"
OUT="$RES/outputs"
GFF="$(basename -- $gff)"
FASTA="$(basename -- $fasta)"  # get base names of files wo paths

echo "$REF"
echo "$OUT"

mkdir -p $REF $OUT

scp "$gff" "$REF/$GFF"
scp "$fasta" "$REF/$FASTA"

# same for full regions (change to whole CDS and 30 nt upstream):
grep -P "\tCDS\t|\tsRNA\t|\tncRNA\t" $gff |\
        awk -F'\t' 'BEGIN { OFS="\t" } \
{if ($7=="-") {$5=$5+30} \
else { $4=$4-30} print $0}'| \
    sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*locus_tag=([^;]+).*)/\1\4\3/' | \
    sed -E 's/([^\t]*\t[^\t]*\t)([^\t]*)(.*gene=([^;]+).*)/\1\2;\4\3/' \
	> "$REF/full_transcripts_$GFF"

# I extract the fasta files from the gff using bedtools:
bedtools getfasta -s -fi $fasta -bed "$REF/full_transcripts_$GFF"  \
	 -name+ -fo "$REF/full_transcripts_$FASTA"

echo "starting if statement"
echo $pna_input

if [ -z  "$pna_input" ];
then
    echo "no PNA put in"
    # check whether user input is correct:
    if [ -z "$target"] || [-z "$length" ] ||  [-z "$bases_before" ] ; then
        echo '\nERROR: Do not specify -t, -l, -b if you submit fasta files of ASOS with -p;' >&2
	echo '       Also do not specify -p if you have specified a target gene with  -t, -l, -b' >&2
        exit 1
    fi
    
    # Now I create a list of all PNAs:
    grep -A 1 $target "$REF/full_transcripts_$FASTA" | \
	sed -E 's/^([A-Z]{46}).*/\1/' > "$REF/targetgene_startreg.fasta"   # select -30 to + 16 region
    # Now I run the python script which I wrote to design PNAs:
    echo $length
    echo $result_id
    python ./scripts/make_pnas.py $length $RES $bases_before 
else
    echo "PNA $pna_input put in"
     # check whether user input is correct:
    if  [ ! -z "$target"] ||  [ ! -z "$length" ] ||   [ ! -z "$bases_before" ]; then
        echo '\nERROR: Do not specify -t, -l, -b if you submit fasta files of ASOS with -p' >&2
	echo '       Also do not specify -p if you have specified a target gene with  -t, -l, -b' >&2
        exit 1
    fi
    python ./scripts/modify_PNAs.py $pna_input $RES
fi


#Now I run seqmap on start regions and whole transcriptome:

seqmap $mismatches "$REF/aso_targets.fasta" "$REF/full_transcripts_$FASTA" \
       "$OUT/offtargets_fulltranscripts.tab" /output_all_matches \
       /forward_strand /output_statistics /available_memory:5000 



# I use awk to determine the mismatch positions:
for NAME in $(ls $OUT/offtargets_*.tab)
do
    echo "$NAME"
    NEWNAME=${NAME%.tab}_sorted.tab
    head -1 $NAME | sed -E "s/(.*)/\\1\tmismatch_positions\tlongest_stretch/" |
	sed -E  's/^trans_id/locus_tag\tgene_name\tstrand/' > $NEWNAME
    echo "$NEWNAME"
    sed 1d $NAME |\
	sort -t "$(printf "\t")"  -k4 |\
	awk -F'\t' 'BEGIN {OFS="\t"; pos=0} 
	{
	    max=length($3)
	    mm="none"
	    stretch=0
	    longest_stretch=0
	     for(i=1; pos == 0 && i <= max; i++) 
	     {
		    v1=substr($3, i, 1) 
		    v2=substr($5, i, 1)
		    if(v1 != v2)
		     {
		       stretch=0
		       if(mm=="none") {mm=i} else {mm=mm ";" i}
		     }
		     else 
		     {
			stretch++
			if(stretch > longest_stretch) {longest_stretch=stretch}
		     }
	     } 
	     $7=mm
	     $8=longest_stretch
	     $2=$2-32
	     print $0
	}' |  sed -E 's/^([^;:]*)::/\1;\1::/'| \
	    sed -E 's/^([^;:]*);([^:]*)::[^\(]*\(([\+\-])\)/\1\t\2\t\3/' >> $NEWNAME
    rm $NAME
    
done

rm -rf $OUT/*_sorted_sorted.tab
python ./scripts/summarize_offtargets.py $OUT 

touch "$RES/$target"

echo "MASON finished"  

