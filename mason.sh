#!/bin/bash

# I start with assigning the flags (user inputs):

result_id=result

while getopts f:g:t:l:m:p:i:b: flag
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
    esac
done


# I print them out to be sure it worked out:
echo "bases_before= $bases_before"
echo "fasta: $fasta";
echo "gff: $gff";

mkdir "./data/$result_id"
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
    # Now I create a list of all PNAs:
    grep -A 1 $target "$REF/full_transcripts_$FASTA" | \
	sed -E 's/^([A-Z]{46}).*/\1/' > "$REF/targetgene_startreg.fasta"   # select -30 to + 16 region
    # Now I run the python script which I wrote to design PNAs:
    echo $length
    echo $result_id
    python ./scripts/make_pnas.py $length $RES $bases_before 
else
    echo "PNA $pna_input put in"
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

