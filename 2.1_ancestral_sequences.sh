##Requirements: PRANK
#Argument 1: Output folder of script 1, including all maf alignments
#Argument 2: Tree file, in nwk format
#Argument 3: Output tag
#Argument 4: FASTA with ORF proteins

for f in $1/orfs/*.maf; 
do 
    echo $(basename $f)
    mkdir -p $1/prank
    prank -d=$f -showanc -showevents -prunetree -prunedata -F -once -t=$2 -o=$1/prank/$(basename $f)
done

[ -f $3.ancestors ] && rm $3.ancestors
[ -f $3.prot.ancestors ] && rm $3.prot.ancestors
[ -f $3.nucl.ancestors ] && rm $3.nucl.ancestors

echo -e 'orf_id\tsp\tev_age\tsyn_age\tgained\tgained_convergent\tlost\tdenovo\tseq' > $3.ancestors
echo -e 'orf_id\tsp\tev_age\tsyn_age\tbranch\tTIS\tmaxORF\tconv\tseq' > $3.prot.ancestors
for f in $1/prank/*best.anc.dnd; 
do 
    python3 2.2_parsing_ancestors.py $f $3 $4
done
