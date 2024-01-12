#Argument 1: BED with ORFs
#Argument 2: Folder with MAF alignments, separated by chrm (chrm1.maf, chrm2.maf...)
#Argument 3: Output folder name
#Argument 4: nwk alignment
#Argument 5: FASTA with ORF proteins

#Requirements: Python3, biopython, BLAST, PHAST (maf_parse), PRANK
echo "Calculating multiple alignments. For huge multiple alignments, each alignment can take several minutes and we recommend to parallelize the searches for multiple ORFs."
python3 1_extract_multiple_alignments.py \
    -b $1 -m $2 -o $3 -f yes

# Calculating ancestral sequences and estimnate intact ORF ancestrally
echo "Calculating ancestral sequences and estimnate intact ORF ancestrally"
bash 2.1_ancestral_sequences.sh $3 $4 $1.ancestors $5

# Running BLAST across orthologous regions to get conservation score
echo "Running BLAST across orthologous regions to get conservation score"
python3 3_sequence_conservation_mammals.py -i $3 -p $5

# Running BLAST across orthologous regions to find specifc proteins
echo "Running BLAST across orthologous regions to find specifc proteins"
python3 4_sequence_specific.py --prot_dir $3 --prot_tar $5
