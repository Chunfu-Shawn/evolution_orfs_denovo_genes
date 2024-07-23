#! /usr/bin/python
__author__="chunfuxiao"
__date__ ="$Jul 23, 2024"

import sys
import os
import pandas as pd
import os.path
import argparse

parser = argparse.ArgumentParser(description="Assign the node of origin based on protein similarities between de novo genes and out-group species")
parser.add_argument("--gene_file", help="Path to de novo gene file",
					default="/home/user/data3/rbase/denovo_tumor/denovo_genes/100_denovo_gene_list/100_denovo_genes.compiled.txt")
parser.add_argument("--blastp_dir", help="Path to amino acid fasta file",
					default="/home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status/outgroup_species")
parser.add_argument("--output_dir", help="Path to amino acid fasta file",
					default="/home/user/data3/rbase/denovo_tumor/denovo_genes/denovo_status")
args = parser.parse_args()

# assign each out-group species to the lineage
outgroup_sp = ["Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Pongo_abelii","Nomascus_leucogenys",
				   "Macaca_mulatta","Macaca_fascicularis","Callithrix_jacchus","Otolemur_garnettii",
				   "Oryctolagus_cuniculus","Mus_musculus","Canis_lupus_familiaris","Loxodonta_africana",
				   "Monodelphis_domestica"]
order = {}
for sp in ["Pan_troglodytes","Pan_paniscus","Gorilla_gorilla","Pongo_abelii","Nomascus_leucogenys"]:
	order[sp] = "hominoid"
for sp in ["Macaca_mulatta","Macaca_fascicularis"]:
	order[sp] = "catarrhini"
for sp in ["Callithrix_jacchus"]: # 普通绒
	order[sp] = "simiiformes"
for sp in ["Otolemur_garnettii"]: # 小耳大婴猴
	order[sp] = "primates"
for sp in ["Oryctolagus_cuniculus","Mus_musculus"]:# 穴兔，小鼠 rodents+primates
	order[sp] = "euarchontoglires"

ferae = ["Canis_lupus_familiaris"] # 家犬
boreoeutheria = ferae

for sp in boreoeutheria: #previous+euarchontoglires
	order[sp] = "boreoeutheria"
for sp in ["Loxodonta_africana"]: # 亚洲象 atlantogenata +boro
	order[sp] = "placentalia"
for sp in ["Monodelphis_domestica"]: # 负鼠 nonplacental+rest
	order[sp] = "mammal"

# main
files = os.listdir(args.blastp_dir)

# read blastp results
data_list = []
for sp in outgroup_sp:
	tmp = pd.read_csv(args.blastp_dir + "/denovo_pep.blastp."+ sp +"_peptides.results.real.txt", sep="\t",
				   names=["Gene ID", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
			  "qlen", "sstart", "send", "slen", "qcovs", "evalue", "bitscore"])
	tmp["Species"] = sp
	data_list.append(tmp)

# merge and count peptides
dataset = pd.concat(data_list,ignore_index = True)
count_pep = dataset.groupby(["Gene ID","Species"]).size().reset_index(name="Count")
gene_list = list(pd.read_csv(args.gene_file, header=0, sep="\t")["Gene ID"])

# assign origin node
gene_sp = []
origin = {}
for gene in gene_list:
	# from hominoid to mammal
	origin[gene] = "human"
	for sp in outgroup_sp:
		gene_sp.append([gene,sp])
		tmp_count_pep = count_pep[(count_pep["Species"]==sp) & (count_pep["Gene ID"]==gene)]
		if not tmp_count_pep["Count"].empty:
			if tmp_count_pep.iloc[0,2] >1: # have homolog
				origin[gene] = order[sp]
gene_origin_df = pd.DataFrame(list(origin.items()), columns=["Gene ID","Outgroup Homolog"])
gene_origin_df.to_csv(args.output_dir + "/outgroup_homolog.peptide_similarity.txt", sep="\t", index=False)

# save peptide count
gene_sp_df = pd.DataFrame(gene_sp, columns=["Gene ID","Species"])
gene_sp_count_df = pd.merge(gene_sp_df,count_pep, how="left", on=["Gene ID","Species"]).fillna(0)
gene_sp_count_df["Count"] = gene_sp_count_df["Count"].astype("int8")
gene_sp_count_df.to_csv(args.output_dir + "/all.similar_peptides.count.txt", sep="\t", index=False)



