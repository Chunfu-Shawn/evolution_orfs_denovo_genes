import sys
import argparse
import os
import glob

__author__ = "Chunfu xiao"
__contributor__="..."
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Chunfu xiao"
__email__ = "chunfushawn@gmail.com"

# Script to aligned ORF against orthologous regions. Requirements: BLAST
# Arguments
parser = argparse.ArgumentParser(description="Arguments for this script")
parser.add_argument("--prot_dir", type=str, help="(Required) Directory including fasta with all protein sequence of ortholog regions (output of 1_extract_multiple_alignments.py)")
parser.add_argument("--prot_tar", type=str, help="(Required) Fasta including all translated sequences encoded by the ORFs")
args = parser.parse_args()

# Main
# assign each species to the lineage
lineages = ["human","hominoid","catarrhini","simiiformes","primates","primatomorpha","euarchontoglires","boreoeutheria","placentalia","mammal"]
order = {"hg38":"human"}
for sp in ("panTro5","panPan2","gorGor5","ponAbe2","nomLeu3"):
	order[sp] = "hominoid"
for sp in ("rheMac8","macFas5","macNem1","papAnu3","manLeu1","cerAty1","chlSab2","nasLar1","rhiRox1","rhiBie1","colAng1","HLpilTep1"):
	order[sp] = "catarrhini"
for sp in ("calJac3","aotNan1","saiBol1","cebCap1"):
	order[sp] = "simiiformes"
for sp in ("tarSyr2","otoGar3","micMur3","proCoq1"):
	order[sp] = "primates"
for sp in ("galVar1", "eulMac1","eulFla1"):
	order[sp] = "primatomorpha"	
for sp in ("tupChi1","jacJac1","micOch1","criGri1","mesAur1","perManBai1","mm10","HLmusCar1","HLmusPah1","rn6","HLmerUng1","nanGal1","HLcasCan1","HLdipOrd2","hetGla2","HLfukDam1","cavPor3","chiLan1","octDeg1","speTri2","HLmarMar1","oryCun2","ochPri3"): #rodents+primates
	order[sp] = "euarchontoglires"

cetartiodactyla =("vicPac2","HLcamFer2","HLcamBac1","HLcamDro1","HLturTru3","HLorcOrc1","HLdelLeu1","lipVex1","phyCat1","balAcu1","HLbalMys1","bosTau8","HLbosInd1","bisBis1","bosMut1","bubBub1","HLoviAri4","HLoviCan1","HLcapHir2","panHod1","HLodoVir1","HLcerEla1","susScr11")
perissodactyla = ("HLequCab3","equPrz1","HLequAsi1","cerSim1")
ferae = ("felCat8","HLaciJub1","panTig1","HLpanPar1","canFam3","HLlycPic1","musFur1","HLenhLut1","HLailFul1","ailMel1","ursMar1","odoRosDiv1","lepWed1","neoSch1","manPen1","HLmanJav1")
chiroptera = ("pteAle1","HLpteVam2","rouAeg1","HLrhiSin1","HLhipArm1","eptFus1","myoDav1","myoBra1","myoLuc2","HLminNat1","HLdesRot1")
eulipotyphla = ("eriEur2","sorAra2","conCri1")
boreoeutheria = cetartiodactyla + perissodactyla + ferae + chiroptera + eulipotyphla

for sp in boreoeutheria: #previous+euarchontoglires
	order[sp] = "boreoeutheria"
for sp in ("loxAfr3","triMan1","HLproCap2","chrAsi1","echTel2","eleEdw1","oryAfe1","dasNov3","HLchoHof2"): #atlantogenata +boro
	order[sp] = "placentalia"
for sp in ("monDom5","sarHar1","HLphaCin1","ornAna2"): #nonplacental+rest
	order[sp] = "mammal"

# Generate all truncated protein if start is a ATG
prot_ort = args.prot_dir + "/orfs/*.fa"
for file in glob.glob(prot_ort):
	s = ""
	species = "hg38"
	outs = open(file.replace(".fa",".trunc.fa"),"w+")
	for line in open(file):
		tru_seq = ""
		new = -1
		if ">" in line:
			s = line.split()[0].replace(">","")
			species = s.split("_")[2]
		else:
			for n,c in enumerate(str(line).rstrip("\n")):
				# if M 
				if (c == "M") and (new == -1) and (n <= 2): # in-frame window 6 nt downstream
					outs.write(">" + s + "_" + str(n) + "\n")
					new = 1
					tru_seq += c
				# if stop
				elif c == "*" and (new == 1):
					tru_seq += c
					new = -1
					outs.write(tru_seq + "\n")
					tru_seq = ""
				# if end of ortholog region
				elif n == len(str(line).rstrip("\n"))-1 and (new == 1):
					tru_seq += c
					new = -1
					outs.write(tru_seq + "\n")
				elif new == 1:
					tru_seq += c
	outs.close()

# Perform blastp ORF peptide against the ortholog peptide sequences
prot_ort = args.prot_dir + "/orfs/*.trunc.fa"
out = open(prot_ort.replace("*.trunc.fa","all.prot_spec.out"),"w+")
out.write("orf_id\tlineage\tfarthest species\tpident\tevalue\tqcovs\n")
for file in glob.glob(prot_ort):
	blastp = {}
	os.system("blastp -query " + args.prot_tar + " -subject " + file +
			" -subject_besthit -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs evalue bitscore\" >  " + 
			file.replace(".fa",".blastp.out"))
	# read blastp results
	for line in open(file.replace(".fa",".blastp.out")):
		# other species
		sseqid = line.split("\t")[1]
		if not sseqid in blastp:
			blastp[sseqid] = []
		pident = float(line.split("\t")[2])
		evalue = float(line.split("\t")[13])
		qcovs = float(line.split("\t")[12])
		blastp[sseqid].append(pident)
		blastp[sseqid].append(evalue)
		blastp[sseqid].append(qcovs)

	# Find farthest species and lineage
	far_sp = ""
	far_idx = 0
	for k,v in blastp.items():
		if (v[0] >= 80) and (v[2] >= 70):
			print("{0}: {1}".format(k,v))
			species = k.split("_")[2]
			i = lineages.index(order[species])
			if i >= far_idx:
				far_sp = k
				far_idx = i

	print("This protein is " + lineages[far_idx] + "-specific.")
	orf_id = far_sp.split("_")[0]
	pident = blastp[far_sp][0]
	evalue = blastp[far_sp][1]
	qcovs = blastp[far_sp][2]
	out.write(orf_id + "\t" + lineages[far_idx] + "\t" + str(far_sp) + "\t" + str(pident) + "\t" + str(evalue) + "\t" + str(qcovs) + "\n")

out.close()