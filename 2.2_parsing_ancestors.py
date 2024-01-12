#! /usr/bin/python
__author__="jruizor"
__date__ ="$Nov 10, 2021 12:30:23 AM$"
'''Check P-site translation of a set of ORFs (p-sites)
'''

import sys
import subprocess
import os
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import os.path

# assign each species to the lineage
order = {}
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

# main
tree = sys.argv[1] #prot tree best.anc.dnd
fastan = sys.argv[1].replace(".dnd",".fas")
name = tree.split("/")[-1].replace(".maf.best.anc.dnd","")
tag = sys.argv[2]
fastap_orf = sys.argv[3]

# generate ancestral protein
fastap = open(sys.argv[1].replace(".dnd",".prot.fas"),"w+")
fa = {}
fan = SeqIO.index(fastan, "fasta")
for n in fan:
	fa[n] = str(Seq(str(fan[n].seq).replace("-","")).translate(cds = False))
	fastap.write(">" + n + "\n" + fa[n] + "\n")
fastap.close()

fline=open(tree).readline().rstrip().split(")")
last_sp = "hg38"
sub = 0
ancestor ={}
for n,line in enumerate(fline):
	if sub != 0:
		sub = sub - 1
		if not last_sp in ancestor:
			ancestor[last_sp] = ["none",[]]
		ancestor[last_sp][1].append(line.split(":")[0])
	else:
		if last_sp != "hg38":
			if not last_sp in ancestor:
				ancestor[last_sp] = ["none",[]]
			ancestor[last_sp][0] = line.split(":")[0]
		try:
			last_sp = line.split(",")[1].split(":")[0].replace("(","")
		except:
			pass
	if n != 0:
		sub = sub + line.count("(")
# print species and ancestor tree
print(ancestor)

outp = open(tag + ".prot.ancestors","a+")
# outp.write("orf_id\tsp\tev_age\tsyn_age\tbranch\tTIS\tmaxORF\tconv\tseq\n")

outh = open(tag + ".ancestors","a+")

outn = open(tag + ".nucl.ancestors","a+")
max_st_cons = "human"
max_sy_cons = "human"
max_pt_cons = "human"
# ORF structure
ev = {}
max_op ={}
tis={}
# protein conservation
blastp={}
cons={}

# perform blastp ORF peptide against the ortholog peptide sequences
# os.system("blastp -query " + fastap_orf + " -subject " + sys.argv[1].replace(".dnd",".prot.fas") + 
# 		" -subject_besthit -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen qcovs evalue bitscore\" >  " + 
# 		fastap_orf.replace(".fa",".blastp.out"))
# # read blastp results
# for line in open(fastap_orf.replace(".fa",".blastp.out")):
# 	sseqid = line.split("\t")[1]
# 	if not sseqid in blastp:
# 		blastp[sseqid] = []
# 	pident = float(line.split("\t")[2])
# 	evalue = float(line.split("\t")[13])
# 	qcovs = float(line.split("\t")[12])
# 	blastp[sseqid].append(pident)
# 	blastp[sseqid].append(evalue)
# 	blastp[sseqid].append(qcovs)
# print(blastp)

# estimate the presence of intact ORF
for last_sp in ancestor:
	tis[last_sp] = "0"
	old_anc = ancestor[last_sp][0]
	if len(str(fa[old_anc]).replace("-","").replace("X","*")) >= 4:
		#Check if methionine is the first position (3), or displaced one or two downstream positions (2), or there is an NTG otherwise (1)
		if str(fa[old_anc]).replace("-","").replace("X","*")[0] == str(fa["hg38"]).replace("-","").replace("X","")[0]:
			tis[last_sp] = "3"
		elif (str(fa[old_anc]).replace("-","").replace("X","*")[1] == str(fa["hg38"]).replace("-","").replace("X","")[0]) or (str(fa[old_anc]).replace("-","").replace("X","*")[2] == str(fa["hg38"]).replace("-","").replace("X","")[0]):
			tis[last_sp] = "2"
		elif str(fan[old_anc].seq).replace("-","")[1:3].upper() == "TG":
			tis[last_sp] = "1"
		else:
			tis[last_sp] = "0"

	#Longest ORF in the stretch
	op = -1
	max_op[last_sp] = 0
	if len(str(fa[old_anc]).replace("-","").replace("X","*")) >= 4:
		for n2 in str(fa[old_anc]).replace("-","").replace("X","*"):
			if (n2 == "M") and (op == -1):
				op = 1
			elif n2 == "*":
				if op > max_op[last_sp]:
					max_op[last_sp] = op
				op = -1
			elif op != -1:
				op += 1
		if op > max_op[last_sp]:
			max_op[last_sp] = op	
		try:
			max_op[last_sp] = float(max_op[last_sp])/len(str(fa["hg38"]).replace("-","").replace("X","").replace("*",""))*100
		except:
			pass

	# >= 70% of the sequence did not contain stop condons truncating the ORF and a conserved ATG start condon
	if (max_op[last_sp] >= 70) and (int(tis[last_sp]) >= 2):
		max_st_cons = order[last_sp]
		ev[last_sp] = "FIXED"
	# stop condons truncating the ORF resulting in < 70% of the ORF
	else:
		ev[last_sp] = "ABSENT"

	# protein conservation from blastp
	# if old_anc in blastp:
	# 	# if identity >= 50 and evalue <= 1e-5 and query coverage per subject (alignment length/subject length) >= 70% and the intact ORF is present
	# 	if (blastp[old_anc][0] >= 50) and (blastp[old_anc][1] <= float(1e-5)) and (blastp[old_anc][2] >= 70) and ev[last_sp] == "FIXED":
	# 		max_pt_cons = order[last_sp]
	# 		cons[old_anc] = "CONSERVED"
	# 	else:
	# 		cons[old_anc] = "DISSIMILAR"
	# else:
	# 	blastp[old_anc] = []
	# 	cons[old_anc] = "DISSIMILAR"
	
	max_sy_cons = order[last_sp]

# ortholog beyond intact ORF
beyond = 0
# whether is a de novo ORF
denovo = 0
for last_sp in ancestor:
	old_anc = ancestor[last_sp][0]
	# if ancestral ortholog beyond intact ORF
	if beyond == 1:
		denovo = 1
	if (order[last_sp] == max_st_cons) or (max_st_cons == "human"):
		beyond = 1

	# Check if ORF is present in some intermediate branch
	count  = [0,0,0]
	# for each species in intermediate branch (sp2)
	for sp2 in ancestor[last_sp][1]:
		tis2 = "0"
		if len(str(fa[sp2]).replace("-","").replace("X","*")) >= 4:
			#Check if methionine is the first position (3), or displaced one or two downstream positions (2), or there is an NTG otherwise (1)
			if str(fa[sp2]).replace("-","").replace("X","*")[0] == str(fa["hg38"]).replace("-","").replace("X","")[0]:
				tis2 = "3"
			elif (str(fa[sp2]).replace("-","").replace("X","*")[1] == str(fa["hg38"]).replace("-","").replace("X","")[0]) or (str(fa[sp2]).replace("-","").replace("X","*")[2] == str(fa["hg38"]).replace("-","").replace("X","")[0]):
				tis2 = "2"
			elif str(fan[old_anc].seq).replace("-","")[1:3].upper() == "TG":
				tis2 = "1"
			else:
				tis2 = "0"

		#Longest ORF in the stretch
		op = -1
		max_op2 = 0
		for n2 in str(fa[sp2]).replace("-","").replace("X","*"):
			if (n2 == "M") and (op == -1):
				op = 1
			if n2 == "*":
				op = -1
			elif op != -1:
				op += 1
		if op > max_op2:
			max_op2 = op	
		try:
			max_op2 = float(max_op2)/len(str(fa["hg38"]).replace("-","").replace("X","").replace("*",""))*100
		except:
			pass

		if (max_op2 >= 70) and (int(tis2) >= 2):
			if ev[last_sp] == "ABSENT":
				ev[last_sp] = "GAINED"
				count[0] += 1
				if beyond == 1:
					count[1] += 1
		else:
			if ev[last_sp] == "FIXED":
				ev[last_sp] = "LOST"	
				count[2] += 1	
		
		# estimate the conservation of the ortholog peptide sequences

	#Write output
	l = name + "\t" + last_sp + "\t" + order[last_sp] + "\t" + order[last_sp] + "\t" + old_anc + "\t" + str(tis[last_sp]) + "\t" + str(max_op[last_sp]) + "\t" + ev[last_sp] + "\t" + str(fa[old_anc]).replace("-","").replace("X","*") + "\n"
	outp.write(l.replace("W\n","\n"))
	outn.write(">" + name + "--" + last_sp + "\n" + str(fan[old_anc].seq) + "\n")	


outh.write(name + "\thg38\t" + str(max_st_cons) + "\t" + str(max_sy_cons) + "\t" + str(count[0]) + "\t" + str(count[1]) + "\t" + str(count[2]) + "\t" + str(denovo) + "\t" + str(fa["hg38"]).replace("-","").replace("X","") + "\n")

outn.write(">" + name + "--hg38\n" + str(fan["hg38"].seq) + "\n")

outp.close()
outn.close()

# os.system("rm " + fastan.split(".best")[0] + "*")
exit(0)


