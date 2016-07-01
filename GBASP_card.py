#! /usr/bin/env python

import Annotation
import prediction_modules as pm
import argparse as ap
import os

parser = ap.ArgumentParser(description="Predict antibiotic resistance from genome, using an antibiotic resistance BLAST database.")
parser.add_argument("-blast",dest="blastFile",help="Result file of BLAST against resistance database, use format '6 qseqid sseqid length qlen slen qstart qend pident'")
parser.add_argument("-proteins",dest="aaFile",help="Fasta file of organism proteome")
parser.add_argument("-organism",dest="org",default="Escherichia coli",help="Genus and species of pathogen (defaults to Escherichia coli)")
parser.add_argument("-genes",dest="genes",action="store_true",default=False,help="Turn on output of resistance gene families")
parser.add_argument("-blac_profile",dest="blacProfile",default="ExpectedProfiles.txt",help="Change the database for expected beta-lactam resistance")

args = parser.parse_args()

aaHandle = open(args.aaFile)

if args.genes:
	geneFile = os.path.splitext(os.path.basename(args.blastFile))[0] + ".genes"
	geneHandle = open(geneFile,'w')

DRUGSdb = "ReferenceFiles/"

betaLactamProFile = DRUGSdb + args.blacProfile

beta_lactams,bLacProfile = pm.import_profile(betaLactamProFile)

wildType = pm.import_wt(args.org,DRUGSdb)
letterCode = {"Susceptible":"S","Intermediate":"I","Unknown":"U","Resistant":"R"}
antibioticLine = "Genome"
predictionLine = os.path.basename(args.blastFile)

betaList = Annotation.import_list(DRUGSdb + "CRD_lactamase_keywords_name.txt")
tetList = Annotation.import_list(DRUGSdb + "CRD_tetracycline_keywords_name.txt")
aminoglyList = Annotation.import_list(DRUGSdb + "CRD_aminoglycoside_keywords_name.txt")
ciproList = Annotation.import_list(DRUGSdb + "CRD_fluoroquinolone_keywords_name.txt")
trsxList = Annotation.import_list(DRUGSdb + "CRD_TRSX_keywords_name.txt")
chlorList = Annotation.import_list(DRUGSdb + "CRD_amphenicol_keywords_name.txt")

listList = [betaList,tetList,aminoglyList,ciproList,trsxList,chlorList]
antibioticList = ["Beta-lactamase","Tetracycline","Aminoglycoside","Fluoroquinolone","Bactrim","Chloramphenicol"]

annotationDict,resfam_positive = Annotation.import_annotations_blast(args.blastFile)

geneList = []

for annotation in annotationDict.keys():
	name = annotationDict[annotation]
	newAnnotation = Annotation.Annotation(name,name,annotation)
	newAnnotation.classify(listList,antibioticList)
	newAnnotation.resfam = resfam_positive[annotation]
	geneList.append(newAnnotation)

betaLactamResGenes = []
ciproResGenes = []
aminoGlyRes = []
tetracyclineResGenes = []
bactrimResGenes = []
chlorResGenes = []

takeNext = 0
matchGene = ""
for line in aaHandle:
	if line.startswith(">"):
		for gene in geneList:
			if gene.unique in line:
				takeNext = 1
				matchGene = gene
				break
	elif takeNext == 1:
		matchGene.aaSeq = line.rstrip("\n")
		takeNext = 0
		matchGene = ""

#beta_lactam_aa = open("AA.beta-lactam.txt",'a')
#cipro_aa = open("AA.cipro.txt",'a')
#aminoglycoside_aa = open("AA.aminoglycoside.txt",'a')
#tetracycline_aa = open("AA.tetracycline.txt",'a')
#trsx_aa = open("AA.bactrim.txt",'a')
#chloramphenicol_aa = open("AA.chloramphenicol.txt",'a')
#other_aa = open("AA.other.txt",'w')

for gene in geneList:
	if gene.antibiotic == "Beta-lactamase":
		betaLactamResGenes.append(gene)
#		beta_lactam_names.write(gene.name + "\n")
		if args.genes:
			geneHandle.write("Beta-lactam\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
	if gene.antibiotic == "Fluoroquinolone":
#		cipro_names.write(gene.name + "\n")
		ciproResGenes.append(gene)
#		cipro_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Aminoglycoside":
#		aminoglycoside_names.write(gene.name + "\n")
		aminoGlyRes.append(gene)
		if args.genes:
			geneHandle.write("Aminoglycoside\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		aminoglycoside_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Tetracycline":
#		tetracycline_names.write(gene.name + "\n")
		tetracyclineResGenes.append(gene)
		if args.genes:
			geneHandle.write("Tetracycline\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		tetracycline_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Bactrim":
#		trsx_names.write(gene.name + "\n")
		bactrimResGenes.append(gene)
		if args.genes:
			geneHandle.write("Trimethoprim\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		trsx_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Chloramphenicol":
#		chloramphenicol_names.write(gene.name + "\n")
		chlorResGenes.append(gene)
		if args.genes:
			geneHandle.write("Chloramphenicol\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		chloramphenicol_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Other":
#		other_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
		pass
bLacResistance,poorHits,hitGenes = pm.beta_lactam(betaLactamResGenes,DRUGSdb,bLacProfile,len(beta_lactams))

antibioticLine = antibioticLine + "\t" + "\t".join(beta_lactams)
try:
	predictionLine = predictionLine + "\t" + "\t".join(bLacResistance)
except TypeError:
	for thing in bLacResistance:
		predictionLine = predictionLine + "\t" + "\t".join(thing)

if (args.org == "Klebsiella pneumoniae") or (args.org == "Enterobacter cloacae") or (args.org == "Enterobacter aerogenes"):
	predictionList = predictionLine.split("\t")
	predictionList[1] = "Resistant"
	predictionLine = "\t".join(predictionList)

#predictionLine = pm.porin_consideration(predictionLine,organism,aaFile,DRUGSdb)

cipPhenotype,cipDescription,cipMutations = pm.ciprofloxacin(ciproResGenes,wildType)

antibioticLine = antibioticLine + "\tCiprofloxacin"
predictionLine = predictionLine + "\t" + cipPhenotype

#trsxPhenotype = pm.bactrim(bactrimResGenes,args.blastFile)
trsxPhenotype = "Susceptible"
if len(bactrimResGenes) > 0:
	trsxPhenotype = "Resistant"
antibioticLine = antibioticLine + "\tTrimethoprim-sulfamethoxazole"
predictionLine = predictionLine + "\t" + trsxPhenotype

agPhenotype = "Susceptible"
if len(aminoGlyRes) > 0:
	agPhenotype = pm.aminoglycoside(aminoGlyRes,"gentamicin",DRUGSdb,True)
antibioticLine = antibioticLine + "\tGentamicin"
predictionLine = predictionLine + "\t" + agPhenotype

tetPhenotype = "Susceptible"
#tetResfamGenes = pm.reduce_to_resfams(tetracyclineResGenes)
if len(tetracyclineResGenes) > 0:
	tetPhenotype = "Resistant"
tetNames = [x.name for x in tetracyclineResGenes]
antibioticLine = antibioticLine + "\tDoxycycline"
predictionLine = predictionLine + "\t" + tetPhenotype

chlorPhenotype = "Susceptible"
chlorResfamGenes = pm.reduce_to_resfams(chlorResGenes)
if len(chlorResfamGenes) > 0:
	chlorPhenotype = "Resistant"
chlorNames = [x.name for x in chlorResfamGenes]
antibioticLine = antibioticLine + "\tChloramphenicol"
predictionLine = predictionLine + "\t" + chlorPhenotype

print(antibioticLine)
print(predictionLine)
