#! /usr/bin/env pyhton

import Annotation_debug as Annotation
import prediction_modules_debug as pm
import sys
import os

annotationFile = sys.argv[1]
aaFile = sys.argv[2]
organism = sys.argv[3]

aaHandle = open(aaFile)

geneFile = os.path.splitext(os.path.basename(annotationFile))[0] + ".genes"
geneHandle = open(geneFile,'w')

DRUGSdb = os.environ['DRUGS_PATH']

betaLactamProFile = DRUGSdb + "ExpectedProfiles.txt"

beta_lactams,bLacProfile = pm.import_profile(betaLactamProFile)

wildType = pm.import_wt(organism,DRUGSdb)
letterCode = {"Susceptible":"S","Intermediate":"I","Unknown":"U","Resistant":"R"}
antibioticLine = "Genome"
predictionLine = os.path.basename(annotationFile)

betaList = Annotation.import_list(DRUGSdb + "beta_lactamase_keywords_2.txt")
tetList = Annotation.import_list(DRUGSdb + "tetracycline_names.txt")
aminoglyList = Annotation.import_list(DRUGSdb + "aminoglycoside_names.txt")
ciproList = Annotation.import_list(DRUGSdb + "ciprofloxacin_names.txt")
trsxList = Annotation.import_list(DRUGSdb + "trimethoprim_sulfamethoxazole_names_2.txt")
chlorList = Annotation.import_list(DRUGSdb + "chloramphenicol_names.txt")

listList = [betaList,tetList,aminoglyList,ciproList,trsxList,chlorList]
antibioticList = ["Beta-lactamase","Tetracycline","Aminoglycoside","Fluoroquinolone","Bactrim","Chloramphenicol"]

annotationDict,resfam_positive = Annotation.import_annotations_Mitch(annotationFile)

geneList = []

for annotation in annotationDict.keys():
	nameDesc = annotationDict[annotation].split("\t")
	newAnnotation = Annotation.Annotation(nameDesc[0],nameDesc[1],annotation)
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
		geneHandle.write("Beta-lactam\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
	if gene.antibiotic == "Fluoroquinolone":
#		cipro_names.write(gene.name + "\n")
		ciproResGenes.append(gene)
#		cipro_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Aminoglycoside":
#		aminoglycoside_names.write(gene.name + "\n")
		aminoGlyRes.append(gene)
		geneHandle.write("Aminoglycoside\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		aminoglycoside_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Tetracycline":
#		tetracycline_names.write(gene.name + "\n")
		tetracyclineResGenes.append(gene)
		geneHandle.write("Tetracycline\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		tetracycline_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Bactrim":
#		trsx_names.write(gene.name + "\n")
		bactrimResGenes.append(gene)
		geneHandle.write("Trimethoprim\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		trsx_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Chloramphenicol":
#		chloramphenicol_names.write(gene.name + "\n")
		chlorResGenes.append(gene)
		geneHandle.write("Chloramphenicol\t{}\t{}\t{}\n".format(gene.unique.rstrip("|"),gene.name,gene.description))
#		chloramphenicol_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Other":
#		other_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
		pass
bLacResistance,poorHits = pm.beta_lactam(betaLactamResGenes,DRUGSdb,bLacProfile,len(beta_lactams))

antibioticLine = antibioticLine + "\t" + "\t".join(beta_lactams)
try:
	predictionLine = predictionLine + "\t" + "\t".join(bLacResistance)
except TypeError:
	for thing in bLacResistance:
		predictionLine = predictionLine + "\t" + "\t".join(thing)

if (organism == "Klebsiella pneumoniae") or (organism == "Enterobacter cloacae") or (organism == "Enterobacter aerogenes"):
	predictionList = predictionLine.split("\t")
	predictionList[1] = "Resistant"
	predictionLine = "\t".join(predictionList)

#predictionLine = pm.porin_consideration(predictionLine,organism,aaFile,DRUGSdb)

cipPhenotype,cipDescription = pm.ciprofloxacin(ciproResGenes,wildType)

antibioticLine = antibioticLine + "\tCiprofloxacin"
predictionLine = predictionLine + "\t" + cipPhenotype

agPhenotype = "Susceptible"
if len(aminoGlyRes) > 0:
	agPhenotype = pm.aminoglycoside(aminoGlyRes,"gentamicin",DRUGSdb)
antibioticLine = antibioticLine + "\tGentamicin"
predictionLine = predictionLine + "\t" + agPhenotype

tetPhenotype = "Susceptible"
tetResfamGenes = pm.reduce_to_resfams(tetracyclineResGenes)
if len(tetResfamGenes) > 0:
	tetPhenotype = "Resistant"
tetNames = [x.name for x in tetResfamGenes]
antibioticLine = antibioticLine + "\tDoxycycline"
predictionLine = predictionLine + "\t" + tetPhenotype

trsxPhenotype = pm.bactrim(bactrimResGenes,annotationFile)
antibioticLine = antibioticLine + "\tTrimethoprim-sulfamethoxazole"
predictionLine = predictionLine + "\t" + trsxPhenotype

chlorPhenotype = "Susceptible"
chlorResfamGenes = pm.reduce_to_resfams(chlorResGenes)
if len(chlorResfamGenes) > 0:
	chlorPhenotype = "Resistant"
chlorNames = [x.name for x in chlorResfamGenes]
antibioticLine = antibioticLine + "\tChloramphenicol"
predictionLine = predictionLine + "\t" + chlorPhenotype

print(antibioticLine)
print(predictionLine)
