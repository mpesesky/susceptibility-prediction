#! /usr/bin/env pyhton

import Annotation
import prediction_modules as pm
import sys
import os

annotationFile = sys.argv[1]
aaFile = sys.argv[2]
organism = sys.argv[3]
#rationalHandle = open("CiproRationale.txt",'a')

aaHandle = open(aaFile)

DRUGSdb = os.environ['DRUGS_PATH']

betaLactamProFile = DRUGSdb + "ExpectedProfiles.txt"

beta_lactams,bLacProfile = pm.import_profile(betaLactamProFile)

wildType = pm.import_wt(organism,DRUGSdb)
letterCode = {"Susceptible":"S","Intermediate":"I","Unknown":"U","Resistant":"R"}
antibioticLine = "Genome"
predictionLine = os.path.basename(annotationFile)

betaList = Annotation.import_list(DRUGSdb + "beta_lactamase_keywords.txt")
tetList = Annotation.import_list(DRUGSdb + "tetracycline_keywords.txt")
aminoglyList = Annotation.import_list(DRUGSdb + "aminoglycoside_keywords.txt")
ciproList = Annotation.import_list(DRUGSdb + "ciprofloxacin_keywords.txt")
trsxList = Annotation.import_list(DRUGSdb + "trimethoprim_sulfamethoxazole_keywords.txt")
chlorList = Annotation.import_list(DRUGSdb + "chloramphenicol_keywords.txt")

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



for gene in geneList:
	if gene.antibiotic == "Beta-lactamase":
		betaLactamResGenes.append(gene)
	if gene.antibiotic == "Fluoroquinolone":
		ciproResGenes.append(gene)
	if gene.antibiotic == "Aminoglycoside":
		aminoGlyRes.append(gene)
	if gene.antibiotic == "Tetracycline":
		tetracyclineResGenes.append(gene)
	if gene.antibiotic == "Bactrim":
		bactrimResGenes.append(gene)
	if gene.antibiotic == "Chloramphenicol":
		chlorResGenes.append(gene)

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
#rationalHandle.write(annotationFile)
#rationalHandle.write("\n")
#rationalHandle.write(cipDescription)
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
