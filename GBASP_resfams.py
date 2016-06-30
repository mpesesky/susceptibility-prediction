#! /usr/bin/env pyhton

import Annotation_debug as Annotation
import prediction_modules_debug as pm
import sys
import os

annotationFile = sys.argv[1] # The output of HMM-based annotations
aaFile = sys.argv[2] # The exact amino-acid sequences used to create the above annotation file (with the same headers)
organism = sys.argv[3] # The genus and species name of the bacterium being analyzed

aaHandle = open(aaFile)

DRUGSdb = os.environ['DRUGS_PATH'] # An environment variable that specifies the directory where all the DRUGS databases
# reside. Must be set in advance.

betaLactamProFile = DRUGSdb + "ExpectedProfiles.txt" # A file containing gene specific antibiotic resistances

beta_lactams,bLacProfile = pm.import_profile(betaLactamProFile) # Create variables to hold the names and gene specific 
# resistance profiles for all beta-lactam antibiotic contained in the database
bacterium = ".".join(os.path.basename(annotationFile).split(".")[:2]) # store the specific name of the bacterium being
# analyzed
gene_file = bacterium + ".resGenes.txt" # Output file that will contain the names of all genes considered in this 
# analysis
geneHandle = open(gene_file,'w') 
geneHandle.write("Genome:\t{}\n".format(bacterium))
geneDict = {} # This dictionary will ultimately contain the gene names granting resistance to eacha ntibiotic for
# output to <gene_file>

wildType = pm.import_wt(organism,DRUGSdb) # Identify and store any wt amino acid sequences that will be needed
#letterCode = {"Susceptible":"S","Intermediate":"I","Unknown":"U","Resistant":"R"} 

# These two lines initialize the prediction output. The antibioticLine is actually the headers for an output table
# with one row, which is the prediction line.
antibioticLine = "Genome" # header for first column of the output table (a.k.a. organism name)
predictionLine = os.path.basename(annotationFile) # Organism name for the output table

# The following lines each store, for the different antibiotic classes, words in the annotation that distinguish 
# each specific type of resistance gene
betaList = Annotation.import_list(DRUGSdb + "beta_lactamase_keywords.txt")
tetList = Annotation.import_list(DRUGSdb + "tetracycline_names.txt")
aminoglyList = Annotation.import_list(DRUGSdb + "aminoglycoside_names.txt")
ciproList = Annotation.import_list(DRUGSdb + "ciprofloxacin_names.txt")
trsxList = Annotation.import_list(DRUGSdb + "trimethoprim_sulfamethoxazole_names.txt")
chlorList = Annotation.import_list(DRUGSdb + "chloramphenicol_names.txt")

listList = [betaList,tetList,aminoglyList,ciproList,trsxList,chlorList] # Easy storage for the keywords
# This list corresponds to the antibiotic categories of the keyword list above
antibioticList = ["Beta-lactamase","Tetracycline","Aminoglycoside","Fluoroquinolone","Bactrim","Chloramphenicol"]

annotationDict,resfam_positive = Annotation.import_annotations_Mitch(annotationFile) # A function that stores the
# annotations from a HMM output file created specifically with my (Mitch's) version of the HMM annotator. If you
# need to be reading these comments then this function probably needs to be updated for your application. The 
# resfam_positive dictionary contains a boolean value for whether any annotations for a gene were from resfams.

geneList = [] # Initializes a list of all genes in the genome as Annotation objects

for annotation in annotationDict.keys(): # Cycling through all of the genes from the annotation file
	nameDesc = annotationDict[annotation].split("\t") # Returns the predicted name and description of the gene
	newAnnotation = Annotation.Annotation(nameDesc[0],nameDesc[1],annotation) # Creates a new Annotation object
# using the predicted gene name, predicted gene description, and gene ID
	newAnnotation.classify(listList,antibioticList) # Determine if the gene is predicted to give resistance to any
# antibiotic, stored in newAnnotation.antibiotic
	newAnnotation.resfam = resfam_positive[annotation] # Mark if the gene was annotated from the resfams database
	geneList.append(newAnnotation)

# The following lists will hold annotation objects representing genes predicted to give resistance to each of the
# tested antibiotic classes
betaLactamResGenes = []
ciproResGenes = []
aminoGlyRes = []
tetracyclineResGenes = []
bactrimResGenes = []
chlorResGenes = []

takeNext = 0 # When reading a FASTA file, indicates if the next sequence line should be captured
matchGene = "" # Variable to hold a reference to a particular Annotation object from one iteration to the next
for line in aaHandle: # This block runs through the amino acid fasta file and finds the amino acid sequences 
# for each Annotation object, storing them in Annotation.aaSeq
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

# In this block, the Annotation objects corresponding to each antibiotic class are added to the respective lists
for gene in geneList:
	if gene.antibiotic == "Beta-lactamase":
		betaLactamResGenes.append(gene)
#		beta_lactam_names.write(gene.name + "\n")
#		beta_lactam_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Fluoroquinolone":
#		cipro_names.write(gene.name + "\n")
		ciproResGenes.append(gene)
#		cipro_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if (gene.antibiotic == "Aminoglycoside") and (gene.resfam): # The second condition is applied because some
#genes that have the Aminoglycoside keywords are only found in pfams and probably don't have resistance function
#		aminoglycoside_names.write(gene.name + "\n")
		aminoGlyRes.append(gene)
#		aminoglycoside_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Tetracycline":
#		tetracycline_names.write(gene.name + "\n")
		tetracyclineResGenes.append(gene)
#		tetracycline_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Bactrim":
#		trsx_names.write(gene.name + "\n")
		bactrimResGenes.append(gene)
#		trsx_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
	if gene.antibiotic == "Chloramphenicol":
#		chloramphenicol_names.write(gene.name + "\n")
		chlorResGenes.append(gene)
#		chloramphenicol_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))
#	if gene.antibiotic == "Other":
#		other_aa.write(">{}\n{}\n".format(gene.name,gene.aaSeq))

# This function makes predictions for all members of the beta-lactam antibiotic class, provided they are a part of
# the ExpectedProfiles.txt table in the DRUGSdb directory
bLacResistance,poorHits,geneDict["Beta-lactamase"] = pm.beta_lactam(betaLactamResGenes,DRUGSdb,bLacProfile,len(beta_lactams))


antibioticLine = antibioticLine + "\t" + "\t".join(beta_lactams) # Adding to the headers to be used in the output
# table

try: # The predictions can be returned in two different formats, depending on whether all beta-lactam resistance
# genes were present in the ExpectedProfiles.txt table or not. Either way this block adds the predictions for
# each beta-lactam antibiotic to the output table
	predictionLine = predictionLine + "\t" + "\t".join(bLacResistance) 
except TypeError:
	for thing in bLacResistance:
		predictionLine = predictionLine + "\t" + "\t".join(thing)

# This might be a little heavy handed, but these three species are all considered to be resistant to Ampicillin
# so I added this in to guarantee those predictions. This did not affect much since any beta-lactamase would also
# give Ampicillin resistance. This should be taken out or more creatively applied to make this script applicable
# to other data sets.
if (organism == "Klebsiella pneumoniae") or (organism == "Enterobacter cloacae") or (organism == "Enterobacter aerogenes"):
	predictionList = predictionLine.split("\t")
	predictionList[1] = "Resistant"
	predictionLine = "\t".join(predictionList)

#predictionLine = pm.porin_consideration(predictionLine,organism,aaFile,DRUGSdb)

# This function predicts resistance to the fluoroquinolones (including ciprofloxacin)
cipPhenotype,cipDescription,geneDict["Fluoroquinolone"] = pm.ciprofloxacin(ciproResGenes,wildType)

antibioticLine = antibioticLine + "\tCiprofloxacin" # Adding "Ciprofloxacin" header to the output table
predictionLine = predictionLine + "\t" + cipPhenotype # Adding the Cip prediction to the output table

agPhenotype = "Susceptible"
if len(aminoGlyRes) > 0: # If there are 0 aminoglycoside resistance genes, then the aminoglycoside prediction will
# be susceptible
	agPhenotype = pm.aminoglycoside(aminoGlyRes,"gentamicin",DRUGSdb)
geneDict["Aminoglycoside"] = [x.name.split("|")[0] for x in aminoGlyRes] # Adding the aminoglycoside resistance
# gene names to the resistance gene name dictionary (geneDict)
antibioticLine = antibioticLine + "\tGentamicin"
predictionLine = predictionLine + "\t" + agPhenotype

tetPhenotype = "Susceptible"
tetResfamGenes = pm.reduce_to_resfams(tetracyclineResGenes) # removes genes not called from resfams.

# There is no function to predict tetracycline or chloramphenicol resistance, if there are any genes matching
# resfams present in the genome predicted to give tetracycline resistance then the prediction is "Resistant"
# otherwise "Susceptible"
if len(tetResfamGenes) > 0:
	tetPhenotype = "Resistant"
tetNames = [x.name for x in tetResfamGenes]
geneDict["Tetracycline"] = [x.split("|")[0] for x in tetNames]
antibioticLine = antibioticLine + "\tDoxycycline"
predictionLine = predictionLine + "\t" + tetPhenotype

# Predict trimethoprim-sulfamethoxazole resistance
trsxPhenotype,geneDict["Bactrim"] = pm.bactrim_2(bactrimResGenes,DRUGSdb)
antibioticLine = antibioticLine + "\tTrimethoprim-sulfamethoxazole"
predictionLine = predictionLine + "\t" + trsxPhenotype

chlorPhenotype = "Susceptible"
chlorResfamGenes = pm.reduce_to_resfams(chlorResGenes)
if len(chlorResfamGenes) > 0:
	chlorPhenotype = "Resistant"
chlorNames = [x.name for x in chlorResfamGenes]
geneDict["Chloramphenicol"] = [x.split("|")[0] for x in chlorNames]
antibioticLine = antibioticLine + "\tChloramphenicol"
predictionLine = predictionLine + "\t" + chlorPhenotype

# Write the antibiotic resistance gene names out to a file
for antibiotic in geneDict.keys():
	geneHandle.write("{}:\t{}\n".format(antibiotic,"\t".join(geneDict[antibiotic])))
geneHandle.close()

# Write the output table to stdout
print(antibioticLine)
print(predictionLine)
