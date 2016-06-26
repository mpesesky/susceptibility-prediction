#! /usr/bin/env python

def import_wt(organism,DRUGSdb):
	"""Input a fasta file and output a dictionary of the contained genes."""
	wtFile = ""
	if organism == "Enterobacter aerogenes":
		wtFile = DRUGSdb + "aerogenes_gyrA_parC.faa"
	elif organism == "Escherichia coli":
		wtFile = DRUGSdb + "coli_gyrA_parC.faa"
	elif organism == "Enterobacter cloacae":
		wtFile = DRUGSdb + "cloacae_wt.faa"
	elif organism == "Klebsiella pneumoniae":
		wtFile = DRUGSdb + "pneumoniae_wt.faa"
	wtHandle = open(wtFile)

	geneDict = {}
	currentHeader = ""
	for line in wtHandle:
		if line.startswith(">"):
			currentHeader = line.split(" ")[0].lstrip(">")
		else:
			geneDict[currentHeader] = line.upper().rstrip("\n")
	return geneDict

def import_profile(proFile):
	"""Convert a expected phenotype profile file into a dictionary of lists."""

	proHandle = open(proFile)
	antibiotics = []
	resistanceGenes = {}

	for line in proHandle:
		if line.startswith("Name"):
			fields = line.rstrip("\n").split("\t")
			antibiotics = fields[1:]
		else:
			fields = line.rstrip("\n").split("\t")
			resistanceGenes[fields[0]] = fields[1:]
	return (antibiotics,resistanceGenes)

def compareSeq(seq,reference,start,end):
	"""Identify mismatching letters between two strings"""
	mutant = []
	for i in range((start-1),(end-1)):
		if seq[i] != reference[i]:
			mutant.append(reference[i] + str(i+1) + seq[i])
	return mutant

def beta_lactam(betaLactamases,DRUGSdb,profile,antibioticsTested):
	"""Predict resistance or susceptibility for a set of beta-lactamases."""
	import Annotation
	from subprocess import Popen,PIPE

	BLASTdb = DRUGSdb + "beta-lactamaseNomenclature"
	blastQuery = ""
	phenotypeList = []

	aaSequences = {}

	for bLac in betaLactamases:
		try:
			x = bLac.aaSeq
		except AttributeError:
			print(bLac.name + "\t" + bLac.unique)
			exit()
		blastLine = ">{}\n{}\n".format(bLac.name,bLac.aaSeq)
		blastQuery = blastQuery + blastLine
	
	cmd = "blastp -db {} -outfmt '6 qseqid sseqid pident length qlen slen' -culling_limit 1".format(BLASTdb)
	BLAST = Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE)
	hitString,err = BLAST.communicate(blastQuery)

	hits = hitString.split("\n")
	poorHits = {}

	for hit in hits[:-1]:
		fields = hit.split("\t")
		try:
			matchLen = int(fields[3])
			queryLen = int(fields[4])
			subLen = int(fields[5])
			pident = float(fields[2])
		except IndexError:
			print(hit)
			exit()

		if (matchLen == queryLen == subLen) and (pident > 95.0):
			try:
				phenotypeList.append(profile[fields[1]])
			except KeyError:
				print("No profile for " + str(fields[1]))
		else:
#			print("Short: " + hit)
			if (pident == 100.0) and (matchLen > 100):
				backUp = fields[1]
			else:
				backUp = fields[0].split("|")[0]
			try:
				phenotypeList.append(profile[backUp])
				print("Short: " + hit)
			except KeyError:
				poorHits[backUp] = hit
	
	combinedList = []
	if len(phenotypeList) == 0:
		fullySusceptible = ['Susceptible'] * antibioticsTested
		return fullySusceptible,poorHits
	elif len(phenotypeList) == 1:
		for index,phenotype in enumerate(phenotypeList[0]):
			if phenotype.strip() == 'R':
				phenotypeList[0][index] = "Resistant"
			elif phenotype.strip() == 'S':
				phenotypeList[0][index] = "Susceptible"
			elif phenotype.strip() == 'V':
				phenotypeList[0][index] = "Unknown"
		return phenotypeList,poorHits
	else:
		for i in range(len(phenotypeList[0])): #To get the number of beta-lactamases tested
			phenotype = 'Susceptible'
			for resGene in phenotypeList: #For each beta-lactamase determine:
				if resGene[i] == 'R':		#Do any of the resistance genes provide resistance
					phenotype = 'Resistant'
			combinedList.append(phenotype)
		return combinedList,poorHits

def closeMatch(int1,int2,distance):
	return (abs(int1-int2) <= distance)

def porin_consideration(phenotypes,organism,aaFile,DRUGSdb):
	"""Helper program to correct for missing porins after beta-lactamase-based phenotypes have been predicted"""
	from subprocess import Popen,PIPE

	porinFile = ""
	if organism == "Enterobacter cloacae":
		return phenotypes
	elif phenotypes[4] == "Susceptible": #Assumes beta-lactams are in the order Ampicillin	Cefazolin	Cefotetan	Ceftazidime	Ceftriaxone
		return phenotypes
	elif organism == "Escherichia coli":
		porinFile = DRUGSdb + "Omps.faa"
	else:
		porinFile = DRUGSdb + "OmpKs.faa"

	cmd = "blastp -query {} -subject {} -outfmt '6 qseqid sseqid pident length qlen slen' -culling_limit 1".format(porinFile,aaFile)
	BLAST = Popen(cmd,shell=True,stdout=PIPE)
	hitString,err = BLAST.communicate()

	hits = hitString.split("\n")
	queries = []
	for hit in hits:
		if hit == "":
			continue
		fields = hit.split("\t")
		query = fields[0]
		perc_ident = float(fields[2])
		matchLength = int(fields[3])
		queryLength = int(fields[4])
		subjectLength = int(fields[5])

		if (perc_ident >= 90.00) and (matchLength >= 100) and (query not in queries):
			queries.append(query)
	
	phenotypeList = phenotypes.split("\t")
	
	if len(queries) < 2:
		for i in range(1,len(phenotypeList)):
			phenotypeList[i] = "Resistant"
	phenotypes = "\t".join(phenotypeList)
	return phenotypes


def ciprofloxacin(cipRelevant,wildTypeQRDR):
	"""Predict susceptibility to the fluoroquinolones."""
	#wildTypeQRDR should be a dictionary with keys as gene names and gene sequences as the values.
	#I believe that only gyrA and parC need to be checked based on current understanding

	phenotype = "Susceptible"
	numResistanceFactors = 0
	rationale = ""

	gyrAFound = 0
	parCFound = 0
	qnrFound = 0
	effluxFound = 0

	for cipGene in cipRelevant:
		if ("gyrA" in cipGene.description) and (cipGene.aaSeq[:1] == wildTypeQRDR["gyrA"][:1]):
			gyrAFound = 1
			mutant = compareSeq(cipGene.aaSeq,wildTypeQRDR["gyrA"],67,107)
			if (len(mutant) != 0):
				numResistanceFactors += len(mutant)
				rationale += "gyrA mutation: {}\n".format(" ".join(mutant))
		elif ("parC" in cipGene.description) and (cipGene.aaSeq[:1] == wildTypeQRDR["parC"][:1]):
			parCFound = 1
			mutant = compareSeq(cipGene.aaSeq,wildTypeQRDR["parC"],67,107)
			if (len(mutant) != 0):
				numResistanceFactors += len(mutant)
				rationale += "parC mutation: {}\n".format(" ".join(mutant))
		elif ("gyrB" in cipGene.description) and (cipGene.aaSeq[:1] == wildTypeQRDR["gyrB"][:1]):
			mutant = compareSeq(cipGene.aaSeq,wildTypeQRDR["gyrB"],427,428)
			if (len(mutant) != 0):
				numResistanceFactors += len(mutant)
				rationale += "gyrB mutation: {}\n".format(" ".join(mutant))
		elif ("Qnr" in cipGene.description) and (qnrFound == 0):
			qnrFound = 1
			numResistanceFactors += 1
			rationale += "qnr found\n"
		elif ("oqx" in cipGene.description) or ("qep" in cipGene.description):
			effluxFound = 1
	if numResistanceFactors >= 2:
		phenotype = "Resistant"
	if ((gyrAFound == 0) or (parCFound == 0)) and (phenotype != "Resistant"):	#If the strain has already been found to be resistant, it does not matter if
																				#We cannot identify the other housekeeping gene 
		phenotype = "Unknown"
		rationale = "parC or gyrA not found"
#	if qnrFound:
#		rationale = rationale + "\nqnr gene found"
#	if effluxFound:	#This will not occur with the current version of resfams, because it does not detect either oqx or qep gene.
#		rationale = rationale + "\nFluoroquinolone efflux pump found"
	
	return (phenotype,rationale)


def aminoglycoside(aminoGlyRes,antibiotic,DRUGSdb):
	"""Predict resistance for a given aminoglycoside antibiotic using the CARD database."""
	
	amglyProfileFile = DRUGSdb + "AminoglycosideResistanceProfile.txt"
	amglys,amglyProfile = import_profile(amglyProfileFile)

	try:
		antibioticLocation = amglys.index(antibiotic)
	except ValueError:
		return "Specified aminoglycoside antibiotic not in aminoglycoside database."

	applicableResGenes = []
	misMatchResGenes = []
	phenotype = "S"

	for gene in aminoGlyRes:
		geneClass = gene.name.split("|")[0]
		genePhenotype = "S"
		try:
			genePhenotype = amglyProfile[geneClass][antibioticLocation]
		except KeyError:
			if phenotype != "R":
				phenotype = "U"
				print(geneClass)
				misMatchResGenes.append(geneClass)
		if genePhenotype == "R":
			phenotype = "R"
			applicableResGenes.append(geneClass)
	uniqueGenes = set(applicableResGenes)
	if phenotype == "R":
		return "Resistant".format(antibiotic," ".join(uniqueGenes))
	elif phenotype == "U":
		return "Unknown".format(antibiotic," ".join(misMatchResGenes))
	else:
		return "Susceptible".format(antibiotic)

def tetracycline(tetRes,DRUGSdb):
	
	return 1

def bactrim(bactrimResGenes,annotationFile):
	"""Predict trimethoprim resistance if organism has multiple different dihydrofolate reductases."""
	from subprocess import Popen,PIPE
	from os import remove
	import os.path

	dhfrList = []
	for gene in bactrimResGenes:
		if "Dihydrofolate reductase" in gene.description:
			dhfrList.append(gene.aaSeq)
	if len(dhfrList) < 2:
		return "Susceptible"
	else:
		query = dhfrList[0]
		subject = ""
		for seq in dhfrList[1:]:
			head = ">DHFR_" + str(dhfrList.index(seq))
			subject = subject + head + "\n" + seq + "\n"
		outfile = "DHFRBLAST_" + os.path.basename(annotationFile)
		outhandle = open(outfile,'w')
		outhandle.write(subject)
		outhandle.close()
		cmd = "blastp -query - -subject {} -outfmt 6".format(outfile)
		blast = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE,stdin=PIPE)
		(BlastOut,Err) = blast.communicate(query)
#		print(BlastOut)
#		print(Err)
		remove(outfile)

		BlastLines = BlastOut.split("\n")
		for line in BlastLines:
			fields = line.split("\t")
			try:
				x = fields[2]
			except IndexError:
				print(BlastOut)
				exit()
			if float(fields[2]) < 95.00:
				return "Resistant"
		return "Susceptible"

def reduce_to_resfams(resGeneList):
	"""Remove Pfam and TIGRFAMs entries from a list of annotations."""

	newList = []
	for gene in resGeneList:
		if gene.resfam:
			newList.append(gene)
	return newList


	
