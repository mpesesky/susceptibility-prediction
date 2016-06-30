#! /usr/bin/env python

def import_wt(organism,DRUGSdb):
	"""Input a fasta file and output a dictionary of the contained genes."""
# Specifically this program is designed to determine the wild type amino acid sequences of
# a few resistance-relevant genes from a curated database.

# This block identifies which fasta file applies to the organism at hand. This will need to
# be expanded as new organisms are tested.
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

# This block puts the wt amino acid sequences into a dictionary with the gene names from
# the header as the keys
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
		if line.startswith("Name"): # This line will is the header of the table
			fields = line.rstrip("\n").split("\t")
			antibiotics = fields[1:]
		else: # The gene by gene entries in the table are added to a dictionary with the
# gene names as the keys
			fields = line.rstrip("\n").split("\t")
			resistanceGenes[fields[0]] = fields[1:]
	return (antibiotics,resistanceGenes)

def compareSeq(seq,reference,start,end):
	"""Identify mismatching letters between two strings"""
# The purpose of this function in context is to check a sequenced gene against the known
# wt sequence
	mutant = [] 
	for i in range((start-1),(end-1)):
		try:
			nuc = seq[i]
		except IndexError:
			return mutant
		if nuc != reference[i]:
			mutant.append(reference[i] + str(i+1) + seq[i])
	return mutant

def beta_lactam(betaLactamases,DRUGSdb,profile,antibioticsTested):
	"""Predict resistance or susceptibility for a set of beta-lactamases."""
	import Annotation_debug as Annotation
	from subprocess import Popen,PIPE

	BLASTdb = DRUGSdb + "beta-lactamaseNomenclature" # BLAST database of specific beta-
#lactmase names
	blastQuery = "" # Initialize the string that will hold the fasta-formated query 
# sequences, to be submitted on the command line (through subprocess)
	phenotypeList = [] # list that will hold the expected resistance profiles for the
# genes identified by the BLAST search
	hitGenes = [] # The sequence-specific names of the genes identified by the BLAST

#	aaSequences = {} 

	for bLac in betaLactamases: # betaLactamases is a list of Annotation objects 
		try: # error check, will kill the program if amino acid sequences were not 
			x = bLac.aaSeq # added in the main DRUGS script
		except AttributeError:
			print(bLac.name + "\t" + bLac.unique)
			exit()
		blastLine = ">{}\n{}\n".format(bLac.name,bLac.aaSeq) # Fasta format
		blastQuery = blastQuery + blastLine # Add this entry to the query
	
# This block performs the command line BLAST, and outputs in a specific tabilar format.
# See BLAST help for an explanation of the flags. A culling limit of 1 ensures only the
# best hit is given when the alignment covers all (or almost all) of the query sequence
	cmd = "blastp -db {} -outfmt '6 qseqid sseqid pident length qlen slen' -culling_limit 1".format(BLASTdb)
	BLAST = Popen(cmd,shell=True,stdin=PIPE,stdout=PIPE)
	hitString,err = BLAST.communicate(blastQuery)

	hits = hitString.split("\n") # List of individual output lines
	poorHits = {} # Will hold information on partial or low identity hits. Not currently
# used by DRUGS.

	for hit in hits[:-1]: # consider each BLAST output line individually
		fields = hit.split("\t") # Break up the table into its component parts
		try:
			query = fields[0] # The query sequence name
			matchLen = int(fields[3]) # The number of bases in the alignment
			queryLen = int(fields[4]) # The number of bases in the query sequence
			subLen = int(fields[5]) # The number of bases in the subject sequence
			pident = float(fields[2]) # The percent identity of the alignment
		except IndexError: # If the BLAST output for some reason does not have enough fields
			print(hit)
			exit()

		if (pident > 95.0) and (matchLen == queryLen == subLen): # if the best hit is a
# reasonably close match to the query
			try: # find the corresponding expected resistance profile
				hitGenes.append(query.split("|")[0])
				phenotypeList.append(profile[fields[1]])
			except KeyError: # if the database does not have a profile for that gene
				print("No profile for " + str(fields[1]))
		else:
			backUp = query.split("|")[0]
			try: # treats the closest hit like it was a high quality, unless it has not match
# in the database
				hitGenes.append(backUp)
				phenotypeList.append(profile[backUp])
			except KeyError:
				poorHits[backUp] = hit
	print("\n".join(poorHits.keys())) # print out the names of genes with partial matches and
# no expected profile to be examined by hand
	print(hits[-1]) # should be a blank line

	combinedList = []
	if len(phenotypeList) == 0: # If the BLAST results were blank
		fullySusceptible = ['Susceptible'] * antibioticsTested # Put "Susceptible" in all fields
		return fullySusceptible,poorHits,hitGenes
	elif len(phenotypeList) == 1: # If there was a single beta-lactamase resistance gene
		for index,phenotype in enumerate(phenotypeList[0]): # Make the output phenotype exactly
# match the expected phenotype for that single gene
			if phenotype.strip() == 'R':
				phenotypeList[0][index] = "Resistant"
			elif phenotype.strip() == 'S':
				phenotypeList[0][index] = "Susceptible"
			elif phenotype.strip() == 'V':
				phenotypeList[0][index] = "Unknown"
		return phenotypeList,poorHits,hitGenes
	else:
		for i in range(len(phenotypeList[0])): #To get the number of beta-lactamases matched
			phenotype = 'Susceptible'
			for resGene in phenotypeList: # For each beta-lactamase determine:
				if resGene[i] == 'R':		# If any of the resistance genes provide resistance
					phenotype = 'Resistant' # Put "Resistant" for that antibiotic in the output
			combinedList.append(phenotype)
		return combinedList,poorHits,hitGenes

def closeMatch(int1,int2,distance): # Not currently used
	return (abs(int1-int2) <= distance)

# Not currently used because it was unreliable in the initial DRUGS dataset (which had only 4
# porin deletions)
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

	phenotype = "Susceptible" # Assume the organism will be susceptible until evidence is dected otherwise
	rationale = "No Resistance Determinants Detected"

	gyrAFound = 0 # Flag determining if the gyrA gene could be found in the draft genome
	parCFound = 0 # Flag determining if the parC gene could be found in the draft genome
	qnrFound = 0 # Flag determining if a qnr gene could be found in the draft genome
	effluxFound = 0 # Flag determining if an oqxAB gene could be found in the draft genome
	mutations = [] # Initialize a list that will keep track of variations from the wild type
	numResistanceFactors = 0 # Initialize a variable that will ennumerate the resistance factors found

	for cipGene in cipRelevant: # cipRelevant is a list of Annotation objects
# First test checks if the found gyrA and wt gyrA start with the same amino acid (M) if not then the found
# gene is missing some of its front amino acids and the position-based definition of the QRDR will not be
# relavent. Subsequent elif blocks for parC and gyrB are the same.
		if ("gyrA" in cipGene.description) and (cipGene.aaSeq[:1] == wildTypeQRDR["gyrA"][:1]):
			gyrAFound = 1
			mutant = compareSeq(cipGene.aaSeq,wildTypeQRDR["gyrA"],67,107) # look for variations in the QRDR
			mutations.extend(["gyrA:"+x for x in mutant]) # Report any variations in mutation notation
			if (len(mutant) != 0): # If at least one mutation was found
				numResistanceFactors += len(mutant) # Count the number of distinct mutations
				rationale = "gyrA mutation: {}".format(" ".join(mutant)) # Rationale for prediction. Not
# output anymore, but still recorded.
		elif ("parC" in cipGene.description) and (cipGene.aaSeq[:1] == wildTypeQRDR["parC"][:1]):
			parCFound = 1
			mutant = compareSeq(cipGene.aaSeq,wildTypeQRDR["parC"],67,107)
			mutations.extend(["parC:"+x for x in mutant])
			if (len(mutant) != 0):
				numResistanceFactors += len(mutant)
				rationale = "parC mutation: {}".format(" ".join(mutant))
		elif ("gyrB" in cipGene.description) and (cipGene.aaSeq[:1] == wildTypeQRDR["gyrB"][:1]):
			mutant = compareSeq(cipGene.aaSeq,wildTypeQRDR["gyrB"],427,428)
			mutations.extend(["gyrB:"+x for x in mutant])
			if (len(mutant) != 0):
				numResistanceFactors += len(mutant)
				rationale = "gyrB mutation: {}".format(" ".join(mutant))
		elif "Qnr" in cipGene.description: # Each qnr gene found is considered equivalent to a QRDR mutation
# in its ability to provide resistance. This is an oversimplification, but clinical levels of resistance are
# only seen if multiple QRDR mutations, or a single QRDR mutation plus a qnr or efflux pump are seen.
			qnrFound = 1
			numResistanceFactors += 1
			rationale += "qnr found\n"
			print(cipGene.name)
			mutations.append("qnr")
		elif ("oqx" in cipGene.description) or ("qep" in cipGene.description):
			effluxFound = 1 # The HMMs for these are not selective enough, so I did not use them in predictions
	if numResistanceFactors >= 2:
		phenotype = "Resistant"
	if ((gyrAFound == 0) or (parCFound == 0)) and (phenotype != "Resistant"):	#If the strain has already been found to be resistant, it does not matter if
																				#We cannot identify the other housekeeping gene 
		phenotype = "Unknown"
		rationale = "parC or gyrA not found"
		mutations.append("Missing")
	if qnrFound:
		rationale = rationale + "\nqnr gene found"
	if effluxFound:	#This will not occur with the current version of resfams, because it does not detect either oqx or qep gene.
		rationale = rationale + "\nFluoroquinolone efflux pump found"
	
	return (phenotype,rationale,mutations)


def aminoglycoside(aminoGlyRes,antibiotic,DRUGSdb,blast):
	"""Predict resistance for a given aminoglycoside antibiotic using the CARD database."""
	
	amglyProfileFile = DRUGSdb + "AminoglycosideResistanceProfile.txt" # Expected resistance profile
	amglys,amglyProfile = import_profile(amglyProfileFile)

# The antibiotic variable specifies the aminoglycoside being tested. This try block just tests if that 
# aminoglycoside is included in the expected profile database
	try:
		antibioticLocation = amglys.index(antibiotic)
	except ValueError:
		return "Specified aminoglycoside antibiotic not in aminoglycoside database."

	applicableResGenes = [] # List of resistance genes in the organism that provide resistance to the
# specified antibiotic
	misMatchResGenes = [] # List of resistance genes in the organism that don't have an expected
# resistance profile
	phenotype = "S" # Prediction is susceptibile until evidence is found for resistance

	for gene in aminoGlyRes: # aminoGlyRes is a list of Annotation objects
		if blast:
			geneClass = gene.description.split("_")[0]
		else:
			geneClass = gene.name.split("|")[0]
		
		genePhenotype = "S" # Initialize variable that will contain the expected phenotype for this gene
# against the specified antibiotic
		try: # find the expected profile for the gene against the specified antibiotic
			genePhenotype = amglyProfile[geneClass][antibioticLocation]
		except KeyError: # If there is not an expected profile, check if the organism has already been
# predicted resistant
			print(geneClass)
			if phenotype != "R":
				phenotype = "U" # if it is susceptible, predict U(nknown) because of this gene
#				print(geneClass)
				misMatchResGenes.append(geneClass)
		if genePhenotype == "R": # If the gene is expected to provide resistance, predict Resistant
			phenotype = "R"
			applicableResGenes.append(geneClass)
	uniqueGenes = set(applicableResGenes) # Remove duplicate aminoglycoside resistance gene families
# This block only returns the prediction
	if phenotype == "R":
		return "Resistant"
	elif phenotype == "U":
		return "Unknown"
	else:
		return "Susceptible"

def tetracycline(tetRes,DRUGSdb):
	# No content because resistance is predicted when any tet resistance gene is found from resfams
	return 1

# Depreciated. Predicted resistance if 2 divergent DHFRs were detected
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

# This version determines which of 20 predetermined DHFR gene clusters each DHFR in the
# organism belongs to. If two or more clusters are represented in the organism then it
# is predicted to be resistant, otherwise susceptible. While conceptually similar to 
# the bactrim function, practically it removes the need for arbitrary cutoffs in
# coverage length and percent identity from a two gene BLAST. The bactrim function
# is more generally applicable, though, since bactrim_2 requires clustering DHFRs from
# the relavent organisms in advance.
def bactrim_2(bactrimResGenes,DRUGSdb):
	from subprocess import Popen,PIPE

	clusterFasta = DRUGSdb + "DHFR_clusters.faa" # Sequences of the clusters
	clusterFile = clusterFasta + ".clstr" # Connects headers from the above fasta with 
# the cluster number.

	clusterHandle = open(clusterFile)
	clusterHeaders = {} # will be used to store the relationship between the header and
# cluster numbers
	currentCluster = "" # stores the cluster number as we loop through the clusterFile
	query = "" # will hold the BLAST query fasta sequences
	totalClusters = [] # will hold the different cluster numbers found in the organism
	phenotype = "" # output

# This block populates a dictionary with the various headers as the keys and the
# corresponding cluster numbers as the values
	for line in clusterHandle:
		if line.rstrip("\n").endswith("*"):
			clusterHeaders[line.split(" ")[1].rstrip(".").lstrip(">")] = currentCluster
		elif line.startswith(">"):
			currentCluster = line.rstrip("\n").lstrip(">")
	clusterHandle.close()
	
# This block creates a query from all of the DHFR sequences identified in the organism
	for gene in bactrimResGenes:
		query = query + ">{}\n{}\n".format(gene.name,gene.aaSeq)
	
# These lines perfom the BLAST and record the output
	cmd = "blastp -query - -subject {} -outfmt 6 -culling_limit 1".format(clusterFasta)
	blast = Popen(cmd,shell=True,stdout=PIPE,stderr=PIPE,stdin=PIPE)
	(blastOut,err) = blast.communicate(query)
	
	blastHits = blastOut.split("\n") # Separates the output into individual lines

	for hit in blastHits[:-1]: # For each BLAST output line
		fields = hit.split("\t")
		subject = fields[1] # The cluster fasta header matching the query
		mlen = int(fields[3]) # The alignment length
		percIdent = float(fields[2]) # The alignment percent identity
		if (mlen > 120) and (percIdent >= 95): # Permissive parameters to check query
# and alignment quality
			totalClusters.append(clusterHeaders[subject]) # Record quality alignments
	totalClusters = set(totalClusters) # Remove redundant clusters

	if len(totalClusters) > 1:
		phenotype = "Resistant"
	else:
		phenotype = "Susceptible"
	return phenotype,totalClusters


def reduce_to_resfams(resGeneList):
	"""Remove Pfam and TIGRFAMs entries from a list of annotations."""

	newList = []
	for gene in resGeneList:
		if gene.resfam:
			newList.append(gene)
	return newList


	
