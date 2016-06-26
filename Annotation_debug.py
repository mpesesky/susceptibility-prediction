#! /usr/bin/env python

#Mitchell Pesesky 2014

def import_list(keywordFile):
	"""Take in a file of a list of keywords and convert to a python list"""
	keywordHandle = open(keywordFile)
	keywordList = []
	for line in keywordHandle:
		keywordList.append(line.rstrip("\n"))
	keywordHandle.close()
	return keywordList

def list_check(someList,annotation):
	"""Check all elements of a list against an annotation for matches"""
	for keyWord in someList:
		if (keyWord.upper() == annotation.name.split("|")[0].upper()):# or (keyWord.upper() in annotation.description.upper()):
			return 1
	return 0 

def import_annotations_blast(annotationFile):
	annotationHandle = open(annotationFile)
	annotations = {}
	resistancePositive = {}

	for line in annotationHandle:
		fields = line.split("\t")
		currentGene = str(fields[0].split(":")[1].split("|")[0] + "|")
		annotations[currentGene] = fields[1]
		resistancePositive[currentGene] = 1
	annotationHandle.close()
	return annotations, resistancePositive

def import_annotations_Mitch(annotationFile):
	"""Read in the annotation file that Mitch's full genome assembly toolbox outputs, return a dictionary of gene_numbers to annotations."""

	annotationHandle = open(annotationFile)
	annotations = {}
	currentGene = ""
	bestAnnotation = 0.0
	resfam_positive = {}

	for line in annotationHandle:
		if line.startswith(">"):
			currentGene = str(line.split(":")[1].split("|")[0] + "|")
			bestAnnotation = 0.0
			resfam_positive[currentGene] = 0
		elif line.startswith("-"):
			fields = line.rstrip("\n").split("\t")
			bitScore = float(fields[-1])
			if "RF" in fields[0].split("|")[1]:
				resfam_positive[currentGene] = 1
			if bitScore > bestAnnotation:
				annotations[currentGene] = fields[0].lstrip("-") + "\t" + fields[1]
				bestAnnotation = bitScore
	annotationHandle.close()
	
	return annotations,resfam_positive

def import_annotations(tabFile):
	"""Read in a .tab annotation output from Resfams, return a dictionary of gene_numbers to annotations."""

	tabHandle = open(tabFile)
	annotations = {}
	currentGene = ""
	bestAnnotation = 0.0

	for line in tabHandle:
		fields = line.split("\t")
		name = fields[0]
		bitScore = fields[-1].rstrip("\t")
		if name != currentGene:
			currentGene = name
			bestAnnotation = bitScore
			annotations[name] = fields[-4] + "\t" + fields[-3]
		elif bitScore > bestAnnotation:
			bestAnnotation = bitScore
			annotations[name] = fields[-4] + "\t" + fields[-3]
	tabHandle.close()

	return annotations


class Annotation:
	"""Scripts and class variables for resistance gene annotations derived from Resfams."""

	def __init__(self,shortName,longName,geneNum):
		"""Tie the decription of the gene's function to it"""
		self.name = shortName
		self.description = longName
		self.unique = geneNum

	
	def classify(self,listList,nameList):
		"""If possible bin the resistance gene into an antibiotic specific list for prediction"""
		
		flagList = []
		for antibioticList in listList:
			flagList.append(list_check(antibioticList,self))

		if sum(flagList) > 1:
			self.antibiotic = "Other"
		elif sum(flagList) == 0:
			self.antibiotic = "None"
		else:
			for i in range(len(flagList)):
				if flagList[i] == 1:
					self.antibiotic = nameList[i]
					break


