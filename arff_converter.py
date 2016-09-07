#! /usr/bin/env python

def import_phenotype(phenotypeFile,genomes):
	"""Retreive the phenotypes to each drug from the phenotype file and place them in the appropriate genome objects"""

	phenotypeHandle = open(phenotypeFile)
	phenotypeOrder = []

	for line in phenotypeHandle:
		if line.startswith("Genome"):
			fields = line.rstrip("\n").split("\t")[1:]
			for antibiotic in fields:
				phenotypeOrder.append(antibiotic)
		else:
			fields = line.rstrip("\n").split("\t")
			genome = ".".join(fields[0].split(".")[:2])
			phenotypes = {}
			for phenotypeIndex in range(1,len(fields)):
				phenotypes[phenotypeOrder[phenotypeIndex-1]] = fields[phenotypeIndex]
			for organism in genomes:
				if genome == organism.genome:
					organism.phenotype = phenotypes
					break
	phenotypeHandle.close()

	return genomes


def collect_attributes(infile,desired):
	"""Create Genomes for each individual organism in the infile, and define the options for each of the 'desired' attributes"""

	options = []
	genomes = []
	inhandle = open(infile)
	currentGenome = ""
	category = ""
	attributes = ""

	for line in inhandle:
		category = line.split(":")[0]
		attributes = line.rstrip("\n").split("\t")[1:]
	
		if category == "Genome":
			currentGenome = Genome(attributes)
		elif category in desired:
			if attributes[0] == "":
				setattr(currentGenome,category,[])
			else:
				setattr(currentGenome,category,attributes)
				options.extend(attributes)

		if category == "Beta-lactamase":
			genomes.append(currentGenome)
		
	inhandle.close()

	options = list(set(options))

	return genomes,options

def zoi_arff(relation,options,genomes,antibiotics,antibioticClasses,outfile):
	"""Output the data in ARFF format."""
	outhandle = open(outfile,'w')

	outhandle.write("@relation \"ZOI_{}\"\n\n".format(relation))
	
	if "Species" in antibiotics:
		outhandle.write("@attribute\tspecies\t{E.aerogenes,E.cloacae,E.coli,K.pneumoniae}\n")
	
	for gene in options:
		gene_out = gene.split(":")[0]
		outhandle.write("@attribute\t\"" + gene + "\"\t{Present,Absent}\n")
	for antibiotic in antibiotics[1:]:
		outhandle.write("@attribute\t\"" + antibiotic + "\" NUMERIC\n")
	
	outhandle.write("\n@data\n")
	for isolate in genomes:
		isolate.solidify()
		dataLine = []
		callSpecies = 0
		if "Species" in antibioticClasses:
			callSpecies = 1
		for attribute in antibioticClasses:
			if attribute == "Species":
				dataLine.append(isolate.Species)
			else:
				present =  getattr(isolate,attribute)
				for gene in options:
					if gene in present:
						dataEntry = "Present"
					else:
						dataEntry = "Absent"
					dataLine.append(dataEntry)
		for attribute in antibiotics[1:]:
			try:
				dataEntry = isolate.phenotype[antibiotic]
			except AttributeError:
				print(isolate.genome)
			dataLine.append(dataEntry)
		outhandle.write(",".join(dataLine))
		outhandle.write("\n")
	outhandle.close()


def make_arff(relation,options,genomes,antibiotics,antibioticClasses,outfile):
	"""Output the data in ARFF format."""
	outhandle = open(outfile,'w')

	outhandle.write("@relation \"{}\"\n\n".format(relation))
	
	if "Species" in antibiotics:
		outhandle.write("@attribute\tspecies\t{E.aerogenes,E.cloacae,E.coli,K.pneumoniae}\n")

	for gene in options:
		gene_out = gene.split(":")[0]
		outhandle.write("@attribute\t\"" + gene + "\"\t{Present,Absent}\n")
#		combinedOptions.append(gene)

	for antibiotic in antibiotics[1:]:
		outhandle.write("@attribute\t\"" + antibiotic + "\"\t{Susceptible,Resistant}\n")
#		combinedOptions.append(antibiotic)

	outhandle.write("\n@data\n")

	for isolate in genomes:
		isolate.solidify()
		dataLine = []
		callSpecies = 0
		if "Species" in antibioticClasses:
			callSpecies = 1
		for attribute in antibioticClasses:
			if attribute == "Species":
				dataLine.append(isolate.Species)
			else:
				present =  getattr(isolate,attribute)
				for gene in options:
					if gene in present:
						dataEntry = "Present"
					else:
						dataEntry = "Absent"
					dataLine.append(dataEntry)
		for attribute in antibiotics[1:]:
			try:
				dataEntry = isolate.phenotype[antibiotic]
			except AttributeError:
				print(isolate.genome)
			dataLine.append(dataEntry)
		outhandle.write(",".join(dataLine))
		outhandle.write("\n")
	outhandle.close()



class Genome:
	def __init__(self,Name):
		self.genome = Name[0]
		self.Species = self.get_species(Name[0])
		self.Fluoroquinolone = []
		self.Bactrim = []
		self.Aminoglycoside = []
		self.Chloramphenicol = []
		self.Tetracycline = []
		self.Beta_lactamase = []

	def get_species(self,name):
		return name.split("_")[0]
	
	def solidify(self):
		genes = []
		genes.extend(self.Fluoroquinolone)
		genes.extend(self.Bactrim)
		genes.extend(self.Aminoglycoside)
		genes.extend(self.Chloramphenicol)
		genes.extend(self.Tetracycline)
		genes.extend(self.Beta_lactamase)
		self.genes = set(genes)

import argparse
parser = argparse.ArgumentParser (description="Change a tab delimited text file of resistance genes to arff format for Weka")
parser.add_argument("infile",metavar="<tabFile>",help="Input tab-delimited text file")
parser.add_argument("phenotypeFile",metavar="<phenotypeFile>",help="Tab-delimited file of the phenotypes for each genome")
parser.add_argument("-o",dest="outfile",default="ResGenes.arff",help="Desired name of output arff file")
parser.add_argument("-attributes",dest="attributes",default="Species,Ampicillin,Cefazolin,Cefotetan,Ceftazidime,Ceftriaxone,Cefepime,Meropenem,Ciprofloxacin,TR-SX,Gentamicin,Doxycycline,Chloramphenicol",help="Comma-separated list of the attributes desired in the arff file")
parser.add_argument("-z",dest="zoi",action='store_true',default=False,help="Input phenotype file is zones of inhibition, rather than \"Susceptible/Resistant\"")
args = parser.parse_args()

antibioticClass = {
	"Species":"Species",
	"Ampicillin":"Beta-lactamase",
	"Cefazolin":"Beta-lactamase",
	"Cefotetan":"Beta-lactamase",
	"Ceftazidime":"Beta-lactamase",
	"Ceftriaxone":"Beta-lactamase",
	"Cefepime":"Beta-lactamase",
	"Meropenem":"Beta-lactamase",
	"Ciprofloxacin":"Fluoroquinolone",
	"TR-SX":"Bactrim",
	"Gentamicin":"Aminoglycoside",
	"Doxycycline":"Tetracycline",
	"Chloramphenicol":"Chloramphenicol"
}

desiredPhen = args.attributes.split(",")

desiredGen = [antibioticClass[x] for x in desiredPhen]

genomes,options = collect_attributes(args.infile,desiredGen)

genomes = import_phenotype(args.phenotypeFile,genomes)

if args.zoi:
	zoi_arff(args.attributes,options,genomes,desiredPhen,desiredGen,args.outfile)
else:
	make_arff(args.attributes,options,genomes,desiredPhen,desiredGen,args.outfile)


