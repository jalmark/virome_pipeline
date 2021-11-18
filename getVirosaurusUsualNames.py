#!/bin/python
# Jalmari Kettunen 27.11.2020
# This script connects accession ID and usual name in Virosaurus reference genomes.
# There are two inputs: samtools idxstats output file and the header file.
# Header file has FASTA headers of the reference genomes, each on its own line.
# Output of this script is a tab-delimited text file which has both Virosaurus usual name and last 3 columns of idxstats output.
# Remember to run this in correct directory!
# Modified to take also multimapping reads into account.

referenceNames = ['90','98','Complete']
for ref in referenceNames:
	print(ref)
	if ref == 'Complete':
		headerFile = open('/users/jk427965/virosaurusCompleteGenomes_humanViruses_headers.txt','r')
	elif ref == '90' or ref == '98':
		headerFile = open('/users/jk427965/virosaurus' + ref + '_humanViruses_headers.txt','r')
	idxstatsFile = open('SRR1073679_' + ref + '_idxstats.txt','r')
	outputFile = open('SRR1073679_' + ref + '_idxstatsUsualNames.txt','w')
	
	accessionDict = {}
	for line in idxstatsFile:
		pieces = line.rstrip().split('\t')
		accession = pieces[0].rstrip(';')
		accessionDict[accession] = pieces[1:]
	for line in headerFile:
		fields = line.lstrip('>').rstrip().split(';')
		candidate = fields[0]
		if candidate in accessionDict:
			if ref == 'Complete':
				usualName = fields[1].replace(' species=','')
			elif ref == '90' or ref == '98':
				usualName = fields[1].replace(' usual name=','')
			constructLine = usualName + '\t' + '\t'.join(accessionDict[candidate])
			outputFile.write(constructLine + '\n')
	headerFile.close()
	idxstatsFile.close()
	outputFile.close()
