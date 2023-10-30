
#Jess Rhodes

from ftplib import FTP
import gzip
import re
import sys
import copy
import os
import subprocess
import operator

def find_closest_gene(geneDictkeys,hit,numberOfGenes,maxDist):
	
	#positions in the gtf file for chromosome, start, and end
	CHROM = 0
	START = 1
	END = 2

	#positions in the hits file for chrom and position
	CHROMH = 0
	POS = 1

	#position of distance in the new list
	DIST = 2

	#empty lists
	distHit = []
	returnList = []


	#loop through the keys passed to the function
	#for every key, calculate distance between it and hit
	for i in geneDictkeys:
		if i[CHROM] == hit[CHROMH]:
			distance = distance_gene(i,hit)
			distanceList = [hit,i,distance]
			distHit.append(distanceList)

	#sort the resulting list of keys by the distance between them and hit
	distHit.sort(key=distance_pos)

	# initialize counter
	count = 0 
	
	#create list of genes for return
	for i in range(0,len(distHit)):
		#if the count is less than the amount of closest genes per locus
		if count < numberOfGenes:
			#if the distance is less than or equal to the maximum distance that can be considered "close"
			if distHit[i][DIST] <= maxDist:
				#if two things have the same distance (eg. a locus is inside two genes)
				#they will not be counted as separate genes, and all will be appended
				#that's why the counter doesn't go up in this situation
				if distHit[i][DIST] == distHit[i+1][DIST]:
					returnList.append(distHit[i])
				#in every other case, append the gene + counter goes up
				else:
					returnList.append(distHit[i])
					count += 1
	return(returnList)


def distance_gene(gene,hit):
	#position of elements in "gene"
	CHROM = 0
	START = 1
	END = 2

	#position of elements in "hit"
	CHROMH = 0
	POS = 1

	#calculate the distance to the start and end of gene
	distanceStart = abs(gene[START]-hit[POS])
	distanceEnd = abs(gene[END]-hit[POS])
	#calculate the end of gene
	lengthGene = abs(gene[END]-gene[START])

	#if both the distance to the start and the end of the gene is less than the length of the gene
	#then the locus is inside the gene, and the distance is zero
	if distanceStart <= lengthGene and distanceEnd <= lengthGene:
		distanceStart = 0
		distanceEnd = 0

	#return the smallest distance
	distanceList = [distanceStart,distanceEnd]
	return(min(distanceList))

def distance_pos(distarray):
	#just a little function as a sort key
	#we can have a little sorting, as a treat
	DIST = 2
	return(distarray[DIST])

#blank dictionary and lists
gtfDict = {}
hitsList = []
hitsDict = {}
sortedHitsList = []
configDict = {}

#config variable
ftpAddress = "ftp.flybase.org"
gtfDirectory = "/genomes/Drosophila_melanogaster/dmel_r5.39_FB2011_07/gff"
gtfFile = "dmel-all-r5.39.gff.gz"
featureType = "gene"

prefixlist = ["DarkvLight"]
chroms = ["2L","2R","3L","3R","X"]

hitsFile = "./BestSnpPositions.ECage.2021.tsv"
hitsNum = 1
maxDist = 1
outputFile = "./TopGenes.ECage.2021.tsv"

chromlist = ["2L","2R","3L","3R","X"]

#if there doesn't exist a temp file for files created during this -
#if not os.path.isdir("tempDmelFiles"):
	#then create it
	#this isn't a option in the config file to prevent code injection
	#that's very unlikely, but better safe than sorry
#	subprocess.call(["mkdir","tempDmelFiles"],shell=True)
#	print("folder tempDmelFiles created")

#if it doesn't already exist, copy gtf file from flybank
ftp = FTP(ftpAddress)
ftp.login()
ftp.cwd(gtfDirectory)
#make the right call string from the gtfFile created
retrstring = "RETR "+gtfFile
ftp.retrbinary(retrstring,open("tempDmelFiles/dmel.gff.gz", 'wb').write)
print("temp file dmel.gff.gz created from flybase gff file")


#gtf file positions
CHROM = 0
FEATURE = 2
START = 3
END = 4
STRAND = 6
GENEINFO = 8

#hits dictionary positions
GENEID = 0
GENESYM = 1
GENEDIST = 2

#hits file positions
POSCHR = 0
PVALUE = 1

#hits file
CHR = 0
POS = 1

#output from closestGene.py positions
HITS = 0
GENEKEY = 1
DIST = 2

#read in the gtf file
with gzip.open("tempDmelFiles/dmel.gff.gz","rt") as file:
	for line in file:
		if line.startswith("##"):
			continue
		elements = line.strip("\n").split("\t")
		#select for elements that are the right feature type
		try:
			if elements[FEATURE] == featureType:
				#make key from chromosome, start, end
				key = (elements[CHROM],int(elements[START]),int(elements[END]))
				#split the gene information up
				geneelements = elements[GENEINFO].split(";")
				modgeneelements = []
				for i in geneelements:
					i = i.strip()
					#get gene id and gene symbol
					if re.match(r'ID=',i):
						i =  i.strip("ID=").strip(" \"").strip("\"").strip("")
						modgeneelements.append(i)
					elif re.match(r'Name=',i):
						i =  i.strip("Name=").strip(" \"").strip("\"").strip("")
						modgeneelements.append(i)
					elif re.match(r'Alias=',i):
						i = i.strip("Alias=").strip(" \"").strip("\"").strip("").split(",") 
						for j in i:
							if j.startswith("FB"):
								modgeneelements.append(j)
				value = modgeneelements
				#key and value - gene info
				gtfDict.update({key:value})
		except:
			continue



#parse hits file
with open(hitsFile,"r") as file:
	count = 0
	for line in file:
		elements = line.strip("\n").split("\t")
		if count>0:
			keyhit = (str(elements[CHR]),int(elements[POS]))
			#append to a list
			hitsList.append(keyhit)
		count+=1

		#go through the positions and chromosomes and find the closest genes
		for i in hitsList:
			#pass the list of genes,the single hit, 
			#the number of hits wanted, and the max distance to closest gene
			cg = find_closest_gene(gtfDict.keys(),i,hitsNum,maxDist)
			#make a blank list for closest genes + distance
			cgList = []
			for j in cg:
				key = j[GENEKEY]
				#make a copy of the info for the closest gene stored in gtfDict
				geneinfohit = copy.deepcopy(gtfDict[key])
				#append the distance onto the info
				#append the info + distance onto the list of closest genes
				cgList.append(geneinfohit)
			#have a dictionary with the hits and their closest genes
			hitsDict.update({i:cgList})

count = 0

with open(outputFile, "w") as wfile:
	wfile.write("CHR"+"\t"+"POS"+"\t"+"ID"+"\t"+"NAME"+"\t"+"ALIAS"+"\n")
	for i in hitsList:
		for j in hitsDict[i]:
			wfile.write(str(i[0])+"\t"+str(i[1]))
			wfile.write("\t"+str(j[0])+"\t"+str(j[1]))
			if len(j)>2:
				wfile.write("\t"+str(j[2]))
				for k in range(3,len(j)):
					wfile.write(";"+str(j[k]))
			wfile.write("\n")

print("output file created")
#cleanup = input("clean up temp file and folder? (y/n)")
#if cleanup.lower() == "y":
	#for a mac/linux machine
	#subprocess.call(["rm","-r","tempDmelFiles"],shell=True)
	
	#for a windows computer
	#subprocess.call(["rmdir", "/s", "tempDmelFiles"],shell=True)





