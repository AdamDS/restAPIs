#!/usr/bin/python
# author: Adam D Scott (amviot@gmail.com)
# first created: 2015*09*28

import sys
import getopt
from biomine.webapi.entrez.pubchemapi import pubchemapi

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-i <inputFile> -o <outputFile>\n"
	helpText += "-c \"Compound\" -d \"Identifier\"\n"
	inputFile = ""
	output = ""
	compound = ""
	identifier = ""
	search = "name"
	try:
		opts, args = getopt.getopt( argv , "h:i:o:c:d:s:" , ["input=" , "output=" , "compound=" , "identifier=" , "search="] )
	except getopt.GetoptError:
		print "biomine ERROR: Command not recognized"
		print( helpText ) 
		sys.exit(2)
	if not opts:
		print "biomine ERROR: Expected flagged input"
		print( helpText ) 
		sys.exit(2)
	for opt, arg in opts:
		#print opt + " " + arg
		if opt in ( "-h" , "--help" ):
			print( helpText )
			sys.exit()
		elif opt in ( "-i" , "--input" ):
			inputFile = arg
		elif opt in ( "-o" , "--output" ):
			output = arg
		elif opt in ( "-c" , "--compound" ):
			compound = arg
		elif opt in ( "-d" , "--identifier" ):
			identifier = arg
		elif opt in ( "-s" , "--search" ):
			search = arg
	return { "input" : inputFile , "output" : output , "compound" : compound , "identifier" : identifier , "search" : search }
	
def checkConnection():
	pubchemInstance = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/synonyms/XML"
	res = requests.get( pubchemInstance )
	if res:
		print "have response"
	else:
		print res.status_code

def readCompounds( inputFile ):
	compounds = []
	if inputFile:
		inFile = open( inputFile , 'r' )
		for line in inFile:
			compounds.append( line.strip() )
	return compounds
	
def main( argv ):
	values = parseArgs( argv )
	inputFile = values["input"]
	outputFile = values["output"]
	compound = values["compound"]
	identifier = values["identifier"]
	search = values["search"]

	results = ""
	compounds = readCompounds( inputFile )
	pubchemInstance = pubchemapi()

	#print search
	if inputFile and outputFile:
		results = pubchemInstance.compoundsSynonyms2File( compounds , outputFile , search = search , identifier = identifier )
	elif inputFile and not outputFile:
		results = pubchemInstance.compoundsSynonyms( compounds , search = search , identifier = identifier )

	#print results

	if compound:
		results = pubchemInstance.compoundSynonyms( compound , search = search , identifier = identifier )

	print results

	#print pubchemInstance.headers
	#print pubchemInstance.data
	#print pubchemInstance.buildURL()
	#for key , value in response.iteritems():
	#	fout.write( key + "\t" + value )
	#if response:
	#	print response.text
	#else:
	#	print response.status_code
	#print pubchemInstance

if __name__ == "__main__":
	main( sys.argv[1:] )
