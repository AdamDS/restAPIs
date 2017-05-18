#!/usr/bin/python
# author: Adam D Scott (amviot@gmail.com)
# first created: 2015*12*06

import sys
import getopt
from biomine.webapi.entrez.entrezapi import entrezapi

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-d \"database\" "
	helpText += "-q \"query\" "
	helpText += "(-o \"outputFormat\")\n"
	inputFile = ""
	output = ""
	database = ""
	query = ""
	try:
		opts, args = getopt.getopt( argv , "h:i:d:q:o:" , ["input=" , "database=" , "query=" , "output="] )
	except getopt.GetoptError:
		print "biomine ERROR: Command not recognized"
		print( helpText ) 
		sys.exit(2)
	if not opts:
		print "biomine ERROR: Expected flagged input"
		print( helpText ) 
		sys.exit(2)
	for opt, arg in opts:
		print opt + " " + arg
		if opt in ( "-h" , "--help" ):
			print( helpText )
			sys.exit()
		elif opt in ( "-i" , "--inputFile" ):
			inputFile = arg
		elif opt in ( "-d" , "--database" ):
			database = arg
		elif opt in ( "-q" , "--query" ):
			query = arg
		elif opt in ( "-o" , "--output" ):
			output = arg
	return { "input" : inputFile , "database" : database , "query" : query , "output" : output }
	
def checkConnection():
	entrezInstance = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=asthma[mesh]+AND+leukotrienes[mesh]+AND+2009[pdat]"
	summaryTest = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=20113659,20074456"
	res = requests.get( entrezInstance )
	if res:
		print "have response"
	else:
		print res.status_code

def readDBSearches( inputFile ):
	searches = []
	dbs = []
	if inputFile:
		inFile = open( inputFile , 'r' )
		for line in inFile:
			fields = line.split( "\t" )
			dbs.append( fields[0].strip() )
			searches.append( fields[1].strip() )
	return { "databases" : dbs , "searches" : searches }
	
def main( argv ):
	values = parseArgs( argv )
	inputFile = values["input"]
	outputFormat = values["output"]
	database = values["database"].lower()
	query = values["query"]

	results = ""
	queries = {}
	if inputFile:
		queries = readDBSearches( inputFile )
	else:
		queries = { "databases" : database , "searches" : [ query ] }
	entrezInstance = entrezapi()

	if inputFile and outputFormat:
		results = entrezInstance.searchPubMed( searches , outputFormat , search = search , identifier = identifier )
	elif inputFile and not outputFormat:
		results = entrezInstance.getPubMedDetails( searches , search = search , identifier = identifier )

	print results

	if query:
		print queries
		if database == "pubmed":
			results = entrezInstance.searchPubMed( queries["searches"] )
		if database == "clinvar":
			results = entrezInstance.searchClinVar( queries["searches"] )

	#results = entrezInstance.getClinicalSignificance( queries["searches"] )
	#entrezInstance.getClinicalSignificance( queries["searches"] )

	print results.text

if __name__ == "__main__":
	main( sys.argv[1:] )
