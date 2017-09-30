#!/usr/bin/python
# author: Adam D Scott (amviot@gmail.com)
# first created: 2016*01*11

import sys
import getopt
from biomine.webapi.exac.exacapi import exacapi
import biomine.variant.variant

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-i <inputFile> -o <outputFile>\n"
	helpText += "-g \"genomic variant\" -t (boolean flag for tsv output)\n"
	inputFile = ""
	output = ""
	genVar = ""
	tsv = False
	try:
		opts, args = getopt.getopt( argv , "h:i:o:g:t" , ["input=" , "output=" , "genVar="] )
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
		elif opt in ( "-g" , "--genVar" ):
			genVar = arg
		elif opt in ( "-t" , "--tsv" ):
			tsv = True
	return { "input" : inputFile , "output" : output , "genVar" : genVar , "tsv" : tsv }
	
def checkConnection():
	exac = "http://exac.broadinstitute.org/variant/22-46615880-T-C"
	res = requests.get( exac )
	if res:
		print "have response"
	else:
		print res.status_code

def main( argv ):
	values = parseArgs( argv )
	inputFile = values["input"]
	outputFile = values["output"]
	genVar = values["genVar"]
	tsv = values["tsv"]

	results = ""
	exacInstance = exacapi()
	print genVar
	g = genVar.split( '-' )
	var = variant.variant( chromosome=g[0] , start=g[1] , reference=g[2] , mutant=g[3] )
	var.reference = 'd'

	aF = exacInstance.getAlleleFrequency( var )
	var.printVariant(' ')
	print "Allele Frequency = " + str(aF)

if __name__ == "__main__":
	main( sys.argv[1:] )
