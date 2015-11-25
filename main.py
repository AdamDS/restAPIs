#!/usr/bin/python
# author: Adam D Scott (amviot@gmail.com)
# first created: 2015*09*28

import sys
import getopt
import requests
import json
import tempfile
from ensemblAPI import ensemblAPI
from requests.auth import HTTPDigestAuth
from xmlutils.xml2json import xml2json
import xml.etree.ElementTree as ET

def parseArgs( argv ):
	helpText = "python main.py" + " "
	helpText += "-i <inputFile> -o <outputFile>\n"
	inputFile = ""
	output = ""
	try:
		opts, args = getopt.getopt( argv , "h:i:o:" , ["input=" , "output="] )
	except getopt.GetoptError:
		print "ADSERROR: Command not recognized"
		print( helpText ) 
		sys.exit(2)
	if not opts:
		print "ADSERROR: Expected flagged input"
		print( helpText ) 
		sys.exit(2)
	for opt, arg in opts:
		print opt + " " + arg
		if opt in ( "-h" , "--help" ):
			print( helpText )
			sys.exit()
		elif opt in ( "-i" , "--input" ):
			inputFile = arg
		elif opt in ( "-o" , "--output" ):
			output = arg
	return { "input" : inputFile , "output" : output }
	
def checkConnection():
	ensemblInstance = "http://rest.ensembl.org/info/ping?content-type=application/json"
	res = requests.get( ensemblInstance )
	if res:
		print "have response"
	else:
		print res.status_code

def readMutations( inputFile ):
	variants = []
	if inputFile:
		inFile = open( inputFile , 'r' )
		for line in inFile:
			fields = line.split( '\t' )
			variants.append( fields[0] + ":" + fields[1] )
	return variants
	
def main( argv ):
	values = parseArgs( argv )
	inputFile = values["input"]
	outputFile = values["output"]

	variants = readMutations( inputFile )
	ensemblInstance = ensemblAPI()
	#ensemblInstance.annotateHGVSList( variants )
	#print ensemblInstance.headers
	#print ensemblInstance.data
	#print ensemblInstance.buildURL()

	ensemblInstance.fOutAnnotateHGVS( outputFile , variants )
	#responses = ensemblInstance.annotateVariants( variants , content = "text/xml" )
	#response = ensemblInstance.annotateVariant( "EGFR:p.L858R" )

	#for key , value in response.iteritems():
	#	fout.write( key + "\t" + value )
	#if response:
	#	print response.text
	#else:
	#	print response.status_code
	#print ensemblInstance

if __name__ == "__main__":
	main( sys.argv[1:] )
