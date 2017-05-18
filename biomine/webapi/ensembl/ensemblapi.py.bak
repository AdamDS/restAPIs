#!/usr/bin/python
# author: Adam D Scott
# first created: 2015*11*23
# citation: McLaren et. al. (doi:10.1093/bioinformatics/btq330)
#	http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btq330
#####
#	Operation					Function				Works With
#####
### Inherited
#	endpoint	"http://rest.ensembl.org"
#	subset		"/vep/human/hgvs/"
#	action		

from biomine.webapi.webapi import webapi
import xml.etree.ElementTree as ET
import json
import time
from biomine.variant.mafvariant import mafvariant
from biomine.variant.vepvariant import vepvariant

class ensemblapi(webapi):
	endpoint = "http://grch37.rest.ensembl.org"
	species = "human"
	hgvsSubset = "/vep/" + species + "/hgvs/"
	regionSubset = "/vep/" + species + "/region/"
	sequenceSubset = "/sequence/region/" + species + "/"
	translationSubset = "/map/translation/"
	blosum = "Blosum62"
	csn = "CSN"	
	compara = "Conservation"
	exac = "exac"
	genesplicer = "GeneSplicer"
	maxentscan = "MaxEntScan"
	updown = "UpDownDistance"
	callback = "callback"
	canonical = "canonical"
	ccds = "ccds"
	dbnsfp = "dbNSFP"
	dbscsnv = "dbscSNV"
	domains = "domains"
	hgvs = "hgvs"
	mirna = "miRNA"
	numbers = "numbers"
	protein = "protein"
	refseq = "xref_refseq"
	def __init__(self,**kwargs):
		subset = kwargs.get("subset",'')
#optional defaults as given by http://rest.ensembl.org/documentation/info/vep_hgvs_get
		self.blosum = kwargs.get( ensemblapi.blosum , False )
		self.csn = kwargs.get( ensemblapi.csn , False )
		self.compara = kwargs.get( ensemblapi.compara , False )
		self.exac = kwargs.get( ensemblapi.exac , False )
		self.genesplicer = kwargs.get( ensemblapi.genesplicer , False )
		self.maxentscan = kwargs.get( ensemblapi.maxentscan , False )
		self.updown = kwargs.get( ensemblapi.updown , 5000 )
		self.callback = kwargs.get( ensemblapi.callback , "" )
		self.canonical = kwargs.get( ensemblapi.canonical , False )
		self.ccds = kwargs.get( ensemblapi.ccds , False )
		self.dbnsfp = kwargs.get( ensemblapi.dbnsfp , "" )
		self.dbscsnv = kwargs.get( ensemblapi.dbscsnv , False )
		self.domains = kwargs.get( ensemblapi.domains , False )
		self.hgvs = kwargs.get( ensemblapi.hgvs , False )
		self.mirna = kwargs.get( ensemblapi.mirna , False )
		self.numbers = kwargs.get( ensemblapi.numbers , False )
		self.protein = kwargs.get( ensemblapi.protein , False )
		self.refseq = kwargs.get( ensemblapi.refseq , False )
		if not subset:
			super(ensemblapi,self).__init__(ensemblapi.endpoint,ensemblapi.hgvsSubset)
		else:
			if ( subset == ensemblapi.hgvsSubset or \
				subset == ensemblapi.translationSubset or \
				subset == ensemblapi.regionSubset ):
				super(ensemblapi,self).__init__(ensemblapi.endpoint,subset)
			else:
				print "biomine ERROR: bad subset. webapi.subset initializing to variant association results"
				super(ensemblapi,self).__init__(ensemblapi.endpoint,ensemblapi.hgvsSubset)

	def doAllOptions( self , **kwargs ):
		"biomine::webapi::ensembl::ensemblapi::doAllOptions"
		self.blosum = True
		self.csn = True
		self.compara = True
		self.exac = True
		self.genesplicer = True
		self.maxentscan = True
		self.canonical = True
		self.ccds = True
		self.dbscsnv = True
		self.domains = True
		self.hgvs = True
		self.mirna = True
		self.numbers = True
		self.protein = True
		self.refseq = True
		self.doOptions( **kwargs )

	def getOptions( self ):
		return {	ensemblapi.blosum : int(self.blosum) ,
					ensemblapi.csn : int(self.csn) ,
					ensemblapi.compara : int(self.compara) ,
					ensemblapi.exac : int(self.exac) ,
					ensemblapi.genesplicer : int(self.genesplicer) ,
					ensemblapi.maxentscan : int(self.maxentscan) ,
					ensemblapi.updown : self.updown ,
					ensemblapi.callback : self.callback ,
					ensemblapi.canonical : int(self.canonical) ,
					ensemblapi.ccds : int(self.ccds) ,
					ensemblapi.dbnsfp : self.dbnsfp ,
					ensemblapi.dbscsnv : int(self.dbscsnv) ,
					ensemblapi.domains : int(self.domains) ,
					ensemblapi.hgvs : int(self.hgvs) ,
					ensemblapi.mirna : int(self.mirna) ,
					ensemblapi.numbers : int(self.numbers) ,
					ensemblapi.protein : int(self.protein) ,
					ensemblapi.refseq : int(self.refseq) ,
		}
	def getOptionsText( self ):
		return {	ensemblapi.blosum : "blosum62" , #unsure
					ensemblapi.csn : "csn" ,
					ensemblapi.compara : "conservation" ,
					ensemblapi.exac : "exac_maf" , #maybe?
					ensemblapi.genesplicer : "gene_splicer" , #unsure
					ensemblapi.maxentscan : "MaxEntScan" ,
					ensemblapi.updown : "UpDownDistance" ,
					ensemblapi.callback : "callback" ,
					ensemblapi.canonical : "canonical" ,
					ensemblapi.ccds : "ccds" ,
					ensemblapi.dbnsfp : "dbnsfp" , #unsure
					ensemblapi.dbscsnv : "dbscSNV" , #unsure
					ensemblapi.domains : "domains" ,
					ensemblapi.hgvs : [ "hgvsc" , "hgvsp" ] ,
					ensemblapi.mirna : "mirna" , #unsure
					ensemblapi.numbers : [ "exon" , "intron" ] ,
					ensemblapi.protein : "protein_id" ,
					ensemblapi.refseq : "refseq_transcript_ids" ,
		}
		
	def setSubset(self,subset):
		self.subset = subset
		self.action = ""
	def resetURL(self):
		self.action = ""

	def beginQuery(self):
		self.action = ""
	def doOptions( self , **kwargs ):
#		print "biomine::webapi::ensembl::ensemblapi::doOptions"
		toData = kwargs.get( 'data' , False )
		self.action = "?"
		options = self.getOptions()
		for option in options:
			if option == ensemblapi.updown or \
			option == ensemblapi.dbnsfp or \
			option == ensemblapi.callback:
				if toData:
					self.addData( option , options[option] )
				else:
					self.action += str(option) + "=" + str(options[option]) + "&"
			elif options[option]:
				if toData:
					self.addData( option , options[option] )
				else:
					self.action += str(option) + "=1&"

	def annotateHGVSScalar2Response( self , hgvsNotated , **kwargs ):
		out = kwargs.get( "content" , 'text/xml' )
		contentInURL = kwargs.get( "inURL" , True )
		self.doOptions()
		self.action = hgvsNotated + "?"
		if contentInURL:
			self.action += "content-type=" + out
			return self.submit()
		else:
			return self.submit( content = out )
	def annotateHGVSArray2Dict( self , hgvsNotatedArray , **kwargs ):
		out = kwargs.get("content",'')
		resultDict = {}
		for var in hgvsNotatedArray:
			self.beginQuery()
			hgvsNotated = var.strip()
			self.annotateHGVSScalar2Response( hgvsNotated , **kwargs )
			resultDict[hgvsNotated] = self.response.text
		return resultDict
	def annotateVariantsPost( self , variants , **kwargs ):
		
#		print "biomine::webapi::ensembl::ensemblapi::annotateVariantsPost"
		doAllOptions = kwargs.get( 'allOptions' , True )
		maxPost = kwargs.get( 'maxPost' , 150 ) #bc error 400 (bad request), 403 (banned), or 504 (gateway/proxy server timeout)
		#maxPost = 400 #https://github.com/ensembl/ensembl-rest/wiki/POST-Requests
		#maxPost = 1000 #http://rest.ensembl.org/documentation/info/vep_region_post
		lengthVariants = len(variants)
		annotatedVariants = {} #dict of vepvariants
		timeAtSubmit = 0
		for i in range(0,lengthVariants,maxPost):
			j = i + maxPost
			if lengthVariants < maxPost:
				j = lengthVariants
			subsetVariants = variants[i:j]
			formattedVariants = []
			nullValue = "."
			delim = " "
			needReferences = self.checkInsertionsReference( subsetVariants , nullValue=nullValue , delim=delim )
			self.fullReset()
			self.setSubset( ensemblapi.regionSubset )
			self.doAllOptions( data=doAllOptions )
			for var in subsetVariants:
				inputVariant = var.ensembl()
				print inputVariant
				formattedVariants.append( inputVariant )
				vepvar = vepvariant( inputVariant=inputVariant , parentVariant=var )
				annotatedVariants[inputVariant] = vepvar
 			#following examples from documentation
			self.addData( "variants" , formattedVariants )
			self.addHeader( "Accept" , "application/json" )
			self.addHeader( "Content-Type" , "application/json" )
			[ annotatedVariants , attempts , success , timeAtSubmit ] = self.trySubmit( timeAtSubmit , **kwargs )
			if not success:
				print "BioMine::webapi::ensembl::ensemblapi Error: cannot access desired XML fields/tags for variants " ,
				print "[" + str(i) + ":" + str(j) + "]"
		return annotatedVariants

	def trySubmit( self , lastSubmitTime , **kwargs ):
		annotatedVariants = {}
		attempts = 0
		keepGoing = True
		timeAtSubmit = lastSubmitTime
		timeSinceSubmit = time.time()
		waitTime = 0
		frequencyLimit = 1/15 #max seconds between requests: 15 requests/second
		for attempts in range(0,30):
			if timeSinceSubmit - timeAtSubmit > waitTime:
				attempts += 1
				timeAtSubmit = time.time()
				self.submit( post=True , **kwargs )
				if self.response.ok and self.response.text:
					root = self.response.json()
					keepGoing = False
					for rootElement in root:
						var = vepvariant()
						var.parseEntryFromVEP( rootElement )
						var.setInputVariant()
						annotatedVariants[var.inputVariant] = var
					return [ annotatedVariants , attempts , True , timeAtSubmit ]
				else:
					try:
						waitTime = self.headers['X-RateLimit-Remaining']
					except:
						print( "BioMine::webapi::ensembl::ensemblapi::trySubmit Error: failed to query due to " + str( self.response.status_code ) + ": " + self.response.reason )
						return [ annotatedVariants , attempts , False , timeAtSubmit ]
		print( "BioMine::webapi::ensembl::ensemblapi::trySubmit Error: failed to query after " + str( attempts ) + " attempts. Stopping due to " + str( self.response.status_code ) + ": " + self.response.reason )
		return [ annotatedVariants , attempts , False , timeAtSubmit ]

	def checkInsertionsReference( self , variants , **kwargs ):
		self.setSubset( ensemblapi.sequenceSubset )
		needReferences = {}
		inputRegions = []
		for var in variants:
			if var.reference == "-":
				inputRegion = var.region()
				needReferences[var.genomicVar()] = var.region()
				inputRegions.append( inputRegion )
		if needReferences:
			self.addData( "regions" , inputRegions )
			inputRegions = []
			self.addHeader( "Accept" , "application/json" )
			self.addHeader( "Content-Type" , "application/json" )
			self.submit( post=True , **kwargs )
			if self.response.ok and self.response.text:
				needReferences = self.updateMissingReferences( variants , needReferences , **kwargs )
			else:
				print "References needed: " ,
				print needReferences
		return needReferences
	def updateMissingReferences( self , variants , needReferences , **kwargs):
		nullValue = kwargs.get( 'nullValue' , '.' )
		delim = kwargs.get( 'delim' , ' ' )
		for var in variants:
			genVar = var.genomicVar()
			if genVar in needReferences:
				root = self.response.json()
				for rootElement in root:
					ID = rootElement.get( 'id' )
					vals = ID.split( ':' )
					if vals[5] == "1":
						vals[5] = "+"
					elif vals[5] == "-1":
						vals[5] = "-"
					inputRegion = vals[2] + ":" \
								+ vals[3] + ".." \
								+ vals[4] + ":" \
								+ vals[5]
					#print var.region() + "\t----\t" + inputRegion
					if var.region() == inputRegion:
						vcfValues = [ vals[2] , vals[3] , nullValue ]
						refSeq = rootElement.get( 'seq' )
						#print refSeq + "  -->  " ,
						#vcfValues.append( refSeq[0] )
						vcfValues.append( nullValue )
						#print vcfValues ,
						#print "  -->  " ,
						vcfValues.append( refSeq[0] + var.alternate )
						#vcfValues.append( var.alternate )
						vcfValues.append( nullValue )
						vcfValues.append( nullValue )
						vcfValues.append( nullValue )
						#print vcfValues
						needReferences[genVar] = delim.join( vcfValues )
		return needReferences	
	
	def getRefSequence( self , variants , **kwargs ):
		for var in variants:
			var
	#def parseJSON( self , json ):
#		print "biomine::webapi::ensembl::ensemblapi::parseJSON - json: " ,
#		print json
		#if type( json ) == list:
#			print "json is a list"
			#for sub in json:
#				print sub
				#if type( sub ) == list or type( sub ) == dict:
					#self.parseJSON( sub )
		#elif type( json ) == dict:
#			print "json is a dict"
			#for sub in sorted(json.keys()):
#				print sub
				#if type( json.get( sub ) ) == list or type( json.get( sub ) ) == dict:
					#self.parseJSON( json.get( sub ) )

#if transcript not in transcripts dict, add new transcript annotation (transcript => [protein change , cdna ])

	def HGVSAnnotationHeader( self ):
		return "\t".join( [ "Gene" , "Mutation" , "Chromosome" , "Start" , "Stop" , "Reference" , "Variant" , "Strand" , "Mutation_Type" ] ) + "\n"
	def HGVSAnnotatedHeader( self , inputFile ):
		header = inputHeader( inputFile )
		return "\t".join( [ header.strip() , "Gene" , "Mutation" , "Chromosome" , "Start" , "Stop" , "Reference" , "Variant" , "Strand" , "Mutation_Type" ] ) + "\n"
	def HGVSErrorHeader( self ):
		return "\t".join( [ "Gene" , "Mutation" , "Error" ] ) + "\n"

	def annotateHGVSScalar2tsv( self , hgvsNotated , **kwargs ):
		out = "text/xml"
		head = kwargs.get( "header" , '' )
		line = kwargs.get( "line" , '' )
		contentInURL = kwargs.get( "inURL" , True )
		self.action = hgvsNotated + "?"
		#self.doOptions()
		if contentInURL:
			self.action += "content-type=" , out
			self.submit()
		else:
			self.submit( content = out )
		#print self.response.text
		geneVariant = hgvsNotated.split( ":" )
		return self.annotateHGVSScalarResponse2tsv( geneVariant , content = out , header = head , line = line )
	def annotateHGVSScalarResponse2tsv( self, geneVariant , **kwargs ):
		head = kwargs.get( "header" , True )
		line = kwargs.get( "line" , '' )
		annotations = ""
		errors = ""
		#print self.response.text
		if head:
			annotations += self.HGVSAnnotationHeader()
			errors += self.HGVSErrorHeader()
		root = self.getXMLRoot()
		print root
		try:
			for result in self.getElement( root , 'data' ):
				if not self.getEntry( result , 'error' ):
					allele = self.getEntry( result , 'allele_string' )
					alleles = ["" , ""]
					if allele:
						alleles = allele.split( '/' )
					start = self.getEntry( result , 'start' )
					if not start:
						start = ""
					stop = self.getEntry( result , 'end' )
					if not stop:
						stop = ""
					chromosome = self.getEntry( result , 'seq_region_name' )
					if not chromosome:
						chromosome = ""
					strand = self.getEntry( result , 'strand' )
					if not strand:
						strand = ""
					consequence = self.getEntry( result , 'most_severe_consequence' )
					if not consequence:
						consequence = ""
					annotations += "\t".join( [ line.strip() , geneVariant[0] , geneVariant[1].strip() , chromosome , start , stop , alleles[0] , alleles[1] , strand , consequence ] ) + "\n"
				else:
					#print response
					annotations = self.nullLine( 10 )
					errors += "\t".join( [ geneVariant[0] , geneVariant[1] , result.get('error') ] ) + "\n"
		except:
			print "biomine::webapi::ensembl::annotateHGVSScalarResponse2tsv Warning: no root= " ,
			print root
			self.errorCheck()
		return { "annotations" : annotations , "errors" : errors }
	def annotateHGVSArray2tsv( self , hgvsNotatedArray , **kwargs ):
		out = "text/xml"
		head = kwargs.get( "header" , '' )
		contentInURL = kwargs.get( "inURL" , True )
		annotations = ""
		errors = ""
		for var in hgvsNotatedArray:
			self.beginQuery();
			hgvsNotated = var.strip()
			print "working on: " + hgvsNotated
			self.annotateHGVSScalar2Response( hgvsNotated , content = out , **kwargs )
			results = self.annotateHGVSScalarResponse2tsv( hgvsNotated.split( ":" ) , content = out , header = head )
			head = False
			annotations += results["annotations"]
			errors += results["errors"]
		return { "annotations" : annotations , "errors" : errors }
	def annotateHGVSDict2tsv( self, resultDict , **kwargs ):
		head = kwargs.get( "header" , True )
		annotations = ""
		errors = ""
		if head:
			annotations += self.HGVSAnnotationHeader()
			errors += self.HGVSErrorHeader()
		for hgvsNotated , response in resultDict.iteritems():
			geneVariant = hgvsNotated.strip().split( ":" )
			output = self.annotateHGVSScalarResponse2tsv( geneVariant , response , header = False )
			annotations += output["annotations"]
			errors += output["errors"]
		return { "annotations" : annotations , "errors" : errors }

	#def annotateHGVSList( self , variants ): #ensembl hasn't made this possible...yet
	#	self.beginQuery();
	#	variantList = []
	#	for var in variants:
	#		variantList.append( var.strip() )
	#	json.dumps( variantList )
	#	self.addHeader( "Accept" , "application/json" )
	#	self.addData( "hgvs_notation" , variantList )
	#	self.doOptions()
	#	return self.submit( content = "text/xml" , data = variantList , post = True )
	#def annotateHGVSList2tsv( self , variantArray ):
	#	resultDict = self.annotateHGVSArray2Dict( variantArray , content = "text/xml" )
	#	return self.annotatedHGVSDict2tsv( resultDict , header = False )
	def annotateHGVSFile( self , inputFile , col1 , col2 , **kwargs ):
		output = kwargs.get( "output" , '' )
		if output:
			#print output
			fout = open( output , 'w' )
			fout.write( self.HGVSAnnotatedHeader( inputFile ) )
			inFile = open( inputFile , 'r' )
			next( inFile )
			for line in inFile:
				columns = line.split( "\t" )
				fields = [ columns[col1] , columns[col2] ]
				#print fields[0] + ":" + fields[1] + "\t\t" + line
				annotated = self.annotateHGVSScalar2tsv( ':'.join( fields ) , **kwargs )
				if annotated:
					fout.write( line.strip() + "\t" + annotated["annotations"].strip() + "\n" )
				else:
					fout.write( line.strip() + "\t" + self.nullLine( 10 ) + "\n" )
	
	def annotateHGVSArray2File( self , variantArray , outputFile , **kwargs ):
		output = self.annotateHGVSArray2tsv( variantArray , **kwargs )
		fout = open( outputFile , 'w' )
		fout.write( self.HGVSAnnotationHeader() )
		for annotation in output["annotations"]:
			fout.write( annotation )
		ferr = open( outputFile + ".err" , 'w' )
		ferr.write( self.HGVSErrorHeader() )
		for error in output["errors"]:
			ferr.write( error )

	def annotateHGVSindel( self , variant , **kwargs ):
		NotImplemented
		#calculate the start and stop c positions
		#for insertions (inframe):
			#get an example set of bases for the amino acid inserted
			#(could try each set and maybe take most damaging, if different outcomes)
		#for deletions (inframe):
			#get sequence of bases from start and stop using Sequence subset
		#for frameshifts:
			#insert/delete one and/or two random bases, assure change is not a stop codon
		#use the sequence and c positions to get VEP annotation
		#annotate each possiblility for a single indel
		#maybe take most damaging if different outcomes

	def annotateVariants( self , variants , **kwargs ):
		if variants:
			self.subset = ensemblapi.regionSubset
			for var in variants:
				if annotated:
					annotations.update( )

	def nullLine( self , columns ):
		nulls = ""
		for i in range(0,columns-1):
			nulls += "NULL\t"
		return nulls.rstrip()

	@classmethod
	def setSpecies( cls , species ):
		cls.species = species
		cls.hgvsSubset = "/vep/" + species + "/hgvs/"
		cls.regionSubset = "/vep/" + species + "/region/"
	def useGRCh38( cls ):
		cls.endpoint = "http://rest.ensembl.org"
	def useGRCh37( cls ):
		cls.endpoint = "http://grch37.rest.ensembl.org"

	@staticmethod
	def inputHeader( inputFile ):
		inFile = open( inputFile , 'r' )
		if inFile:
			line = next(inFile).decode()
			#print line
			return line
		else:
			return ""
