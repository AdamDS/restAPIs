#!/usr/bin/python
# author: Adam D Scott
# first created: 2015*12*06
#####
#	Operation	Function						Returns
#	esearch		term queries					id(s)
#	esummary	id search						entry details
#	efetch		download file of id				file
#	elink		linked entries in db A & B		dictionary of query keys & entry value
#	epost		id search						dictionary of query keys & entry value
#####
# Test sites
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=Leu858Arg
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=16609
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gene&db=protein&id=194680922,50978626,28558982,9507199,6678417&linkname=protein_gene&cmd=neighbor_history
### Inherited
#	endpoint	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
#	subset		"/esearch.fcgi"
#	action		"?db=pubmed&term=astma&usehistory=y"

import xml.etree.ElementTree as ET
from WebAPI.restAPI import restAPI
import variant

class entrezAPI(restAPI):
	endpoint = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
	esearch = "esearch.fcgi?"
	esummary = "esummary.fcgi?"
	elink = "elink.fcgi?"
	epost = "epost.fcgi?"
	efetch = "efetch.fcgi?"
	pubmed = "PubMed"
	clinvar = "ClinVar"
	protein = "Protein"
	gene = "Gene"
	dbsnp = "dbSNP"
	dbvar = "dbVar"
	omim = "OMIM"
	defaultGroup = "grabBag"
	largestBatchSize = 500
	def __init__(self,**kwargs):
		subset = kwargs.get("subset",'')
		if not subset:
			super(entrezAPI,self).__init__(entrezAPI.endpoint,entrezAPI.esearch)
		else:
			if ( subset == entrezAPI.esearch or subset == entrezAPI.esummary or
				 subset == entrezAPI.elink or subset == entrezAPI.epost or 
				 subset == entrezAPI.efetch ):
				super(entrezAPI,self).__init__(entrezAPI.endpoint,subset)
			else:
				print "ADSERROR: bad subset. restAPI.subset initializing to variant association results"
				super(entrezAPI,self).__init__(entrezAPI.endpoint,entrezAPI.hgvs)
		self.database = kwargs.get("database",entrezAPI.clinvar)
		self.queries = { entrezAPI.defaultGroup : "" }
		self.uids = []
		#self.linkname = ""
		#self.databaseFrom = ""
		#self.command = ""
		self.query_key = ""
		self.web_env = ""
		self.retmax = ""
		self.retstart = ""
		self.rettype = ""
		self.retmode = ""
		self.usehistory = ""
		self.assembly = "GRCh37"

	def addQuery( self , term , **kwargs ):
		field = kwargs.get( "field" , "" )
		group = kwargs.get( "group" , entrezAPI.defaultGroup )
		query = ""
		#print "Adding query: " + term + " " + field + " " + group , 
		if group in self.queries: #the group already has a query going
			#print "\tGroup=" + group + " exists!"
			if group != entrezAPI.defaultGroup: #not the default group
				#print "\t\tCurrent query is=" + self.queries[group]
				query = self.queries[group] + "+AND+"
			else: #anything goes in the default group
				#print "\t\tCurrent query is=" + self.queries[group]
				query = self.queries[group] + "+OR+"
		#print "\tcurrrent query is: " + query + "."
		query += term
		if field:
			query += "[" + field + "]"
		#print "\tNew query for group=" + group + " is: " + query
		tempq = { group : query }
		self.queries.update( tempq )
	def buildQuery( self ):
		query = ""
		for group in self.queries:
			if self.queries[group]:
				query += "(" + self.queries[group] + ")+OR+"
		return query
	def addID( self , uid ):
		self.uids.append( uid )
	
	def resetQuery( self ):
		self.queries = {}
	def resetID( self ):
		self.uids = []

	def buildSearchAction( self ):
		query = self.buildQuery()
		self.action = '&'.join( [ "db=" + self.database , "term=" + query ] )
		self.actionReturnParameters( )
		return self.action
	def buildWebEnvAction( self ):
		self.action = '&'.join( [ "db=" + self.database , "WebEnv=" + self.web_env , 
		"query_key=" + self.query_key ] )
		self.usehistory = ""
		self.actionReturnParameters( )
		return self.action
	def buildSummaryAction( self , ids ):
		self.uid = ','.join( ids )
		self.action = '&'.join( [ self.database , self.uid ] )
		return self.action
	def actionReturnParameters( self ):
		if self.rettype:
			self.action += "&rettype=" + self.rettype
		if self.retmode:
			self.action += "&retmode=" + self.retmode
		if self.retstart:
			self.action += "&retstart=" + str(self.retstart)
		if self.retmax:
			self.action += "&retmax=" + str(self.retmax)
		if self.usehistory:
			self.action += "&usehistory=" + self.usehistory
		return self.action
	def setRetmax( self , total ):
		if total > entrezAPI.largestBatchSize:
			self.retmax = entrezAPI.largestBatchSize
		else:
			self.retmax = total
	
	def doBatch( self , batchSize ):
		self.usehistory = "y"
		action = self.buildSearchAction( )
		#print action
		url = self.buildURL()
		response = self.submit()
		root = self.getXMLroot()
		for webenv_generator in root.iter( "WebEnv" ):
			self.web_env = webenv_generator.text
 		for querykey_generator in root.iter( "QueryKey" ):
			self.query_key = querykey_generator.text
		totalRecords = 0
		for count_generator in root.iter( 'Count' ):
			if totalRecords == 0:
				totalRecords = count_generator.text
		self.setRetmax( totalRecords )
		#print "WebEnv = " + self.web_env
		#print "QueryKey = " + self.query_key
		#print "Return max = " + str(self.retmax)
		#print "Batch size = " + str(batchSize)
		self.subset = entrezAPI.esummary
		self.action = self.buildWebEnvAction()
		self.buildURL()
		summaryResponse = self.submit()
		variants = self.getClinVarVariantEntry()
		traits = self.getClinVarTraitEntry()
		#self.getClinVarClinicalEntry()

	def getClinVarVariantEntry( self ):
		print "\tgetClinVarVariantEntry"
		root = self.getXMLroot()
		entries = {}
		for DocumentSummary in root.findall( 'DocumentSummary' ):
			print "\t\tdocsum"
			uid = DocumentSummary.attrib["uid"]
			print uid
			var = variant.variant()
			for variation in DocumentSummary.findall( 'variation' ):
				for variation_xref in DocumentSummary.findall( 'variation_xref' ):
					if variation_xref.find( 'db_source' ) == entrezAPI.dbsnp:
						var.dbSNP = variation_xref.find( 'db_id' )
						print var.dbSNP
					if variation_xref.find( 'db_source' ) == entrezAPI.omim:
						var.omim = variation_xref.find( 'db_id' )
						print var.omim
				for assembly_set in variation.findall( 'assembly_set' ):
					assembly_name = assembly_set.find( 'assembly_name' ).text
					print var.assembly_name
					if assembly_name.text == self.assembly:
						var.chromosome = assembly_set.find( 'chr' ).text
						var.start = assembly_set.find( 'start' ).text
						var.stop = assembly_set.find( 'stop' ).text
						var.mutant = assembly_set.find( 'alt' ).text
						var.reference = assembly_set.find( 'ref' ).text
						#assembly_acc_ver = assembly_set.find( 'assembly_acc_ver' )
			gene = DocumentSummary.find( 'gene' )
			var.gene = gene.find( 'symbol' ).text
			var.strand = gene.find( 'strand' ).text
			entries[uid] = var
			var.printVariant()
		return entries

	def getClinVarTraitEntry( self ):
		print "\tgetClinVarTraitEntry"
		root = self.getXMLroot()
		entries = {}
		for DocumentSummary in root.findall( 'DocumentSummary' ):
			uid = DocumentSummary.attrib["uid"]
			var.reset()
			#print uid
			for trait in DocumentSummary.findall( 'trait' ):
				trait_name = trait.find( 'trait_name' ).text
				trait_xrefs = {}
				entries[uid] = { trait_name : trait_xrefs }
				for trait_xref in trait.findall( 'trait_xref' ):
					db_source = trait_xref.find( 'db_source' ).text
					db_id = trait_xref.find( 'db_id' ).text
					txr = {}
					if trait_name in trait_xrefs:
						txr = entries[trait_name]
					txr.update( { db_source : db_id } )
					trait_xrefs.update( { trait_name : txr } )
				entries.update( { uid : trait_xrefs } )
		return entries

	def getClinVarClinicalEntry( self ):
		root = self.getXMLroot()
		entries = {}
		for DocumentSummary in root.findall( 'DocumentSummary' ):
			uid = DocumentSummary.attrib["uid"]
			var.reset()
			#print uid
			for clinical_significance in DocumentSummary.findall( 'clinical_significance' ):
				description = clinical_significance.find( 'description' ).text
				review_status = clinical_significance.find( 'review_status' ).text
				entries[uid]["description"] = description
				entries[uid]["review_status"] = review_status
		return entries

	def searchPubMed( self , query ):
		self.subset = entrezAPI.esearch
		self.database = entrezAPI.pubmed
		self.action = self.buildSearchAction( query )
		print self.subset
		print self.database
		print self.action
		return self.submit( )
	def searchClinVar( self , query , **kwargs ):
		self.subset = entrezAPI.esearch
		self.database = entrezAPI.clinvar
		self.action = self.buildSearchAction( query )
		print self.subset
		print self.database
		print self.action
		return self.submit( )
	
	def getClinicalSignificance( self , queries ):
		self.subset = entrezAPI.esearch
		#search ClinVar
		self.database = entrezAPI.clinvar
		query = self.buildSearchAction( queries )
		self.submit() 
		print self.response.text
		#parse the Id field within IdList
		ids = []
		root = ET.fromstring( self.response.text )
		for thisID in root.iter( 'Id' ):
			print thisID
			ids.append( thisID )
		#string ids together to search esummary
		searchIDs = ','.join( ids )
		print searchIDs
		#parse each DocumentSummary for clinical_significance
		#parse clinical_signficance for description
	
	def getXMLroot( self ):
		#print self.response.text
		return ET.fromstring( self.response.text )
