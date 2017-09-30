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
from biomine.webapi.webapi import webapi
from biomine.variant.variant import variant
from biomine.variant.mafvariant import mafvariant
from biomine.variant.clinvarvariant import clinvarvariant
import re

class entrezapi(webapi):
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
			super(entrezapi,self).__init__(entrezapi.endpoint,entrezapi.esearch)
		else:
			if ( subset == entrezapi.esearch or subset == entrezapi.esummary or
				 subset == entrezapi.elink or subset == entrezapi.epost or 
				 subset == entrezapi.efetch ):
				super(entrezapi,self).__init__(entrezapi.endpoint,subset)
			else:
		#		print "biomine ERROR: bad subset. restapi.subset initializing to variant association results"
				super(entrezapi,self).__init__(entrezapi.endpoint,entrezapi.hgvs)
		self.database = kwargs.get("database",entrezapi.clinvar)
		self.queries = { entrezapi.defaultGroup : "" }
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
		group = kwargs.get( "group" , entrezapi.defaultGroup )
		query = ""
		#print "Adding query: " + term + " " + field + " " + group , 
		if group in self.queries: #the group already has a query going
			#print "\tGroup=" + group + " exists!"
			if group != entrezapi.defaultGroup: #not the default group
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
	def prepQuery( self , userVariants ): #expects a list of variants
		self.resetQuery()
		for var in userVariants:
			thisGroup = var.uniqueVar()
			self.addQuery( str(var.gene) , field="gene" , group=thisGroup )
			self.addQuery( str(var.chromosome) , field="chr" , group=thisGroup )
			self.addQuery( str(var.start) + ":" + str(var.stop) , field="chrpos37" , group=thisGroup )
			self.addQuery( "human" , field="orgn" , group=thisGroup )
			#self.addQuery( var.variantClass , "vartype" )
			#self.addQuery( var.referencePeptide + var.positionPeptide + var.alternatePeptide , "Variant name" )
			#var.referencePeptide , var.positionPeptide , var.alternatePeptide
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
		try:
			self.action = '&'.join( [ "db=" + self.database , "WebEnv=" + self.web_env , 
			"query_key=" + self.query_key ] )
			self.usehistory = ""
			self.actionReturnParameters( )
			return self.action
		except:
			print "entrez Error: can't use WebEnv"
			pass
	def buildSummaryAction( self , ids ):
		self.uid = ','.join( ids )
		self.action = "db=" + self.database + "&id=" + self.uid
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
	def setRetmax( self , total , **kwargs ):
		summaryBatchSize = kwargs.get( 'summaryBatchSize' , entrezapi.largestBatchSize )
		if total > summaryBatchSize:
	#		print str(total) + " > " + str(summaryBatchSize)
			self.retmax = summaryBatchSize
		else:
	#		print str(total) + " <= " + str(summaryBatchSize)
			self.retmax = total
	
	def doBatch( self , summaryBatchSize ):
		self.usehistory = "y"
		action = self.buildSearchAction( )
		url = self.buildURL()
#		print url
		response = self.submit()
		root = self.getXMLRoot()
		variants = {}
		if root != "<None/>":
			self.web_env = self.getEntry( root , 'WebEnv' )
			self.query_key = self.getEntry( root , 'QueryKey' )
			totalRecords = 0
			totalRecords = self.getEntry( root , 'Count' )
			try:
				if int(totalRecords) > 0:
			#		print totalRecords + "records found"
					self.setRetmax( int(totalRecords) , summaryBatchSize=summaryBatchSize )
					self.subset = entrezapi.esummary
					for self.retstart in range( 0 , int(totalRecords) , self.retmax ):
						self.action = self.buildWebEnvAction()
						url = self.buildURL()
				#		print url
						summaryResponse = self.submit()
						variants.update( self.getClinVarEntry() )
						#print "These are the ClinVar variants: " 
						#print variants
			except:
				print "CharGer Warning: No ClinVar records in this batch\n"
		return variants
	
	def parseClinVarTitle( self , DocumentSummary ):
#		print "biomine::webapi::entrez::entrezapi::parseClinVarTitle - " ,
		title = self.getEntry( DocumentSummary , 'title' )
#		print title
		codonPos = ""
		peptideRef = ""
		peptidePos = ""
		peptideAlt = ""
		var = mafvariant()
		residueMatches = re.search( "\((p\.\w+)\)" , title )
#		print "peptide variant: " ,
#		print residueMatches
		if residueMatches:
			hgvsp = residueMatches.groups()[-1]
#			print hgvsp
			[ peptideRef , peptidePos , peptideAlt ] = var.splitHGVSp( hgvsp )
		codonMatches = re.search( "(c\.\d+)" , title )
#		print "codon variant: " ,
#		print codonMatches
		if codonMatches:
			hgvsc = codonMatches.groups()[-1]
#			print hgvsc
			[ codonRef , codonPos , codonAlt ] = var.splitHGVSc( hgvsc )
		return { "title" : title , \
		"referencePeptide" : peptideRef , \
		"positionPeptide" : peptidePos , \
		"alternatePeptide" : peptideAlt , \
		"positionCodon" : codonPos }

	def getClinVarEntry( self ):
#		print "\twebapi::entrez::entrezapi::getClinVarEntry"
		root = self.getXMLRoot()
		variants = {}
		for DocumentSummary in root.iter( 'DocumentSummary' ):
			uid = DocumentSummary.attrib["uid"]
			var = clinvarvariant( uid=uid )
			self.getClinVarVariantEntry( var , DocumentSummary )
			self.getClinVarTraitEntry( var , DocumentSummary )
			self.getClinVarClinicalEntry( var , DocumentSummary )
			if var.genomicVar():
				variants[var.genomicVar()] = var
			else:
				variants["null"] = var
				print "entrez Warning: could not set ClinVar variant " + uid
		return variants
	def getClinVarVariantEntry( self , var , DocumentSummary ):
#		print "\twebapi::entrez::entrezapi::getClinVarVariantEntry"
		titleDetails = self.parseClinVarTitle( DocumentSummary )
#		print titleDetails
		var.referencePeptide = titleDetails["referencePeptide"]
		var.positionPeptide = titleDetails["positionPeptide"]
		var.alternatePeptide = titleDetails["alternatePeptide"]
		var.positionCodon = titleDetails["positionCodon"]
		for variation in DocumentSummary.iter( 'variation' ):
			for variation_xref in DocumentSummary.iter( 'variation_xref' ):
				dbs = self.getEntry( variation_xref , 'db_source' )
				if dbs == entrezapi.dbsnp:
					var.dbSNP = self.getEntry( variation_xref , 'db_id' )
					#print "dbSNP rs" + var.dbSNP
				#if dbs == entrezapi.omim:
					#var.omim = self.getEntry( variation_xref , 'db_id' )
					#print "OMIM " + var.omim
			for assembly_set in variation.iter( 'assembly_set' ):
				assembly_name = self.getEntry( assembly_set , 'assembly_name' )
				if assembly_name == self.assembly:
					var.chromosome = assembly_set.find( 'chr' ).text
					var.start = assembly_set.find( 'start' ).text
					var.stop = assembly_set.find( 'stop' ).text
					var.alternate = assembly_set.find( 'alt' ).text
					var.reference = assembly_set.find( 'ref' ).text
					#assembly_acc_ver = assembly_set.find( 'assembly_acc_ver' )
		for gene in DocumentSummary.iter( 'gene' ):
			var.gene = self.getEntry( gene , 'symbol' )
			var.strand = self.getEntry( gene , 'strand' )
		if not var.genomicVar():
			print "entrez Warning: no ClinVar variant entry"
	def getClinVarTraitEntry( self , var , DocumentSummary ):
#		print "\twebapi::entrez::entrezapi::getClinVarTraitEntry"
		try:
			for trait in DocumentSummary.iter( 'trait' ):
				trait_name = self.getEntry( trait , 'trait_name' )
				trait_xrefs = {}
				var.trait = { trait_name : trait_xrefs }
				for trait_xref in trait.iter( 'trait_xref' ):
					db_source = self.getEntry( trait_xref , 'db_source' )
					db_id = self.getEntry( trait_xref , 'db_id' )
					txr = {}
					if trait_name in trait_xrefs:
						txr = trait_xrefs[trait_name]
					txr.update( { db_source : db_id } )
					trait_xrefs.update( { trait_name : txr } )
					var.trait.update( { trait_name : trait_xrefs } )
		except:
			print "entrez Warning: no ClinVar trait entry"
			pass
	def getClinVarClinicalEntry( self , var , DocumentSummary ):
#		print "\twebapi::entrez::entrezapi::getClinVarClinicalEntry"
		try:
			for clinical_significance in DocumentSummary.iter( 'clinical_significance' ):
				var.clinical["description"] = self.getEntry( clinical_significance , 'description' ).strip()
				var.clinical["review_status"] = self.getEntry( clinical_significance , 'review_status' ).strip()
		except:
			print "entrez Warning: no ClinVar clinical entry"
			pass
	def getClinVarPubMedIDs( self , var , DocumentSummary ):
		try:
			return var.linkPubMed()
		except:
			print "entrez Error: no ClinVar uid"
			pass
			

	def searchPubMed( self , query ):
		self.subset = entrezapi.esearch
		self.database = entrezapi.pubmed
		self.action = self.buildSearchAction( query )
		#print self.subset
		#print self.database
		#print self.action
		return self.submit( )
	def searchClinVar( self , query , **kwargs ):
		self.subset = entrezapi.esearch
		self.database = entrezapi.clinvar
		self.action = self.buildSearchAction( query )
		#print self.subset
		#print self.database
		#print self.action
		return self.submit( )
	
	def getSNPGlobalMAF( self , rs ):
		self.addQuery( rs )
		self.subset = entrezapi.esummary
		self.submit()
		root = self.getXMLRoot()
		entries = {}
		for DocSum in root.iter( 'DocSum' ):
			print ""
