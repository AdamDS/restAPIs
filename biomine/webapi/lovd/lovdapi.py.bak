#!/usr/bin/python
# author: Adam D Scott
# first created: 2016*1*11
#####
#	Operation					Function				Works With
#####
### Inherited
#	endpoint	"exac.broadinstitute.org"
#	subset		"/variant/"
#	action		

from biomine.webapi.webapi import webapi
import xml.etree.ElementTree as ET
import json
import AdvancedHTMLParser
import biomine.variant

class exacapi(webapi):
	endpoint = "http://exac.broadinstitute.org"
	variant = "/variant/"
	def __init__(self,**kwargs):
		subset = kwargs.get("subset",'')
		if not subset:
			super(exacapi,self).__init__(exacapi.endpoint,exacapi.variant)
		else:
			if ( subset == exacapi.variant ):
				super(exacapi,self).__init__(exacapi.endpoint,subset)
			else:
				print "biomine ERROR: bad subset. webapi.subset initializing to variant association results"
				super(exacapi,self).__init__(exacapi.endpoint,exacapi.variant)

	def setSubset(self,subset):
		self.subset = subset
		self.action = ""
	def resetURL(self):
		self.action = ""
	def buildQuery( self , var ):
		self.action = '-'.join( [ var.chromosome , var.start , var.reference , var.mutant ] )

	def getPage( self , var ):
		self.buildQuery( var )
		self.submit()
		return self.parseHTMLResponse()
	
	def getAlleleFrequencies( self , variants ):
		entries = {}
		for var in variants:
			genVar = var.genomicVar()
			entries[genVar] = self.getAlleleFrequency( var )
		return entries
		
	def getAlleleFrequency( self , var ):
		page = self.getPage( var )
		if page:
			dl = page.getElementsByTagName( 'dl' )
			if dl:
				dt = dl.getElementsByTagName( 'dt' )
				if dt:
					j = 0
					for i in range( 0 , len( dt ) ):
						if dt[i].innerHTML == "Allele Frequency":
							j = i
					dd = dl.getElementsByTagName( 'dd' )
					if dd:
						alleleFrequency = dd[j].innerHTML
						return float(alleleFrequency)
		return
		print "biomine Error: could not get allele frequency from exac"
