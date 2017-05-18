#!/usr/bin/python
# ExAC Variant - 
# author: Adam D Scott (adamscott@wustl.edu)
# version: v0.0 - 2016*05*20
# AC=38,2;AC_AFR=1,0;AC_AMR=0,0;AC_Adj=5,0;AC_EAS=0,0;AC_FIN=0,0;AC_Het=5,0,0;AC_Hom=0,0;AC_NFE=3,0;AC_OTH=0,0;AC_SAS=1,0;AF=3.200e-04,1.684e-05;AN=118764;AN_AFR=2366;AN_AMR=2222;AN_Adj=28902;AN_EAS=1822;AN_FIN=600;AN_NFE=13932;AN_OTH=268;AN_SAS=7692;

from biomine.variant.clinvarvariant import clinvarvariant
from biomine.variant.vepvariant import vepvariant
from biomine.variant.mafvariant import mafvariant
from biomine.variant.variant import variant
import re

class exacvariant( variant ):
	def __init__( self , **kwargs ):
		super( exacvariant , self ).__init__( **kwargs )
		self.counts = populationmeasures( valueType = "AC" )
		self.numbers = populationmeasures( valueType = "AN" )
		self.frequency = 0
		self.vepvariant = None
		self.clinvarvariant = None
		self.mafvariant = None
		self.POS = None
		self.REF = None
		self.ALT = None
		self.allele = None
		self.QUAL = None
		self.FILTER = None
		self.INFO = None
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( exacvariant , self ).copyInfo( aParentVariant )
	def copyInfo( self , copy ):
		super( exacvariant , self ).copyInfo( copy )
		self.counts = copy.counts
		self.numbers = copy.numbers
		self.frequency = copy.frequency
		self.vepvariant = copy.vepvariant
		self.clinvarvariant = copy.clinvarvariant
		self.mafvariant = copy.mafvariant
		self.POS = self.POS
		self.REF = self.REF
		self.ALT = self.ALT
		self.allele = self.allele
		self.QUAL = self.QUAL
		self.FILTER = self.FILTER
		self.INFO = self.INFO
	#def print( self , **kwargs ):
	#	delim = kwargs.get( 'delim' , "\t" )
	#	minimal = kwargs.get( 'minimal' , False )
	#	if not minimal:
	#		super( variant , self ).print( **kwargs ) ,
	#	print delim ,
	#	print self.counts ,
	#	print delim ,
	#	print self.numbers ,
	#	print delim ,
	#	print self.frequency

	def setCounts( self , info , **kwargs ):
		self.counts.setPopulation( info , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.afr , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.amr , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.eas , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.fin , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.nfe , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.oth , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.sas , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.adj , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.het , **kwargs )
		self.counts.setPopulation( info , pop = populationmeasures.hom , **kwargs )
	def setNumbers( self , info , **kwargs ):
		self.numbers.setPopulation( info , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.afr , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.amr , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.eas , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.fin , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.nfe , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.oth , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.sas , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.adj , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.het , **kwargs )
		self.numbers.setPopulation( info , pop = populationmeasures.hom , **kwargs )
	def getPopulationCount( self , **kwargs ):
		pop = kwargs.get( 'pop' , '' )
		return self.counts.getPopulation( pop = pop )
	def getPopulationNumber( self , **kwargs ):
		pop = kwargs.get( 'pop' , '' )
		return self.numbers.getPopulation( pop = pop )
	def setFrequency( self , info , **kwargs ):
		whichAlternate = kwargs.get( 'alti' , 0 )
		frequency = info.get( "AF" , 0 )
		if type( frequency ) is list:
			self.frequency = float( frequency[whichAlternate] )
		else:
			self.frequency = float( frequency )
	def getFrequency( self , info , **kwargs ):
		return self.frequency

class populationmeasures( object ):
	afr = "AFR"
	amr = "AMR"
	eas = "EAS"
	fin = "FIN"
	nfe = "NFE"
	oth = "OTH"
	sas = "SAS"
	adj = "Adj"
	het = "Het"
	hom = "Hom"
	count = "AC"
	number = "AN"
	def __init__( self , **kwargs ):
		self.total = 0
		self.afr = 0
		self.amr = 0
		self.eas = 0
		self.fin = 0
		self.nfe = 0
		self.oth = 0
		self.sas = 0
		self.adj = 0
		self.het = 0
		self.hom = 0
		self.valueType = kwargs.get( "valueType" , "AC" )
	def __str__( self , **kwargs ):
		delim = kwargs.get( 'delim' , "\t" )
		string = str( self.total ) + delim ,
		str( self.afr ) + delim ,
		str( self.amr ) + delim ,
		str( self.eas ) + delim ,
		str( self.fin ) + delim ,
		str( self.nfe ) + delim ,
		str( self.oth ) + delim ,
		str( self.sas ) + delim ,
		str( self.adj )
		if self.valueType == "AC":
			string += delim + str( self.het ) + delim ,
			str( self.hom )
		return string
	def __all__( self , **kwargs ):
		annotated = kwargs.get( 'annotated' , False )
		asStr = kwargs.get( 'str' , False )
		if annotated:
			return [ self.annotated( ) , \
					 self.annotated( pop=populationmeasures.adj ) , \
					 self.annotated( pop=populationmeasures.afr ) , \
					 self.annotated( pop=populationmeasures.amr ) , \
					 self.annotated( pop=populationmeasures.eas ) , \
					 self.annotated( pop=populationmeasures.fin ) , \
					 self.annotated( pop=populationmeasures.nfe ) , \
					 self.annotated( pop=populationmeasures.oth ) , \
					 self.annotated( pop=populationmeasures.sas ) , \
					 self.annotated( pop=populationmeasures.het ) , \
					 self.annotated( pop=populationmeasures.hom ) ]
		else:
			if asStr:
				return [ str( self.total ) , \
						 str( self.afr ) , \
						 str( self.amr ) , \
						 str( self.eas ) , \
						 str( self.fin ) , \
						 str( self.nfe ) , \
						 str( self.oth ) , \
						 str( self.sas ) , \
						 str( self.adj ) , \
						 str( self.het ) , \
						 str( self.hom ) ]
			else:
				return [ self.total , \
						 self.afr , \
						 self.amr , \
						 self.eas , \
						 self.fin , \
						 self.nfe , \
						 self.oth , \
						 self.sas , \
						 self.adj , \
						 self.het , \
						 self.hom ]
	def annotated( self , **kwargs ):
		pop = kwargs.get( 'pop' , '' )
		delim = kwargs.get( 'delim' , "=" )
		annoName = self.fieldName( pop ) + delim
		if pop == populationmeasures.afr:
			return annoName + str( self.afr )
		if pop == populationmeasures.amr:
			return annoName + str( self.amr )
		if pop == populationmeasures.eas:
			return annoName + str( self.eas )
		if pop == populationmeasures.fin:
			return annoName + str( self.fin )
		if pop == populationmeasures.nfe:
			return annoName + str( self.nfe )
		if pop == populationmeasures.oth:
			return annoName + str( self.oth )
		if pop == populationmeasures.sas:
			return annoName + str( self.sas )
		if pop == populationmeasures.adj:
			return annoName + str( self.adj )
		if self.valueType == populationmeasures.count:
			if pop == populationmeasures.het:
				return annoName + str( self.het )
			if pop == populationmeasures.hom:
				return annoName + str( self.hom )
		return annoName + str( self.total )
	def __str__( self , **kwargs ):
		delim = kwargs.get( 'delim' , '\t' )
		print delim.join( map( str , self.__all__() ) )
	def fieldName( self , pop ):
		if pop:
			return self.valueType + "_" + pop
		else:
			return self.valueType
	def setPopulation( self , info , **kwargs ):
		pop = kwargs.get( 'pop' , '' )
		name = self.fieldName( pop )
		whichAlternate = kwargs.get( 'alti' , 0 )
		valuesStr = info.get( name , 0 )
		value = 0
		if type( valuesStr ) is list:
			value = valuesStr[whichAlternate]
		else:
			value = valuesStr
		if pop == populationmeasures.afr:
			self.afr = value
		if pop == populationmeasures.amr:
			self.amr = value
		if pop == populationmeasures.eas:
			self.eas = value
		if pop == populationmeasures.fin:
			self.fin = value
		if pop == populationmeasures.nfe:
			self.nfe = value
		if pop == populationmeasures.oth:
			self.oth = value
		if pop == populationmeasures.sas:
			self.sas = value
		if pop == populationmeasures.adj:
			self.adj = value
		if self.valueType == populationmeasures.count:
			if pop == populationmeasures.het:
				self.het = value
			if pop == populationmeasures.hom:
				self.hom = value
		if not pop:
			self.total = value
	def getPopulation( self , **kwargs ):
		pop = kwargs.get( 'pop' , '' )
		name = self.fieldName( pop )
		if pop == populationmeasures.afr:
			return self.afr
		if pop == populationmeasures.amr:
			return self.amr
		if pop == populationmeasures.eas:
			return self.eas
		if pop == populationmeasures.fin:
			return self.fin
		if pop == populationmeasures.nfe:
			return self.nfe
		if pop == populationmeasures.oth:
			return self.oth
		if pop == populationmeasures.sas:
			return self.sas
		if pop == populationmeasures.adj:
			return self.getAdjusted( **kwargs )
		if self.valueType == populationmeasures.count:
			if pop == populationmeasures.het:
				return self.getHeterozygous( **kwargs )
			if pop == populationmeasures.hom:
				return self.getHomozygous( **kwargs )
		return self.getTotal( **kwargs )
	def getAdjusted( self , **kwargs ):
		return self.adj
	def getHeterozygous( self , **kwargs ):
		return self.het
	def getHomozygous( self , **kwargs ):
		return self.hom
	def getTotal( self , **kwargs ):
		return self.total
