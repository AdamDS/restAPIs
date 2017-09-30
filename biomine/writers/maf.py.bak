#!/usr/bin/python
# ExAC Parser - 
# author: Adam D Scott (adamscott@wustl.edu)
# version: v0.0 - 2016*05*18
# AC=38,2;AC_AFR=1,0;AC_AMR=0,0;AC_Adj=5,0;AC_EAS=0,0;AC_FIN=0,0;AC_Het=5,0,0;AC_Hom=0,0;AC_NFE=3,0;AC_OTH=0,0;AC_SAS=1,0;AF=3.200e-04,1.684e-05;AN=118764;AN_AFR=2366;AN_AMR=2222;AN_Adj=28902;AN_EAS=1822;AN_FIN=600;AN_NFE=13932;AN_OTH=268;AN_SAS=7692;

from biomine.variant.clinvarvariant import clinvarvariant
from biomine.variant.vepvariant import vepvariant
from biomine.variant.vepconsequencevariant import vepconsequencevariant
from biomine.variant.mafvariant import mafvariant
from biomine.variant.variant import variant
import re

class maf(object):
	''' Example usage:
			vars = mafvariants()
			vars.getInputData( vcf=file )
	'''
	columnHeader = [ \
		"Hugo_Symbol" , \
		"Entrez_Gene_Id" , \
		"Center" , \
		"NCBI_Build" , \
		"Chromosome" , \
		"Start_Position" , \
		"End_Position" , \
		"Strand" , \
		"Variant_Classification" , \
		"Variant_Type" , \
		"Reference_Allele" , \
		"Tumor_Seq_Allele1" , \
		"Tumor_Seq_Allele2" , \
		"dbSNP_RS" , \
		"dbSNP_Val_Status" , \
		"Tumor_Sample_Barcode" , \
		"Matched_Norm_Sample_Barcode" , \
		"Match_Norm_Seq_Allele1" , \
		"Match_Norm_Seq_Allele2" , \
		"Tumor_Validation_Allele1" , \
		"Tumor_Validation_Allele2" , \
		"Match_Norm_Validation_Allele1" , \
		"Match_Norm_Validation_Allele2" , \
		"Verification_Status" , \
		"Validation_Status" , \
		"Mutation_Status" , \
		"Sequencing_Phase" , \
		"Sequence_Source" , \
		"Validation_Method" , \
		"Score" , \
		"BAM_File" , \
		"Sequencer" ]
	def __init__( self , **kwargs ):
		self.customHeader = kwargs.get( 'customHeader' , [] )
		self.filename = kwargs.get( 'filename' , "" )

	def write( self , **kwargs ):
		fh = open( outputFile , 'w' )
		customHeaders = kwargs.get( 'customHeader' , [] )
		fh.write( delim.join( mafvariant.columnHeader ) )
		if customHeaders:
			fh.write( "\t" + delim.join( customHeaders ) )
		fh.write( "\n" )
		na = "NA"
		for var in mafvariants:
			fields = [ \
				getMember( var.gene ) , \
				na , \
				na , \
				getMember( var.chromosome ) , \
				getMember( var.start ) , \
				getMember( var.stop ) , \
				getMember( var.strand ) , \
				getClassification( var ) , \
				getMember( var.reference ) , \
				getMember( var.reference ) , \
				getMember( var.alternate ) , \
				getMember( var.dbsnp ) , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
				na , \
			]
			line = "\t".join( fields )
			fh.write( line + "\n" )
		fh.close()
		print "finished generating .maf\n"

	@staticmethod
	def getMember( self , member , **kwargs ):
		na = kwargs.get( 'na' , "NA" )
		if member:
			return member
		return na	
	def getClassification( self , var , **kwargs ):
		classes = [ "Frame_Shift_Del" , \
					"Frame_Shift_Ins" , \
					"In_Frame_Del" , \
					"In_Frame_Ins" , \
					"Missense_Mutation" , \
					"Nonsense_Mutation" , \
					"Silent" , \
					"Splice_Site" , \
					"Translation_Start_Site" , \
					"Nonstop_Mutation" , \
					"RNA" , \
					"Targeted_Region" ]
		if type( var ) is mafvariant:
			modAlt = len( var.alternate ) % 3
			modRef = len( var.reference ) % 3
			if var.reference == "-":
				if var.alternate == "-": #->-; shouldn't happen
					print "WARNING: no reference or alternate. Classification set to Silent."
					return classes[6]
				else:
					if modAlt == 0: #->xyz; In_Frame_Ins
						return classes[3]
					else: #->x or ->xy; Frame_Shift_Ins
						return classes[1]
			else:
				if var.alternate == "-":
					if modRef == 0: #xyz>-; In_Frame_Del
						return classes[2]
					else: #x>- or xy>-; Frame_Shift_Del
						return classes[0]
				else:
					if modRef == 0:
						if modAlt == 0: #xyz>xyz; in frame complex indel - In_Frame_Del
							return classes[2]
						else: #xyz>x or xyz>xy; frame shift complex indel - Frame_Shift_Del
							return classes[0]
					elif len( var.referencePeptide ) == 1 and len( var.alternatePeptide ) == 1: #x>y; SNP
						if var.referencePeptide == var.alternatePeptide: #Silent
							return classes[6]
						else: #Missense_Mutation
							return classes[4]
		else:
			if var.reference == "-":
				if var.alternate == "-": #->-; shouldn't happen
					print "WARNING: no reference or alternate. Classification set to Silent."
					return classes[6]
				else:
					if mod( len( var.alternate ) , 3 ) == 0: #->xyz; In_Frame_Ins
						return classes[3]
					else: #->x or ->xy; Frame_Shift_Ins
						return classes[1]
			else:
				if var.alternate == "-":
					if mod( len( var.reference ) , 3 ) == 0: #xyz>-; In_Frame_Del
						return classes[2]
					else: #x>- or xy>-; Frame_Shift_Del
						return classes[0]
				else:
					if mod( len( var.reference ) , 3 ) == 0:
						if mod( len( var.alternate ) , 3 ) == 0: #xyz>xyz; in frame complex indel
							return classes[2]
						else:
							return classes[0]
						return classes[2]
					elif len( var.referencePeptide ) == 1 and len( var.alternatePeptide ) == 1:
						if var.referencePeptide == var.alternatePeptide:
							return classes[6]
						else:
							return classes[4]
