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
from biomine.variant.exacvariant import exacvariant
import pysam
import vcf
import re
import os.path
import gc

class exacparser(object):
	''' Example usage:
			vars = exacvariants()
			vars.getInputData( vcf=file )
	'''
	severeRank = [ 	"transcript_ablation" , \
					"splice_acceptor_variant" , \
					"splice_donor_variant" , \
					"stop_gained" , \
					"frameshift_variant" , \
					"stop_lost" , \
					"start_lost" , \
					"transcript_amplification" , \
					"inframe_insertion" , \
					"inframe_deletion" , \
					"missense_variant" , \
					"protein_altering_variant" , \
					"splice_region_variant" , \
					"incomplete_terminal_codon_variant" , \
					"stop_retained_variant" , \
					"synonymous_variant" , \
					"coding_sequence_variant" , \
					"mature_miRNA_variant" , \
					"5_prime_UTR_variant" , \
					"3_prime_UTR_variant" , \
					"non_coding_transcript_exon_variant" , \
					"intron_variant" , \
					"NMD_transcript_variant" , \
					"non_coding_transcript_variant" , \
					"upstream_gene_variant" , \
					"downstream_gene_variant" , \
					"TFBS_ablation" , \
					"TFBS_amplification" , \
					"TF_binding_site_variant" , \
					"regulatory_region_ablation" , \
					"regulatory_region_amplification" , \
					"feature_elongation" , \
					"regulatory_region_variant" , \
					"feature_truncation" , \
					"intergenic_variant" ]
	def __init__( self , **kwargs ):
		self.variants = kwargs.get( 'variants' , [] )
		#self.clinvarVariants = kwargs.get( 'clinvarVariants' , {} )
		self.vcfHeaderInfo = kwargs.get( 'vcfHeaderInfo' , [] )
		self.vcfKeyIndex = kwargs.get( 'vcfKeyIndex' , {} )

### Retrieve input data from user ###
	def getInputData( self , **kwargs ):
		vcfFile = kwargs.get( 'vcf' , "" )
		vepDone = False
		exacDone = False
		if vcfFile:
			vepDone = self.readVCF( vcfFile , **kwargs )
			exacDone=False # currently only has 1000G
		return ( vepDone , exacDone )

	def readVCF( self , inputFile , **kwargs ):
		( inFile , outFH ) = self.initializeFile( inputFile , **kwargs )
		self.getVCFMetaData( inFile , **kwargs )
		kwargs['outFH'] = outFH
		self.getVCFVariants( inFile , **kwargs )
	def initializeFile( self , inputFile , **kwargs ):
		writeToFile = kwargs.get( 'writeToFile' , "biomine.parsed.exac.default.tsv" )
		outFH = open( writeToFile , 'a' )
		if ( outFH.closed ):
			print "file handle closed: " + writeToFile
		inFile = None;
		if ( re.match( "\.gz" , inputFile ) ):
			inFile = vcf.Reader( open( inputFile , 'r' ) , compressed=True )
		else:
			inFile = vcf.Reader( open( inputFile , 'r' ) )
		return ( inFile , outFH )
	def getVCFMetaData( self , inFile , **kwargs ):
		vepDone = False
		metadata = inFile.metadata
		for pairs in metadata:
			if pairs == 'VEP': #TODO: make this a function, passed as an argument to getVCFMetaData
				print "This .vcf has VEP annotations!"
				infos = inFile.infos
				for info_ID in infos.items():
					if info_ID[0] == "CSQ": #CSQ tag the VEP annotation, probably means consequence
						csq = info_ID
						Info = csq[1] #Info(...)
						if Info:
							desc = Info[3] #Consequence type...Format: Allele|Gene|...
							keysString = desc.split( "Format: " )[1]
							self.vcfHeaderInfo = keysString.split( "|" )
							self.vcfKeyIndex = {}
							i = 0
							for key in self.vcfHeaderInfo:
								self.vcfKeyIndex[key] = i
								i = i + 1
	def getVCFVariants( self , inFile , **kwargs ): 
		for record in inFile:
			self.collectAlternates( record , **kwargs )
		self.finalWriteVariants( **kwargs )
	def collectAlternates( self , record , **kwargs ):
		chrom = record.CHROM
		reference = record.REF
		alternates = record.ALT
		alti = -1
		for alternate in alternates:
			alti += 1
			begin = record.start
			end = record.end
			( begin , end , ref , alt ) = self.getPositionsAndAlleles( record , alti , **kwargs )
			self.getVariant( record , begin , end , ref , alt , alti , **kwargs )
	def getVariant( self , record , begin , end , ref , alt , alti , **kwargs ):
		var = exacvariant( \
			chromosome = record.CHROM , \
			start = begin , \
			stop = end , \
			dbsnp = record.ID , \
			reference = ref , \
			alternate = alt , \
		)
		self.getCSQ( var , record , begin , end , ref , alt , alti , **kwargs )
		var.POS = record.POS
		var.REF = record.REF
		var.ALT = record.ALT[alti]
		var.allele = alti
		var.QUAL = record.QUAL
		var.FILTER = record.FILTER
		var.INFO = record.INFO
		self.setWriteVariant( var , **kwargs )
	def getPositionsAndAlleles( self , record , alti , **kwargs ):
		begin = record.start + 1 #1-base beginning of ref
		end = record.end #0-base ending of ref
		ref = record.REF
		alt = record.ALT[alti]
		if record.is_indel: #indel
			positionOfMismatch = self.forwardMismatchDetection( ref , alt , **kwargs )
			revPositionOfMismatch = self.reverseMismatchDetection( ref , alt , **kwargs )
			( begin , end , ref , alt ) = self.determineOverlap( positionOfMismatch , revPositionOfMismatch , ref , alt , begin , alti , **kwargs )
		elif record.is_deletion:
			begin = start + 1
			end = stop
			ref = ref[1:] #assumes only one base overlap
			alt = "-"
		return ( begin , end , ref , alt )
	def forwardMismatchDetection( self , reference , alternate , **kwargs ):
		refBases = str( reference )
		altBases = str( alternate )
		lengthOfReference = len( refBases )
		lengthOfAlternate = len( altBases )
		positionOfMismatch = 0
		if ( lengthOfReference > 1 or lengthOfAlternate > 1 ):
			referenceIndex = 0
			alternateIndex = 0
			referenceBase = refBases[ referenceIndex ]
			alternateBase = altBases[ alternateIndex ]
			mismatch = False
			while ( not mismatch and ( referenceIndex < lengthOfReference and alternateIndex < lengthOfAlternate ) ):
				if ( referenceBase == alternateBase ): #match & haven't mismatched
					positionOfMismatch += 1
				else: #mismatch
					mismatch = True
				#print '  '.join( ( "forward: " , str( positionOfMismatch ) , str( referenceIndex ) , str( referenceBase ) , \
				#	str( alternateIndex ) , str( alternateBase ) , str( mismatch ) ) )
				referenceIndex += 1
				alternateIndex += 1
				if ( referenceIndex < lengthOfReference ):
					referenceBase = refBases[ referenceIndex ]
				else:
					referenceBase = "-"
				if ( alternateIndex < lengthOfAlternate ):
					alternateBase = altBases[ alternateIndex ]
				else:
					alternateBase = "-"
		return positionOfMismatch
	def reverseMismatchDetection( self , reference , alternate , **kwargs ):
		refBases = str( reference )
		altBases = str( alternate )
		lengthOfReference = len( refBases )
		lengthOfAlternate = len( altBases )
		revPositionOfMismatch = lengthOfAlternate
		referenceIndex = lengthOfReference - 1
		alternateIndex = lengthOfAlternate - 1
		referenceBase = refBases[ referenceIndex ]
		alternateBase = altBases[ alternateIndex ]
		mismatch = False
		while ( not mismatch and ( referenceIndex >= 0 and alternateIndex >= 0 ) ):
			if ( referenceBase == alternateBase ): #match & haven't mismatched
				revPositionOfMismatch -= 1
			else: #mismatch
				mismatch = True
			#print '  '.join( ( "reverse: " , str( revPositionOfMismatch ) , str( referenceIndex ) , str( referenceBase ) , \
			#	str( alternateIndex ) , str( alternateBase ) , str( mismatch ) ) )
			referenceIndex -= 1
			alternateIndex -= 1
			if ( referenceIndex > 0 ):
				referenceBase = refBases[ referenceIndex ]
			else:
				referenceBase = "-"
			if ( alternateIndex > 0 ):
				alternateBase = altBases[ alternateIndex ]
			else:
				alternateBase = "-"
		return revPositionOfMismatch
	def errprint( self , case , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end ):
		print "case" + str( case ) + ": " + '\t'.join( ( str( begin ) , ref , str( positionOfMismatch ) , str( end ) , alt , str( revPositionOfMismatch ) ) )
	def errfollowup( self , case , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end ):
		print "followup" + str( case ) + ": " + '\t'.join( ( str( begin ) , ref , str( positionOfMismatch ) , str( end ) , alt , str( revPositionOfMismatch ) ) )
	def determineOverlap( self , positionOfMismatch , revPositionOfMismatch , reference , alternate , start , alti , **kwargs ):
		ref = reference
		alt = str( alternate )
		begin = start
		end = start
		if alt == "None":
			alt = None
		( begin , end , ref , alt ) = self.overlapLogic( positionOfMismatch , revPositionOfMismatch , ref , alt , begin , end )
		return ( begin , end , ref , alt )

	def mismatchDetection( self , reference , alternate , **kwargs ):
		reverse = kwargs.get( 'reverse' , False )
		refBases = str( reference )
		altBases = str( alternate )
		if ( reverse ):
			refBases = str( reference[::-1] )
			altBases = str( alternate[::-1] )
		lengthOfReference = len( refBases )
		lengthOfAlternate = len( altBases )
		positionOfMismatch = 0
		if ( lengthOfReference > 1 or lengthOfAlternate > 1 ):
			referenceIndex = 0
			alternateIndex = 0
			referenceBase = refBases[ referenceIndex ]
			alternateBase = altBases[ alternateIndex ]
			mismatch = False
			while ( not mismatch and ( referenceIndex < lengthOfReference and alternateIndex < lengthOfAlternate ) ):
				if ( referenceBase == alternateBase ): #match & haven't mismatched
					positionOfMismatch += 1
				else: #mismatch
					mismatch = True
				#print '  '.join( ( "forward: " , str( positionOfMismatch ) , str( referenceIndex ) , str( referenceBase ) , \
				#	str( alternateIndex ) , str( alternateBase ) , str( mismatch ) ) )
				referenceIndex += 1
				alternateIndex += 1
				if ( referenceIndex < lengthOfReference ):
					referenceBase = refBases[ referenceIndex ]
				else:
					referenceBase = "-"
				if ( alternateIndex < lengthOfAlternate ):
					alternateBase = altBases[ alternateIndex ]
				else:
					alternateBase = "-"
		return ( positionOfMismatch , refBases[ positionOfMismatch : ] , altBases[ positionOfMismatch : ] )

	def overlapLogic2( self , positionOfMismatch , revPositionOfMismatch , ref , alt , begin , end ):
		lengthOfReference = len( ref )
		lengthOfAlternate = len( alt )
		if ( lengthOfReference < lengthOfAlternate ): #C>CGAGA|1,5 , CCCCT>CCCCTCCCT , CCCCCACCCCA>CCCCCACCCCACCCCA , CCC>CCCCC , ATATAT>ATGCATATAT|-:GCAT
			ref = "-"
			begin += positionOfMismatch + 1
			end = begin + 1
			if ( positionOfMismatch < revPositionOfMismatch ):
				alt = alt[ positionOfMismatch : lengthOfAlternate ] #1:5 , 5:9 , 11:15
			elif ( positionOfMismatch > revPositionOfMismatch ):
				alt = alt[ revPositionOfMismatch : revPositionOfMismatch + 1 ]
			else:
				alt = alt[ revPositionOfMismatch : revPositionOfMismatch + 1 ]
		elif ( lengthOfReference > lengthOfAlternate ): #CCCCCACCCCA>C , TAA>TA , GTCCTCCTCGCCC>GTCCTCGCCC , AAAAA>AA|2,4
			ref = ref[ positionOfMismatch : revPositionOfMismatch ]
			begin += positionOfMismatch + 1
			alt = "-"
			end = begin + lengthOfReference
		else: #TAA>TAG , TAA>TGA , TAA>GAA , TAA>GCA , TAA>TGC , TAA>GAG
			ref = ref[ positionOfMismatch : revPositionOfMismatch ]
			begin += positionOfMismatch + 1
			alt = alt[ positionOfMismatch : revPositionOfMismatch ]
			end = begin + 1
		return ( begin , end , ref , alt )

	def overlapLogic( self , positionOfMismatch , revPositionOfMismatch , ref , alt , begin , end ):
		lengthOfReference = len( ref )
		lengthOfAlternate = len( alt )
		if ( ( positionOfMismatch - 1 ) < ( revPositionOfMismatch + 1 ) ): #C>CGAGA , CG>CGAGA
			if ( positionOfMismatch == 1 and revPositionOfMismatch == 1 ): #CCCCCACCCCA>C
				self.errprint( 1 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
				ref = ref[ positionOfMismatch : ]
				alt = "-"
				begin += positionOfMismatch
				self.errfollowup( 1 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
			elif ( positionOfMismatch == 0 and revPositionOfMismatch == 0 ):
				self.errprint( 2 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
				ref = ref[ positionOfMismatch ]
				alt = alt[ revPositionOfMismatch ]
				self.errfollowup( 2 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
			elif ( positionOfMismatch == revPositionOfMismatch ): #CCCCT>CCCCTCCCT
				self.errprint( 3 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
				ref = "-"
				alt = alt[ revPositionOfMismatch : ]
				begin += positionOfMismatch - 1
				end = begin + 1
				self.errfollowup( 3 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
			else:
				if ( ( positionOfMismatch < revPositionOfMismatch ) ): #G>GCACACA , CCCCT>CCCCTCCCTCCC
					self.errprint( 4 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
					ref = "-"
					begin += positionOfMismatch
					alt = alt[ positionOfMismatch : revPositionOfMismatch ]
					end = begin + 1
					self.errfollowup( 4 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
				else: #TAA>TA , TCCC>T
					self.errprint( 5 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
					ref = ref[ revPositionOfMismatch : positionOfMismatch ]
					begin += positionOfMismatch + 1
					alt = "-"
					end += lengthOfReference
					self.errfollowup( 5 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
		elif ( positionOfMismatch > revPositionOfMismatch ): #TT>CGAGA, p==0, revPositionOfMismatch==1
			if ( len( ref ) >= len( alt ) ): #CCCTCCTCCTCGT>CCCTCCTCGT , GTCCTCCTCGCCC>GTCCTCGCCC
				self.errprint( 6 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
				ref = ref[ positionOfMismatch : lengthOfReference - lengthOfAlternate + positionOfMismatch ]
				alt = "-"
				begin += positionOfMismatch + 1
				end = begin - len( ref ) + 1
				self.errfollowup( 6 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
			else: #CCCCCACCCCA>CCCCCACCCCACCCCA
				self.errprint( 7 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
				ref = "-"
				alt = alt[ positionOfMismatch : ]
				begin += positionOfMismatch -1
				end = begin + 1
				self.errfollowup( 7 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
		else:
			self.errprint( 8 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
			ref = "-"
			alt = alt[ positionOfMismatch : ] 
			end += 1
			self.errfollowup( 8 , ref , alt , positionOfMismatch , revPositionOfMismatch , begin , end )
		return ( begin , end , ref , alt )

	def getCSQ( self , var , record , begin , end , ref , alt , alti , **kwargs ):
		info = record.INFO
		csq = info.get( 'CSQ' , "noCSQ" )
		preVEP = []
		if not csq == "noCSQ":
			vepDone = True
			exacDone = True
			var.vepVariant = vepvariant()
			for thisCSQ in csq:
				values = thisCSQ.split( "|" )
				var.vcfInfo = values
				aas = [None , None] 
				if self.getVCFKeyIndex( values , "Amino_acids" ): #8 => Amino_acids
					aas = self.getVCFKeyIndex( values , "Amino_acids" ).split("/") 
					if len( aas ) > 1:
						aas[0] = mafvariant().convertAA( aas[0] )
						aas[1] = mafvariant().convertAA( aas[1] )
					else:
						#28 => HGVSc
						#29 => HGVSp
						hgvsp = self.getVCFKeyIndex( values , "HGVSp" ).split( ":" )
						changep = None
						if len( hgvsp ) > 1:
							changep = re.match( "p\." , hgvsp[1] )
						if changep:
							aas = mafvariant().splitHGVSp( hgvsp[1] )
							aas[0] = mafvariant().convertAA( aas[0] )
							aas[2] = mafvariant().convertAA( aas[2] )
						else:
							aas.append( None )
							needVEP = True
							preVEP.append( var )
				exons = [None , None]
				if self.getVCFKeyIndex( values , "EXON" ): #25 => EXON
					exons = self.getVCFKeyIndex( values , "EXON" ).split( "/" )
					if len( exons ) == 1:
						exons.append(None)
				introns = [None , None]
				if self.getVCFKeyIndex( values , "INTRON" ): #26 => INTRON
					introns = self.getVCFKeyIndex( values , "INTRON" ).split( "/" )
					if len( introns ) == 1:
						introns.append(None)
				siftStuff = [None , None]
				if self.getVCFKeyIndex( values , "SIFT" ):
					siftStuff = self.getVCFKeyIndex( values , "SIFT" ).split( "(" ) 
					if len( siftStuff ) == 1:
						siftStuff.append( None )
					else:
						siftStuff[1] = siftStuff[1].rstrip( ")" )
				polyPhenStuff = [None , None]
				if self.getVCFKeyIndex( values , "PolyPhen" ):
					polyPhenStuff = self.getVCFKeyIndex( values , "PolyPhen" ).split( "(" ) 
					if len( polyPhenStuff ) == 1:
						polyPhenStuff.append( None )
					else:
						polyPhenStuff[1] = polyPhenStuff[1].rstrip( ")" )

				vcv = vepconsequencevariant( \
					chromosome = record.CHROM , \
					start = begin , \
					stop = end , \
					dbsnp = record.ID , \
					reference = ref , \
					alternate = alt , \
					gene_id=self.getVCFKeyIndex( values , "Gene" ) , \
					transcriptCodon=self.getVCFKeyIndex( values , "Feature" ) , \
					consequence_terms=self.getVCFKeyIndex( values , "Consequence" ).split( "&" ) , \
					positionCodon=self.getVCFKeyIndex( values , "cDNA_position" ) , \
					positionPeptide=self.getVCFKeyIndex( values , "Protein_position" ) , \
					referencePeptide=aas[0] , \
					alternatePeptide=aas[1] , \
					strand=self.getVCFKeyIndex( values , "STRAND" ) , \
					gene=self.getVCFKeyIndex( values , "SYMBOL" ) , \
					gene_symbol_source=self.getVCFKeyIndex( values , "SYMBOL_SOURCE" ) , \
					hgnc_id=self.getVCFKeyIndex( values , "HGNC_ID" ) , \
					biotype=self.getVCFKeyIndex( values , "BIOTYPE" ) , \
					canonical=self.getVCFKeyIndex( values , "CANONICAL" ) , \
					ccds=self.getVCFKeyIndex( values , "CCDS" ) , \
					transcriptPeptide=self.getVCFKeyIndex( values , "ENSP" ) , \
					predictionSIFT=siftStuff[0] , \
					scoreSIFT=siftStuff[1] , \
					predictionPolyphen=polyPhenStuff[0] , \
					scorePolyphen=polyPhenStuff[1] , \
					exon=exons[0] , \
					totalExons=exons[1] , \
					intron=introns[0] , \
					totalIntrons=introns[1] , \
				)

				var.alleleFrequency = self.getVCFKeyIndex( values , "GMAF" )
				var.vepVariant.consequences.append( vcv )
			self.determineMostSevere( var , **kwargs )
			self.setAlleleMeasures( var , info , **kwargs )
		return None

	def determineMostSevere( self , var , **kwargs ):
		mostSevere = None
		rankMostSevere = 10000
		mostSevereCons = exacparser.severeRank[-1]
		for cons in var.vepVariant.consequences:
			for term in cons.terms:
				if term in exacparser.severeRank:
					rank = exacparser.severeRank.index( term )
				else:
					rank = 10000
				if rank < rankMostSevere:
					mostSevere = cons
					rankMostSevere = rank
					mostSevereCons = term
				elif rank == rankMostSevere:
					if cons.canonical:
						mostSevere = cons
		self.setMostSevere( var , mostSevere , mostSevereCons , **kwargs )
		return None
			
	def setWriteVariant( self , var , **kwargs ):
		writeToFile = kwargs.get( 'writeToFile' , "" )
		writeFunction = kwargs.get( 'writeFunction' , None )
		bufferSize = kwargs.get( 'bufferSize' , 1000 )
		if ( writeToFile and writeFunction ):
			if ( len( self.variants ) > bufferSize ):
				print gc.get_count()
				print "have " + str( len( self.variants ) ) + "; now write"
				for v in self.variants:
					writeFunction( v , **kwargs )
				print gc.get_count()
				print gc.garbage
				self.variants = []
				print "flushed variants (" + str( len( self.variants ) ) + ")"
				self.variants.append( var )
			else:
				self.variants.append( var )
		else:
			self.variants.append( var )
		return None

	def finalWriteVariants( self , **kwargs ):
		writeToFile = kwargs.get( 'writeToFile' , "" )
		writeFunction = kwargs.get( 'writeFunction' , None )
		bufferSize = kwargs.get( 'bufferSize' , 1000 )
		outFH = kwargs.get( 'outFH' , None )
		if ( len( self.variants ) > 0 and writeToFile and writeFunction ):
			print "last " + str( len( self.variants ) ) + "; now write"
			for v in self.variants:
				writeFunction( v , **kwargs )
			if ( bufferSize > 0 ):
				self.variants = []
				print "flushed variants (" + str( len( self.variants ) ) + ")"
		if ( type( outFH ) is 'file' ):
			if ( not outFH.closed ):
				outFH.close()
				print "biomine closed"
		return None
	
	def getVCFKeyIndex( self , values , field ):
		if field in self.vcfKeyIndex:
			return values[self.vcfKeyIndex[field]]
		return None

	def setAlleleMeasures( self , var , info , **kwargs ):
		self.setCounts( var , info , **kwargs )
		self.setNumbers( var , info , **kwargs )
		self.setFrequency( var , info , **kwargs )

	@staticmethod
	def setMostSevere( var , mostSevere , mostSevereCons , **kwargs ):
		if mostSevere:
			var.gene = mostSevere.gene
			var.referencePeptide = mostSevere.referencePeptide
			var.positionPeptide = mostSevere.positionPeptide
			var.alternatePeptide = mostSevere.alternatePeptide
			var.transcriptPeptide = mostSevere.transcriptPeptide
			var.transcriptCodon = mostSevere.transcriptCodon
			var.positionCodon = mostSevere.positionCodon
			var.vepVariant.mostSevereConsequence = mostSevereCons
			var.variantClass = mostSevereCons
	
	@staticmethod
	def setCounts( var , info , **kwargs ):
		var.setCounts( info , **kwargs )
	
	@staticmethod
	def setNumbers( var , info , **kwargs ):
		var.setNumbers( info , **kwargs )
	
	@staticmethod
	def setFrequency( var , info , **kwargs ):
		var.setFrequency( info , **kwargs )

	@staticmethod
	def writeField( out , value , last=False ):
		if ( last ):
			out.write( str( value ) + "\n" )
		else:
			out.write( str( value ) + "\t" )

	def searchVCF( self , **kwargs ):
		vcf = kwargs.get( 'vcf' , "" )
		self.variants = kwargs.get( 'queries' , [] )
		entries = {}
		with pysam.VariantFile( vcf , "r" ) as vf:
			for var in self.variants:
				entries[var.genomicVar()] = 1
				here = ':'.join( [ str( var.chromosome ) , str( var.start ) , str( var.stop ) ] )
				for hit in vf.fetch( region = here ):
					for i in range( 0 , len( hit.alts ) ):
						evar = variant( chromosome = hit.contig , \
										start = hit.start , \
										dbsnp = hit.id , \
										reference = hit.alleles[0] , \
										alternate = hit.alts[i] , \
									)
						if var.sameGenomicVariant( evar ):
							entries[var.genomicVar()] =  hit.info["AF"]

