from biomine.variant.mafvariant import mafvariant
from biomine.variant.vepconsequencevariant import vepconsequencevariant
#{
#	"allele_string": "G/A",
#	"assembly_name": "GRCh37",
#	"colocated_variants": [
#		{
#			"allele_string": "G/A",
#			"end": 140534527,
#			"id": "COSM1312758",
#			"phenotype_or_disease": 1,
#			"seq_region_name": "7",
#			"somatic": 1,
#			"start": 140534527,
#			"strand": 1
#		}
#	],
#	"end": 140534527,
#	"id": null,
#	"input": "7 140534527 . G A . . .",
#	"most_severe_consequence": "missense_variant",
#	"seq_region_name": "7",
#	"start": 140534527,
#	"strand": 1,
#	"transcript_consequences": [
#		{
#			"biotype": "retained_intron",
#			"cdna_end": 392,
#			"cdna_start": 392,
#			"consequence_terms": [
#				"non_coding_transcript_exon_variant",
#				"non_coding_transcript_variant"
#			],
#			"exon": "3/3",
#			"gene_id": "ENSG00000157764",
#			"gene_symbol": "BRAF",
#			"gene_symbol_source": "HGNC",
#			"hgnc_id": 1097,
#			"hgvsc": "ENST00000469930.1:n.392C>T",
#			"impact": "MODIFIER",
#			"strand": -1,
#			"transcript_id": "ENST00000469930",
#			"variant_allele": "A"
#		},
#		{
#			"amino_acids": "S/L",
#			"biotype": "protein_coding",
#			"canonical": 1,
#			"ccds": "CCDS5863.1",
#			"cdna_end": 447,
#			"cdna_start": 447,
#			"cds_end": 386,
#			"cds_start": 386,
#			"codons": "tCa/tTa",
#			"consequence_terms": [
#				"missense_variant"
#			],
#			"domains": [
#				{
#					"db": "hmmpanther",
#					"name": "PTHR23257"
#				},
#				{
#					"db": "hmmpanther",
#					"name": "PTHR23257"
#				}
#			],
#			"exon": "3/18",
#			"gene_id": "ENSG00000157764",
#			"gene_symbol": "BRAF",
#			"gene_symbol_source": "HGNC",
#			"hgnc_id": 1097,
#			"hgvsc": "ENST00000288602.6:c.386C>T",
#			"hgvsp": "ENSP00000288602.6:p.Ser129Leu",
#			"impact": "MODERATE",
#			"polyphen_prediction": "benign",
#			"polyphen_score": 0.003,
#			"protein_end": 129,
#			"protein_id": "ENSP00000288602",
#			"protein_start": 129,
#			"refseq_transcript_ids": [
#				"NM_004333.4"
#			],
#			"sift_prediction": "tolerated",
#			"sift_score": 0.2,
#			"strand": -1,
#			"transcript_id": "ENST00000288602",
#			"variant_allele": "A"
#		},
#		{
#			...
#       }
#       ]
#   },

class vepvariant(mafvariant):
	def __init__(self , **kwargs):
		super(vepvariant,self).__init__(**kwargs)
		self.inputVariant = kwargs.get('inputVariant',"")
		self.mostSevereConsequence = kwargs.get('mostSevereConsequence',"")
		self.consequences = kwargs.get('consequences',[])
		self.colocations = kwargs.get('colocations',[])
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( vepvariant , self ).copyInfo( aParentVariant )
	def copyInfo( self , copy ):
		super( vepvariant , self ).copyInfo( copy )
		self.inputVariant = copy.inputVariant
		self.mostSevereConsequence = copy.mostSevereConsequence
		self.consequences = copy.consequences
		self.colocations = copy.colocations
	def fillMissingInfo( self , copy ):
		#print "Variant.vepvariant::fillMissingInfo" ,
		super( vepvariant , self ).fillMissingInfo( copy )
		if not self.inputVariant:
			try:
				self.inputVariant = copy.inputVariant
			except:
				pass
		if not self.mostSevereConsequence:
			try:
				self.mostSevereConsequence = copy.mostSevereConsequence
			except:
				pass
		if not self.consequences:
			try:
				self.consequences = copy.consequences
			except:
				print( "BioMine::variant::vepvariant::fillMissingInfo Warning: no consequences with which to fill" )
				pass
			for consequence in self.consequences:
				if consequence.canonical:
					if consequence.geneSymbolSource == "HGNC":
						super( vepvariant , self ).fillMissingInfo( copy )
					mafvariant.fillMissingInfo( self , consequence )
					#print( str( consequence.terms[0] ) )
					self.variantClass = consequence.terms[0]
		if not self.colocations:
			try:
				self.colocations = copy.colocations
			except:
				pass

	def setInputVariant( self , **kwargs ):
		asVCF = kwargs.get( 'vcf' , True )
		if asVCF:
			self.inputVariant = self.vcf( delim=' ' )
		else:
			self.inputVariant = self.genomicVar()

	def printVariant(self,delim , **kwargs ):
		onlyThisVariant = kwargs.get( 'minimal' , False )
		print( "vepvariant {" )
		if not onlyThisVariant:
			super(vepvariant,self).printVariant( delim , **kwargs )
		print "vepvariant: { " ,
		if self.inputVariant:
			print "inputVariant=" ,
			print self.inputVariant + delim ,
		if self.mostSevereConsequence:
			print "mostSevereConsequence=" , 
			print str( self.mostSevereConsequence ) + delim ,
		if self.consequences:
			print "consequences= [" ,
			for cons in sorted(self.consequences):
				cons.printVariant(delim,**kwargs)
			print "]" + delim ,
		if self.colocations:
			print "colocations= [" ,
			for anno in self.colocations:
				print str(anno) + ", " ,
			print "]" + delim ,
		print " }"
		print( "}" )
	def attr(self):
		attributes = super(vepvariant,self).attr()
		if self.inputVariant:
			attributes.append(self.inputVariant)
		if self.mostSevereConsequence:
			attributes.append(self.mostSevereConsequence)
		if self.consequences:
			attributes.append(self.consequences)
		if self.colocations:
			attributes.append(self.colocations)
		return attributes

	def printConsequencesProteogenomicVar( self ):
		print self.proteogenomicVar()
		for consequence in self.consequences:
			print self.genomicVar() + ", " ,
			print consequence.codingHGVS()

	def parseEntryFromVEP( self , rootElement ):
		''' Expect rootElement as JSON (dict) '''
#		print "biomine::variant::vepvariant::parseEntryFromVEP"
		self.chromosome = rootElement.get( 'seq_region_name' )
		self.start = rootElement.get( 'start' )
		self.stop = rootElement.get( 'end' )
		allele_string = rootElement.get( 'allele_string' ).split('/')
		self.reference = allele_string[0]
		if len( allele_string ) > 1:
			self.alternate = allele_string[1]
		else:
			self.alternate = allele_string[0]
		self.strand = rootElement.get( 'strand' )
		self.assembly = rootElement.get( 'assembly_name' )
		self.mostSevereConsequence = rootElement.get( 'most_severe_consequence' )
		transcriptConsequences = rootElement.get( 'transcript_consequences' )
		self.setTranscriptConsequences( transcriptConsequences )
		colocatedVariants = rootElement.get( 'colocated_variants' )
	def setColocatedVariants( self , colocatedVariants ):
		''' Expect colocatedVariants as dict from JSON '''
		try:
			for colocated in colocatedVariants:
				otherVar = vepcolocatedvariant( parentVariant = self )
				otherVar.parseColocatedVariant( colocated )
				self.colocations.append( otherVar )
		except:
			print( "BioMine Warning: Cannot set colocated variants - no colocations of " + self.genomicVar() )
			pass
			
	def setTranscriptConsequences( self , transcriptConsequences ):
		''' Expect transcriptConsequences as dict from JSON '''
#		print "biomine::variant::vepvariant::setTranscriptConsequence"
		#if not transcriptConsequences:
		#	return
		try:
			for consequence in transcriptConsequences: #list of dict's
				otherVar = vepconsequencevariant( parentVariant=self )
				otherVar.parseTranscriptConsequence( consequence )
				self.consequences.append( otherVar )
				if otherVar.canonical:
					self.gene = otherVar.gene
		except:
			print( "BioMine Warning: Cannot set transcript consequences - no consequence for " + self.genomicVar() )
			pass
