from biomine.variant.mafvariant import mafvariant
#{
#	"biotype": "retained_intron",
#	"cdna_end": 392,
#	"cdna_start": 392,
#	"consequence_terms": [
#		"non_coding_transcript_exon_variant",
#		"non_coding_transcript_variant"
#	],
#	"exon": "3/3",
#	"gene_id": "ENSG00000157764",
#	"gene_symbol": "BRAF",
#	"gene_symbol_source": "HGNC",
#	"hgnc_id": 1097,
#	"hgvsc": "ENST00000469930.1:n.392C>T",
#	"impact": "MODIFIER",
#	"strand": -1,
#	"transcript_id": "ENST00000469930",
#	"variant_allele": "A"
#},
#{
#   "amino_acids": "S/L",
#   "biotype": "protein_coding",
#   "canonical": 1,
#   "ccds": "CCDS5863.1",
#   "cdna_end": 447,
#   "cdna_start": 447,
#   "cds_end": 386,
#   "cds_start": 386,
#   "codons": "tCa/tTa",
#   "consequence_terms": [
#	   "missense_variant"
#   ],
#   "domains": [
#	   {
#		   "db": "hmmpanther",
#		   "name": "PTHR23257"
#	   },
#	   {
#		   "db": "hmmpanther",
#		   "name": "PTHR23257"
#	   }
#   ],
#   "exon": "3/18",
#   "gene_id": "ENSG00000157764",
#   "gene_symbol": "BRAF",
#   "gene_symbol_source": "HGNC",
#   "hgnc_id": 1097,
#   "hgvsc": "ENST00000288602.6:c.386C>T",
#   "hgvsp": "ENSP00000288602.6:p.Ser129Leu",
#   "impact": "MODERATE",
#   "polyphen_prediction": "benign",
#   "polyphen_score": 0.003,
#   "protein_end": 129,
#   "protein_id": "ENSP00000288602",
#   "protein_start": 129,
#   "refseq_transcript_ids": [
#	   "NM_004333.4"
#   ],
#   "sift_prediction": "tolerated",
#   "sift_score": 0.2,
#   "strand": -1,
#   "transcript_id": "ENST00000288602",
#   "variant_allele": "A"
#},

class vepconsequencevariant(mafvariant):
	def __init__(self , **kwargs):
		super(vepconsequencevariant,self).__init__(**kwargs)
		self.biotype = kwargs.get('biotype',"")
		self.terms = kwargs.get('consequence_terms',[])
		self.exon = kwargs.get('exon',None)
		self.totalExons = kwargs.get('totalExons',None)
		self.intron = kwargs.get('intron',None)
		self.totalIntron = kwargs.get('totalIntrons',None)
		self.ensg = kwargs.get('gene_id',"")
		self.geneSymbolSource = kwargs.get('gene_symbol_source',"")
		self.hgnc = kwargs.get('hgnc_id',"")
		self.impact = kwargs.get('impact',"")
		self.positionCDS = kwargs.get('positionCDS',None)
		self.referenceCodons = kwargs.get('referenceCodons',"")
		self.alternateCodons = kwargs.get('alternateCodons',"")
		self.domains = kwargs.get('domains',{})
		self.predictionPolyphen = kwargs.get('predictionPolyphen',"")
		self.scorePolyphen = kwargs.get('scorePolyphen',None)
		self.predictionSIFT = kwargs.get('predictionSIFT',"")
		self.scoreSIFT = kwargs.get('scoreSIFT',None)
		### Optional annotations
		self.blosum = kwargs.get( 'blosum' , None )
		self.csn = kwargs.get( 'csn' , None )
		self.compara = kwargs.get( 'compara' , None )
		self.exac = kwargs.get( 'exac' , None )
		self.genesplicer = kwargs.get( 'genesplicer' , None )
		self.maxentscan = kwargs.get( 'maxentscan' , [] )
		#self.updown = kwargs.get( 'updown' , None )
		#self.callback = kwargs.get( 'callback' , None )
		self.canonical = kwargs.get( 'canonical' , False )
		self.ccds = kwargs.get( 'ccds' , None )
		self.dbnsfp = kwargs.get( 'dbnsfp' , None )
		self.dbscsnv = kwargs.get( 'dbscsnv' , None )
		#self.mirna = kwargs.get( 'mirna' , None )
		#self.numbers = kwargs.get( 'numbers' , None )
		#self.protein = kwargs.get( 'protein' , None )
		self.transcriptsRefSeq = kwargs.get('transcriptsRefSeq',[])
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( vepconsequencevariant , self ).copyInfo( aParentVariant )

	def printVariant(self,delim , **kwargs ):
		onlyThisVariant = kwargs.get( 'minimal' , False )
		if not onlyThisVariant:
			super(vepconsequencevariant,self).printVariant( delim , **kwargs )
		print "vepconsequencevariant: { " ,
		if self.biotype:
			print "biotype=" ,
			print self.biotype + delim ,
		if self.terms:
			print "terms= [" ,
			for cons in self.terms:
				print cons + delim ,
			print "]" + delim ,
		if self.exon:
			print "exon= " + self.exon + delim ,
		if self.totalExons:
			print "totalExons= " + self.totalExons + delim ,
		if self.intron:
			print "intron= " + self.intron+ delim ,
		if self.ensg:
			print "Ensembl_gene_id= " + self.ensg + delim ,
		if self.geneSymbolSource:
			print "gene_symbol_source= " + self.geneSymbolSource + delim ,
		if self.hgnc:
			print "HGNC_ID= " + str( self.hgnc ) + delim ,
		if self.impact:
			print "impact= " + self.impact + delim ,
		if self.maxentscan:
			print "Max_Ent_Scan[ref,alt,diff]= [ " ,
			for val in self.maxentscan:
				print str( val ) , "," ,
			print " ]" + delim ,
		print "canonical= " + str( self.canonical ) + delim ,
		if self.positionCDS:
			print "positionCDS= " + str(self.positionCDS) + delim ,
		if self.referenceCodons:
			print "referenceCodons= " + self.referenceCodons + delim ,
		if self.alternateCodons:
			print "alternateCodons= " + self.alternateCodons + delim ,
		if self.domains:
			print "domains= {" ,
			for domain in self.domains.keys():
				print str(domain) ,
				print " : " ,
				print str(self.domains[domain]) ,
				print ", " ,
			print "}" + delim ,
		if self.predictionPolyphen:
			print "predictionPolyphen= " + self.predictionPolyphen + delim ,
		if self.scorePolyphen:
			print "scorePolyphen= " + str(self.scorePolyphen) + delim ,
		if self.transcriptsRefSeq:
			print "transcriptsRefSeq= [" ,
			for ids in self.transcriptsRefSeq:
				print ids + delim ,
			print "]" + delim ,
		if self.predictionSIFT:
			print "predictionSIFT= " + self.predictionSIFT + delim ,
		if self.scoreSIFT:
			print "scoreSIFT= " + str(self.scoreSIFT) + delim ,
		print " }"
	def attr(self):
		attributes = super(vepconsequencevariant,self).attr()
		if self.biotype:
			attributes.append(self.biotype)
		if self.terms:
			attributes.append(self.terms)
		if self.exon:
			attributes.append(self.exon)
		if self.totalExons:
			attributes.append(self.totalExons)
		if self.intron:
			attributes.append(self.intron)
		if self.totalIntrons:
			attributes.append(self.totalIntrons)
		if self.ensg:
			attributes.append(self.ensg)
		if self.impact:
			attributes.append(self.impact)
		return attributes

	def parseTranscriptConsequence( self , consequence ):
		''' Expect consequence type as dict from JSON '''
		#print "biomine::variant::vepconsequencevariant::parseTranscriptConsequence"
		if "amino_acids" in consequence:
			amino_acids = consequence.get( 'amino_acids' ).split('/')
			self.referencePeptide = amino_acids[0]
			if len( amino_acids ) > 1:
				self.alternatePeptide = amino_acids[1]
			else:
				self.alternatePeptide = amino_acids[0]
		if "biotype" in consequence:
			self.biotype = consequence.get( 'biotype' )
		if "canonical" in consequence:
			canonical = consequence.get( 'canonical' , False )
			if canonical:
				self.canonical = True
		if "ccds" in consequence:
			self.ccds = consequence.get( 'ccds' )
		if "cdna_start" in consequence:
			self.positionCodon = consequence.get( 'cdna_start' )
		if "cds_start" in consequence:
			self.positionCDS = consequence.get( 'cds_start' )
		if "codons" in consequence:
			codons = consequence.get( 'codons' ).split('/')
			self.referenceCodons = codons[0]
			if len( codons ) > 1:
				self.alternateCodons = codons[1]
			else:
				self.alternateCodons = codons[0]
		if "conservation" in consequence:
			self.compara = consequence.get( 'conservation' )
		if "maxentscan_ref" in consequence:
			self.maxentscan.append( consequence.get( 'maxentscan_ref' ) )
		if "maxentscan_alt" in consequence:
			self.maxentscan.append( consequence.get( 'maxentscan_alt' ) )
		if "maxentscan_diff" in consequence:
			self.maxentscan.append( consequence.get( 'maxentscan_diff' ) )
		if "gene_splicer" in consequence:						#unsure about string
			self.genesplicer = consequence.get( 'gene_splicer' )
		if "consequence_terms" in consequence:
			terms = consequence.get( 'consequence_terms' )
			for term in terms:
				self.terms.append( term )
		if "domains" in consequence:
			domains = consequence.get( 'domains' )
			for domain in domains:
				self.domains[domain["db"]] = domain["name"]
		if "exon" in consequence:
			exonOfExons = consequence.get( 'exon' ).split('/')
			self.exon = exonOfExons[0]
			if len( exonOfExons ) > 1:
				self.totalExons = exonOfExons[1]
			else:
				self.totalExons = exonOfExons[0]
		if "intron" in consequence:
			intronOfIntrons = consequence.get( 'intron' ).split('/')
			self.intron = intronOfIntrons[0]
			if len( intronOfIntrons ) > 1:
				self.totalIntrons = intronOfIntrons[1]
			else:
				self.totalIntrons = intronOfIntrons[0]
		if "gene_id" in consequence:
			self.ensg = consequence.get( 'gene_id' )
		if "gene_symbol" in consequence:
			self.gene = consequence.get( 'gene_symbol' )
		if "gene_symbol_source" in consequence:
			self.geneSymbolSource = consequence.get( 'gene_symbol_source' )
		if "hgnc_id" in consequence:
			self.hgnc = consequence.get( 'hgnc_id' )
#		if "hgvsc" in consequence:
#			hgvsc = consequence.get( 'hgvsc' )
#			[ self.reference , self.positionCodon , self.alternate ] = self.splitHGVSc( hgvsc , xDot="[nc]\." )
#		if "hgvsp" in consequence:
#			hgvsp = consequence.get( 'hgvsp' )
#			[ self.referencePeptide , self.positionPeptide , self.alternatePeptide ] = self.splitHGVSp( hgvsp )
		if "impact" in consequence:
			self.impact = consequence.get( 'impact' )
		if "polyphen_prediction" in consequence:
			self.predictionPolyphen = consequence.get( 'polyphen_prediction' )
		if "polyphen_score" in consequence:
			self.scorePolyphen = consequence.get( 'polyphen_score' )
		if "protein_id" in consequence:
			self.transcriptPeptide = consequence.get( 'protein_id' )
		if "protein_start" in consequence:
			self.positionPeptide = consequence.get( 'protein_start' )
		if "refseq_transcript_ids" in consequence:
			refseqIDs = consequence.get( 'refseq_transcript_ids' )
			for refseqID in refseqIDs:
				self.transcriptsRefSeq.append( refseqID )
		if "sift_prediction" in consequence:
			self.predictionSIFT = consequence.get( 'sift_prediction' )
		if "sift_score" in consequence:
			self.scoreSIFT = consequence.get( 'sift_score' )
		if "strand" in consequence:
			self.strand = consequence.get( 'strand' )
		if "transcript_id" in consequence:
			self.transcriptCodon = consequence.get( 'transcript_id' )
		if "variant_allele" in consequence:
			self.alternate = consequence.get( 'variant_allele' )
