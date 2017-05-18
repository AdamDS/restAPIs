from biomine.variant.mafvariant import mafvariant
#<colocated_variants id="CS114003" allele_string="HGMD_MUTATION" end="7578556" phenotype_or_disease="1" seq_region_name="17" start="7578556" strand="1"/>
#<colocated_variants id="TP53_g.12362A>G" allele_string="T/C" end="7578556" seq_region_name="17" start="7578556" strand="1"/>
#<colocated_variants id="TP53_g.12362A>T" allele_string="T/A" end="7578556" seq_region_name="17" start="7578556" strand="1"/>
#<colocated_variants id="TP53_g.12362del" allele_string="T/-" end="7578556" seq_region_name="17" start="7578556" strand="1"/>
#<colocated_variants id="TP53_g.12362A>C" allele_string="T/G" end="7578556" seq_region_name="17" start="7578556" strand="1"/>
#<colocated_variants id="COSM45672" allele_string="T/C" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM46049" allele_string="T/G" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM45658" allele_string="T/-" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM21585" allele_string="T/A" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM562615" allele_string="T/A" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM318167" allele_string="T/C" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3675526" allele_string="T/G" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM318166" allele_string="T/C" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM562616" allele_string="T/A" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3675529" allele_string="T/G" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3675528" allele_string="T/G" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3388230" allele_string="T/C" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3717681" allele_string="T/A" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM1649367" allele_string="T/A" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM1646843" allele_string="T/C" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3675527" allele_string="T/G" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM318168" allele_string="T/C" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM3675530" allele_string="T/G" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>
#<colocated_variants id="COSM562617" allele_string="T/A" end="7578556" phenotype_or_disease="1" seq_region_name="17" somatic="1" start="7578556" strand="1"/>

class vepcolocatedvariant(mafvariant):
	def __init__(self , **kwargs):
		super(vepcolocatedvariant,self).__init__(**kwargs)
		self.ID = kwargs.get('biotype',"")
		self.phenotypeOrDisease = kwargs.get('exon',None)
		self.somatic= kwargs.get('totalExons',None)
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( vepcolocatedvariant , self ).copyInfo( aParentVariant )

	def printVariant(self,delim , **kwargs ):
		onlyThisVariant = kwargs.get( 'minimal' , False )
		if not onlyThisVariant:
			super(vepcolocatedvariant,self).printVariant( delim , **kwargs )
		print "vepcolocatedvariant: { " ,
		if self.ID:
			print "ID=" ,
			print self.ID + delim ,
		if self.phenotypeOrDisease:
			print "phenotypeOrDisease= " + self.phenotypeOrDisease + delim ,
		if self.somatic:
			print "somatic= " + self.somatic
		print " }"
	def attr(self):
		attributes = super(vepcolocatedvariant,self).attr()
		if self.ID:
			attributes.append(self.ID)
		if self.phenotypeOrDisease:
			attributes.append(self.phenotypeOrDisease)
		if self.somatic:
			attributes.append(self.somatic)
		return attributes

	def parseColocatedVariant( self , colocated , vepvar ):
		''' Expect colocated type as dict from JSON '''
#		print "biomine::variant::vepcolocatedvariant::parseColocatedVariant"
		allele_string = colocated.get( 'allele_string' )
		self.copyInfo( vepvar , 'variant' )
		if allele_string == "HGMD_MUTATION":
			self.ID = colocated.get( 'id' )
			self.chromosome = colocated.get( 'chromosome' )
			self.start = colocated.get( 'start' )
			self.stop = colocated.get( 'stop' )
			self.strand = colocated.get( 'strand' )
		else:
			alleles = allele_string.split('/')
			self.reference = allele_string[0]
			if len( amino_acids ) > 1:
				self.alternatePeptide = amino_acids[1]
			else:
				self.alternatePeptide = amino_acids[0]
