#1-based closed coordinate system
#https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+%28MAF%29+Specification
class variant(object):
	NULL = "."
	def __init__(self , **kwargs):
		self.gene = kwargs.get('gene',"")
		self.chromosome = kwargs.get('chromosome',None)
		if self.chromosome:
			self.cleanChromosome()
		self.start = kwargs.get('start',None)
		self.stop = kwargs.get('stop',None)
		self.reference = kwargs.get('reference',"-")
		self.alternate = kwargs.get('alternate',"-")
		strand = kwargs.get('strand',"+")
		self.setStrand( strand )
		self.sample = kwargs.get('sample',None)
		self.assembly = kwargs.get('assembly',None)
		self.dbsnp = kwargs.get('dbsnp',None)
	#def setReference( self , start ):
		#if str(reference) == "0" or not start:
			#self.start = "-"
	def setStrand( self , strand ):
		if strand == -1 or strand == "-":
			self.strand = "-"
		else:
			self.strand = "+"
#		print self.strand
	
	def getReference( self ):
		if self.reference == variant.NULL:
			return variant.NULL
		return self.reference
	def checkIfRefAltStrand( self , k ):
		if ( k == "reference" ):
			return True
		if ( k == "alternate" ):
			return True
		if ( k == "strand" ):
			return True
		return False

	def nonzeroRefAltStrand( self , k ):
		if ( k == "reference" ):
			if ( self.__dict__.get( k ) != "-" ):
				return True
		elif ( k == "alternate" ):
			if ( self.__dict__.get( k ) != "-" ):
				return True
		elif ( k == "strand" ):
			if ( self.__dict__.get( k ) != "+" ):
				return True
		return False
	def __nonzero__( self ):
		for k , v in self.__dict__.iteritems():
			if ( self.checkIfRefAltStrand( k ) ):
				if ( self.nonzeroRefAltStrand( k ) ):
					return True
			else:
				if ( bool( v ) ):
					return True
		return False
		
	def copyInfo( self , copy ):
		self.gene = copy.gene
		self.chromosome = copy.chromosome
		self.start = copy.start
		self.stop = copy.stop
		self.reference = copy.reference
		self.alternate = copy.alternate
		self.setStrand( copy.strand )
		self.sample = copy.sample
		self.assembly = copy.assembly
		self.dbsnp = copy.dbsnp
	def fillMissingInfo( self , copy ):
		#print "Variant.variant::fillMissingInfo" ,
		if not self.gene:
			self.gene = copy.gene
		if not self.chromosome:
			self.chromosome = copy.chromosome
		if not self.start:
			self.start = copy.start
		if not self.stop:
			self.stop = copy.stop
		if not self.reference or ( self.reference == "-" and copy.reference != "-" ):
			self.reference = copy.reference
		if not self.alternate or ( self.alternate == "-" and copy.alternate != "-" ):
			self.alternate = copy.alternate
		if not self.strand:
			self.setStrand( copy.strand )
		if not self.sample:
			self.sample = copy.sample
		if not self.assembly:
			self.assembly = copy.assembly
		if not self.dbsnp:
			self.dbsnp = copy.dbsnp
			
	def printVariant(self,delim , **kwargs ):
		print "variant { " ,
		if self.gene:
			print "gene= " ,
			print self.gene + delim ,
		if self.chromosome:
			print "chromosome= " ,
			print str(self.chromosome) + delim ,
		if self.start:
			print "start= " ,
			print str(self.start) + delim ,
		if self.stop:
			print "stop= " ,
			print str(self.stop) + delim ,
		if self.reference:
			print "reference= " ,
			print self.reference + delim ,
		if self.alternate:
			print "alternate= " ,
			print self.alternate + delim ,
		if self.strand:
			print "strand= " ,
			print str(self.strand) + delim ,
		if self.dbsnp:
			print "dbsnp= " ,
			print "rs" + self.dbsnp + delim ,
		if self.sample:
			print "sample= " ,
			print self.sample + delim ,
		if self.assembly:
			print "assembly= " ,
			print self.assembly + delim ,
		print " }"
	def attr(self):
		attributes = []
		if self.gene:
			attributes.append(self.gene)
		if self.chromosome:
			attributes.append(self.chromosome)
		if self.start:
			attributes.append(self.start)
		if self.stop:
			attributes.append(self.stop)
		if self.reference:
			attributes.append(self.reference)
		if self.alternate:
			attributes.append(self.alternate)
		if self.strand:
			attributes.append(self.strand)
		if self.dbsnp:
			attributes.append(self.dbsnp)
		if self.sample:
			attributes.append(self.sample)
		return attributes
	def genomicVar( self ):
		out = str(self.gene) + ":" \
		+ str(self.chromosome) + ":" \
		+ str(self.start) + "-" \
		+ str(self.stop) \
		+ str(self.reference) \
		+ ">" + str(self.alternate)
		return out
	def HGVSg( self ):
		out = str(self.chromosome) + ":g." \
			+ str(self.start)
		if self.reference != "-" and self.alternate != "-":
			if self.start == self.stop: #snp
				out += str(self.reference) + ">" \
				+ str(self.alternate)
			else: #complex
				out += "_" + str(self.stop) \
				+ "del" + str(self.reference) \
				+ "ins" + str(self.alternate)
		elif self.reference == "-" and self.alternate != "-":
			out += "_" + str(self.stop) \
			+ "ins" + str(self.alternate)
		elif self.reference != "-" and self.alternate == "-":
			out += "_" + str(self.stop) \
			+ "del" + str(self.reference)
		return out
	def vcf( self , **kwargs ):
		delim = kwargs.get( 'delim' , ' ' )
		null = kwargs.get( 'null' , "." )
		ref = self.reference
		if ref == "-":
			ref = null
		alt = self.alternate
		if alt == "-":
			alt = null
		dbsnp = self.dbsnp
		if not dbsnp:
			dbsnp = null
		return delim.join( [ self.chromosome , \
							str( self.start ) , \
							dbsnp , ref , alt , \
							null , null , null ] )
	def ensembl( self , **kwargs ):
		#http://useast.ensembl.org/info/docs/tools/vep/vep_formats.html#input
		delim = kwargs.get( 'delim' , ' ' )
		null = kwargs.get( 'null' , "-" )
		ref = self.reference
		alt = self.alternate
		lenRef = len( self.reference )
		lenAlt = len( self.alternate )
		stop = self.start
		if ref == "-": #insertion
			stop = int( self.start ) - 1
		elif alt == "-": #deletion
			stop = int( self.start ) + lenRef
		refalt = ref + "/" + alt
		return delim.join( [ self.chromosome , str( self.start ) , str( stop ) , refalt , str( self.strand ) ] )
		
	def region( self ):
		return str( self.chromosome ) + ":" \
			+ str( self.start ) + ".." \
			+ str( self.stop ) + ":" \
			+ str( self.strand )
	def vcfLine2Variant( self , record , **kwargs ):
		#http://pyvcf.readthedocs.org/en/latest/api.html#vcf-model-record
		self.chromosome = record.CHROM
		self.cleanChromosome()
		#start/stop should be 0-base, half-open [)
		self.start = record.POS
		self.reference = record.REF
		alternates = record.ALT
		if len(alternates) > 1:
			self.alternate = alternates[0] #TODO should get ALL alternates
		self.dbsnp = record.ID
		self.stop = self.start #assume SNP
		if self.reference == "-": #insertion
			self.stop += 1
		else: #deletion, xNP, complex
			if len( self.reference ) > 1:
				self.stop = len( self.reference ) + self.start - 1
	def cleanChromosome( self ):
		''' Get the chromosome number in case chr or Chr is present'''
		if type( self.chromosome ) is str:
			chrom = self.chromosome.lower()
			clean = chrom.replace( "chr" , "" )
			self.chromosome = clean.upper()
	def mafLine2Variant( self , line , **kwargs ):
##		print "variant::mafLine2Variant - " ,
		fields = line.split( "\t" )
		self.gene = fields[0]	#1	Hugo_Symbol
		self.assembly = fields[3] #4	NCBI_Build
		self.chromosome = fields[4]	#5	Chromosome
		self.start = fields[5]	#6	Start_Position
		self.stop = fields[6]	#7	End_Position
		self.setStrand( fields[7] )	#8	Strand
		self.reference = fields[10]	#11	Reference_Allele
		if str(self.reference) == "0" or not self.reference:
			self.reference = "-"
		self.alternate = fields[11] if fields[11] != fields[10] else fields[12]	#12	Tumor_Seq_Allele1	#13	Tumor_Seq_Allele2
		if str(self.alternate) == "0" or not self.alternate:
			self.alternate = "-"
		self.dbsnp = fields[13]
		self.sample = fields[15]
	def uniqueVar( self ):
		if self.sample:
			return self.sample + "::" + self.genomicVar()
		else:
			return "nosample::" + self.genomicVar()


	def setStopFromReferenceAndAlternate( self ):
		if len( self.reference ) > 1 or len( self.alternate ) > 1 \
		or self.reference == "-" or self.alternate == "-" \
		or self.reference == variant.NULL or self.alternate == variant.NULL:
			if self.reference == "-" or self.reference == variant.NULL:
				self.stop = self.start
			else:
				self.stop = int( self.start ) + len( self.alternate )
				if self.reference[0] == self.alternate[0]:
					self.stop += 1
		else:
			self.stop = self.start

	def determineLengthOfIndel( self ):
		lengthOfIndel = 0
		endDeletion = 0
		if self.isIndel():
			lenRef = len(self.reference)
			lenAlt = len(self.alternate)
			lengthOfIndel = lenAlt - lenRef
#           if lenRef > 1:
#               print "is del"
		return lengthOfIndel

	def isIndel( self ):
		if self.variantType == "INS" or self.variantType == "DEL":
			return True
		return False

	def compareVariants( self , otherVariant ):
		common = 0
		attributes = self.attr()
		otherAttributes = otherVariant.attr()
		for attr in attributes:
			for otherAttr in otherAttributes:
				if attr == otherAttr:
					common += 1
		percentMatch = common / len( attributes )
		return percentMatch

	def readVariants(self,inputFile):
		variants = []
		if inputFile:
			inFile = open( inputFile , 'r' )
			for line in inFile:
				fields = line.split( '\t' )
				variants.append( fields[0] + ":" + fields[1] )
		return variants
	def reset( self ):
		self = variant.__init()

	def sameGenomicVariant( self , otherVar ):
#		print "sameGenomicVariant - " ,
#		print str( self.alternate ) + " vs " + str( otherVar.alternate ) ,
		if self.sameGenomicReference( otherVar ):
			if otherVar.alternate == self.alternate:
#				print "True"
				return True
#		print "False"
		return False
	def sameGenomicReference( self , otherVar ):
#		print "sameGenomicReference - " ,
#		print str( self.reference ) + " vs " + str( otherVar.reference ) ,
		if self.sameGenomicPosition( otherVar ):
			if otherVar.reference == self.reference:
#				print "True"
				return True
#		print "False"
		return False
	def sameGenomicPosition( self , otherVar ):
#		print "sameGenomicPosition - " ,
#		print str( self.chromosome ) + ":" + str( self.start ) + " vs " + str( otherVar.chromosome ) + ":" + str( otherVar.start ) ,
		if otherVar.chromosome == self.chromosome and \
			str(otherVar.start) == str(self.start):
#			print "True"
			return True
#		print "False"
		return False
	def removeOverlapFromReferenceAndAlternate( self ):
		updatedRef = ""
		updatedAlt = ""
		loverlap = 0
		ref = self.reference
		alt = self.alternate
		for i in range( 0 , min( [ len( ref ) , len( alt ) ] ) ):
			if ref[i] == alt[i]:
				loverlap = i
				next
			else:
				break
		if loverlap < len( ref ) - 1:
			updatedRef = ref[loverlap+1:]
		if loverlap < len( alt ) - 1:
			updatedAlt = alt[loverlap+1:]
		self.setStartStopForOverlapFix( loverlap , \
											 updatedRef , updatedAlt )
		self.reference = self.nullCheck( updatedRef )
		self.alternate = self.nullCheck( updatedAlt )
	def setStartStopForOverlapFix( self , loverlap , ref , alt ):
		start = self.start + loverlap + 1
		stop = self.start
		if len( ref ) < len( alt ) and len( ref ) == 0: #insertion
			start -= 1
			stop = start + 1
		elif len( ref ) > len( alt ) and len( alt ) == 0: #deletion
			stop = start + len( ref ) - 1
		elif len( ref ) > 0 and len( alt ) > 0: #complex
			stop = start + len( alt )
		else: #dunno
			print( "biomine::variant::setStartStopForOverlapFix warning:" \
			+ " Unusual ref alt start stop: " + ref + " " \
			+ alt + " " + str( self.start ) + " " + str( self.stop ) \
			+ ", reference and alternate not cleaned" )
		self.start = start
		self.stop = stop
	@staticmethod
	def nullCheck( val ):
		if not val:
			val = "-"
		return val
'''
mu = variant(gene="BRAF",chromosome=7,start=12345,stop=123456,reference="AT",alternate="GC")
mu.printVariant('\t')
ou = variant(gene="BRAF",chromosome=7,start=12345,stop=123456,reference="AT",alternate="GC")
ou.printVariant('\t')
pu = variant(chromosome=7,start=12345,stop=123456,reference="AT",alternate="GC")
pu.printVariant('\t')
qu = variant(chromosome=7,start=12345,reference="AT",alternate="GC")
qu.printVariant('\t')

print qu.attr()
'''
