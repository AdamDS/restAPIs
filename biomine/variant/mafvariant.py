import re
from biomine.variant.variant import variant
from urllib2 import unquote

class mafvariant(variant):
	nonCoding = re.compile( "([NULL|\*|\+|\-])" )
	upstream = re.compile( "(\*)" )
	intronic = re.compile( "([\+|\-])" )
	inversion = re.compile( "(inv)" )
	duplication = re.compile( "(dup)" )
	deletion = re.compile( "(del)" )
	insertion = re.compile( "(ins)" )
	indel = re.compile( "(del.*ins)" )
	multiple = re.compile( "(_)" )
	varSymbols = { "fs" : "fs" , "x" : "*" , "ter" : "*" , "=" : "" }
	revVarSymbols = dict( ( v , k ) for k , v in varSymbols.iteritems() )
	toShort = { "ala" : "A" , "arg" : "R" , "asn" : "N" , "asp" : "D" , \
				"cys" : "C" , "gln" : "Q" , "glu" : "E" , "gly" : "G" , \
				"his" : "H" , "ile" : "I" , "leu" : "L" , "lys" : "K" , \
				"met" : "M" , "pro" : "P" , "phe" : "F" , "ser" : "S" , \
				"thr" : "T" , "trp" : "W" , "tyr" : "Y" , "val" : "V" }
	toLong = dict( ( v , k ) for k , v in toShort.iteritems() )
	def __init__(self , **kwargs):
		super(mafvariant,self).__init__(**kwargs)
		self.referencePeptide = kwargs.get('referencePeptide',"")
		self.positionPeptide = kwargs.get('positionPeptide',None)
		self.alternatePeptide = kwargs.get('alternatePeptide',"")
		self.transcriptPeptide = kwargs.get('transcriptPeptide',None)
		self.positionCodon = kwargs.get('positionCodon',None)
		self.transcriptCodon = kwargs.get('transcriptCodon',None)
		self.variantClass = kwargs.get('variantClass',None)
		self.variantType = kwargs.get('variantType',None)
		self.disease = kwargs.get('disease',None)
		aParentVariant = kwargs.get( 'parentVariant' , None )
		if aParentVariant:
			super( mafvariant , self ).copyInfo( aParentVariant )
	def copyInfo( self , copy ):
		super( mafvariant , self ).copyInfo( copy )
		try:
			self.referencePeptide = copy.referencePeptide
		except:
			pass
		try:
			self.positionPeptide = copy.positionPeptide
		except:
			pass
		try:
			self.alternatePeptide = copy.alternatePeptide
		except:
			pass
		try:
			self.transcriptPeptide = copy.transcriptPeptide
		except:
			pass
		try:
			self.positionCodon = copy.positionCodon
		except:
			pass
		try:
			self.transcriptCodon = copy.transcriptCodon
		except:
			pass
		try:
			self.variantClass = copy.variantClass
		except:
			pass
		try:
			self.variantType = copy.variantType
		except:
			pass
		try:
			self.disease = copy.disease
		except:
			pass
	def fillMissingInfo( self , copy ):
		#print "Variant.mafvariant::fillMissingInfo" ,
		super( mafvariant , self ).fillMissingInfo( copy )
		if not self.referencePeptide:
			try:
				self.referencePeptide = copy.referencePeptide
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no referencePeptide with which to fill" )
		if not self.positionPeptide:
			try:
				self.positionPeptide = copy.positionPeptide
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no positionPeptide with which to fill" )
		if not self.alternatePeptide:
			try:
				self.alternatePeptide = copy.alternatePeptide
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no alternatePeptide with which to fill" )
		if not self.transcriptPeptide:
			try:
				self.transcriptPeptide = copy.transcriptPeptide
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no transcriptPeptide with which to fill" )
		if not self.positionCodon:
			try:
				self.positionCodon = copy.positionCodon
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no positionCodon with which to fill" )
		if not self.transcriptCodon:
			try:
				self.transcriptCodon = copy.transcriptCodon
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no transcriptCodon with which to fill" )
		if not self.variantClass:
			try:
				self.variantClass = copy.variantClass
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no variantClass with which to fill" )
		if not self.variantType:
			try:
				self.variantType = copy.variantType
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no variantType with which to fill" )
		if not self.disease:
			try:
				self.disease = copy.disease
			except:
				print( "BioMine::variant::mafvariant::fillMissingInfo Warning: no disease with which to fill" )
	def __nonzero__( self ):
		for k , v in self.__dict__.iteritems():
			if ( self.checkIfRefAltStrand( k ) ):
				if ( self.nonzeroRefAltStrand( k ) ):
					return True
			else:
				if ( bool( v ) ):
					return True
		return False

	def printVariant(self,delim , **kwargs ):
		onlyThisVariant = kwargs.get( 'minimal' , False )
		print( "mafvariant {" )
		if not onlyThisVariant:
			super(mafvariant,self).printVariant(delim , **kwargs )
			print "mafvariant: { " ,
			if self.referencePeptide:
				print "referencePeptide= " ,
				print self.referencePeptide + delim ,
			if self.positionPeptide:
				print "positionPeptide= " ,
				print str(self.positionPeptide) + delim ,
			if self.alternatePeptide:
				print "alternatePeptide= " ,
				print self.alternatePeptide + delim ,
			if self.transcriptPeptide:
				print "transcriptPeptide= " ,
				print self.transcriptPeptide + delim ,
			if self.positionCodon:
				print "positionCodon= " ,
				print str(self.positionCodon) + delim ,
			if self.transcriptCodon:
				print "transcriptCodon= " ,
				print self.transcriptCodon + delim ,
			if self.variantClass:
				print "variantClass= " ,
				print self.variantClass + delim ,
			if self.variantType:
				print "variantType= " ,
				print self.variantType + delim ,
			if self.disease:
				print "disease= " ,
				print self.disease + delim ,
			print " }"
		print( "}" )
	def attr(self):
		attributes = []
		super(mafvariant,self).attr()
		if self.positionPeptide:
			attributes.append(self.positionPeptide)
		if self.referencePeptide:
			attributes.append(self.referencePeptide)
		if self.alternatePeptide:
			attributes.append(self.alternatePeptide)
		if self.transcriptPeptide:
			attributes.append(self.transcriptPeptide)
		if self.positionCodon:
			attributes.append(self.positionCodon)
		if self.transcriptCodon:
			attributes.append(self.transcriptCodon)
		if self.disease:
			attributes.append(self.disease)
		return attributes
			
#	def fromHGVSc( self , hgvsc ):
#		pass

	def mafLine2Variant( self , line , **kwargs ):
#		print "biomine::variant::mafvariant::mafLine2Variant - " , 
		super(mafvariant,self).mafLine2Variant( line , **kwargs )
		codonColumn = kwargs.get( 'codon' , 48 )
		deltaPeptideColumn = kwargs.get( 'peptideChange' , 49 )
		fields = line.split( "\t" )
		self.variantClass = fields[8]	#9	Variant_Classification
		self.variantType = fields[9]	#10	Variant_Type
		self.splitHGVSc( fields[int(codonColumn)] )
		self.splitHGVSp( fields[int(deltaPeptideColumn)] ) #################################### Custom field, not reliable in general

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

### Peptide
	def samePeptideReference( self , otherVar ):
#		print "samePeptideReference - " ,
		if otherVar.referencePeptide == self.referencePeptide and \
			otherVar.positionPeptide == self.positionPeptide: #same genomic position & reference
				return True
		return False
	def samePeptideChange( self , otherVar ):
#		print "samePeptideChange - " ,
		if self.samePeptideReference( otherVar ):
			if otherVar.alternatePeptide == self.alternatePeptide:
				return True
		return False
	def hgvspIsIndel( self , hgvsp ):
		pattern = re.compile( "([a-zA-Z]+?)([0-9]+?)([a-zA-Z\*]+)" )
		change = pattern.match( hgvsp )
		if change:
			fsPattern = re.compile( "(fs)" )
			fsMatch = fsPattern.match( hgvsp )
			if fsMatch:
				return 2
			else:
				return 1
		return 0
	def HGVSp( self ):
#		print "biomine::variant::mafvariant::HGVSp"
		hgvsp = ""
		if self.transcriptPeptide:
			hgvsp += str(self.transcriptPeptide)
		elif self.gene:
			hgvsp += str(self.gene)
		hgvsp += ":p."
		if self.referencePeptide:
			hgvsp += str(self.referencePeptide)
		if self.positionPeptide:
			hgvsp += str(self.positionPeptide)
		if self.alternatePeptide:
			hgvsp += str(self.alternatePeptide)
		return hgvsp
	def splitHGVSp( self , hgvsp ):
##		print "biomine::variant::mafvariant::splitHGVSp - "
		hgvsp=unquote(hgvsp)
		ref = ""
		pos = ""
		mut = ""
		if not hgvsp:
			return [ ref , pos , mut ]
		shgvsp = hgvsp.split( ":" )
		#print( hgvsp )
		if len( shgvsp ) > 1:
			self.transcriptPeptide = shgvsp[0]
			hgvsp = shgvsp[1]
		else:
			hgvsp = shgvsp[0]
		#print( str( self.transcriptPeptide ) + "  :  " + hgvsp )
		isNon = self.hgvspIsNonCoding( hgvsp )
		if not isNon:
#TODO handle splice variants: NM_000030.2:c.424-2A>G  NP_000021.1:p.Gly_142Gln145del
			changep = re.match( "p\.([a-zA-Z\*]+?)([0-9\?]+?)([a-zA-Z]{1,3}|[\*\?])(ext[0-9\*\?]*)*" , hgvsp )
			changee = re.match( "(e)([0-9]+?)([\+\-][0-9]+?)" , hgvsp )
			unknown = re.match( "p\.[\?\=\|(\=\)|0|0\?]" , hgvsp )
			if changep: #peptide
				ref = changep.group( 1 )
				pos = changep.group( 2 )
				mut = changep.group( 3 )
				ref = self.convertAA( ref )
				mut = self.convertAA( mut )
				if ( changep.group( 4 ) ):
					mut += changep.group( 4 )
				self.referencePeptide = ref
				self.positionPeptide = pos
				self.alternatePeptide = mut
			elif changee: #intronic
				#print( changee.groups() )
				changes = changee.groups()
				ref = changee.group( 1 )
				pos = changee.group( 2 )
				mut = changee.group( 3 )
				ref = self.convertAA( ref )
				mut = self.convertAA( mut )
				self.referencePeptide = ref
				self.positionPeptide = pos
				self.alternatePeptide = mut
			elif unknown:
				pass
			else:
				print "biomine::variant::mafvariant Warning: could not find amino acid change or intronic change"
				print "  Hint: Is the input amino acid change column correct?"
				print "    Problem variant: " ,
				print( self.proteogenomicVar() + "  --  " + hgvsp )
		else:
			parts = hgvsp.split('\.')
			pos = parts[-1]
		return [ ref , pos , mut ]
	def convertAA( self , pep ):
		#print "mafvariant::convertAA - " + pep
		lowPep = pep.lower()
		pattern = re.compile( "(fs)" )
		fs = pattern.match( lowPep )
		if ( lowPep in mafvariant.toShort ):
			return mafvariant.toShort[lowPep]
		if ( lowPep in mafvariant.varSymbols ):
			return mafvariant.varSymbols[lowPep]
		if ( pep not in mafvariant.toLong and pep not in mafvariant.revVarSymbols ):
			print( "biomine warning: " + str( pep ) + " not found in conversion tables" )
		return pep
	def hgvspIsNonCoding( self , hgvsp ):
#		print "biomine::variant::mafvariant::hgvspIsNonCoding - " ,
		pattern = re.compile( "(NULL)" )
		match = pattern.match( hgvsp )
		if match:
			return True
		return False

### Codon
	def hgvscIsNonCoding( self , hgvsc ):
#		print "biomine::variant::mafvariant::hgvscIsNonCoding - " ,
		match = mafvariant.nonCoding.match( hgvsc )
		if match:
			return True
		return False

	def hgvscIsIndel( self , hgvsc ):
#		print "biomine::variant::mafvariant::hgvscIsIndel - " ,
		match = mafvariant.indel.match( hgvsc )
		if match: #then its a complex indel
			return 2
		match = mafvariant.insertion.match( hgvsc )
		if match:
			return 1
		match = mafvariant.deletion.match( hgvsc )
		if match:
			return -1
		return 0

	def hgvscIsInversion( self , hgvsc ):
		match = mafvariant.inversion.match( hgvsc )
		if match:
			return True
		return False

	def hgvscIsDuplication( self , hgvsc ):
		match = mafvariant.duplication.match( hgvsc )
		if match:
			return True
		return False

	def hgvscIsMultiple( self , hgvsc ):
		match = mafvariant.multiple.match( hgvsc )
		if match:
			return True
		return False

	def HGVSc( self ):
#		print "biomine::variant::mafvariant::HGVSc"
		hgvsc = ""
		if self.transcriptCodon:
			hgvsc += str(self.transcriptCodon)
		elif self.gene:
			hgvsc += str(self.gene)
		hgvsc += ":c."
		if self.positionCodon:
			hgvsc += str(self.positionCodon).replace( "-" , "_" )
		if self.reference != "-" and self.alternate != "-":
			if self.start == self.stop: #snp
				hgvsc += str(self.reference) + ">" \
				+ str(self.alternate)
			else: #complex
				hgvsc += "del" + str(self.reference) \
				+ "ins" + str(self.alternate)
		elif self.reference == "-" and self.alternate != "-":
			hgvsc += "ins" + str(self.alternate)
		elif self.reference != "-" and self.alternate == "-":
			hgvsc += "del" + str(self.reference)
		return hgvsc
	def HGVSct( self ):
		return str(self.transcriptCodon) + ":c." \
		+ str( self.positionCodon ).replace( "-" , "_" ) \
		+ str(self.reference) + "_" \
		+ str(self.alternate)
	def codingHGVS( self ):
#		print "biomine::variant::mafvariant::codingHGVS"
		return self.HGVSc() + '::' + self.HGVSp()
	def proteogenomicVar( self ):
#		print "biomine::variant::mafvariant::proteogenomicVar"
		return self.genomicVar() + "::" + self.codingHGVS()
	def uniqueProteogenomicVar( self ):
#		print "biomine::variant::mafvariant::uniqueProteogenomicVar"
		if self.sample:
			return self.sample + "::" + self.proteogenomicVar()
		else:
			return "nosample::" + self.proteogenomicVar()

	def splitHGVSc( self , hgvsc , xDot="c\." , override = False ):
		#print "biomine::variant::mafvariant::splitHGVSc - " ,
		ref = ""
		pos = ""
		mut = ""
		if not hgvsc:
			return [ ref , pos , mut ]
		shgvsc = hgvsc.split( ":" )
		if len( shgvsc ) > 1:
			self.transcriptCodon = shgvsc[0]
			hgvsc = shgvsc[1]
		else:
			hgvsc = shgvsc[0]
		pattern = re.compile( ".*" + xDot + "(.*)" )
		change = pattern.match( hgvsc )
		#print( change )
		if change:
			hgvsc = change.group( 1 )
			#print( hgvsc )
			isInd = self.hgvscIsIndel( hgvsc ) 
			isNon = self.hgvscIsNonCoding( hgvsc ) 
			isDup = self.hgvscIsDuplication( hgvsc )
			isInv = self.hgvscIsInversion( hgvsc )
			isMul = self.hgvscIsMultiple( hgvsc )
			posOnly = self.hasCodonPositionOnly( hgvsc )
			#print( ' | '.join( [ str( hgvsc ) , str( isNon ) , str( isMul ) , str( posOnly ) , str( isInd ) , str( isDup ) , str( isInv ) ] ) )
			if not posOnly:
				if isInd != 0:
					if isInd == -1: #simple deletion
						return self.splitDeletionHGVSc( hgvsc , noncoding = isNon , multiple = isMul , override = override )
					elif isInd == 1: #simple insertion
						return self.splitInsertionHGVSc( hgvsc , noncoding = isNon , multiple = isMul , override = override )
					elif isInd == 2: #then its a complex indel
						return self.splitComplexHGVSc( hgvsc , noncoding = isNon , multiple = isMul , override = override )
				elif isDup:
					return self.splitDuplicationHGVSc( hgvsc , noncoding = isNon , multiple = isMul , override = override )
				elif isInv:
					return self.splitInversionHGVSc( hgvsc , noncoding = isNon , multiple = isMul , override = override )
				else: #then its a snv
					return self.splitSNVHGVSc( hgvsc , noncoding = isNon , override = override )
			else:
				if isMul:
					startStop = hgvsc.split('_')
					pos = startStop[0]
				else:
					pos = hgvsc
				self.positionCodon = pos
			return [ "" , pos , "" ]
		else:
			print "biomine::variant::mafvariant Warning: could not find HGVS codon change"
			print "  Hint: Is the input codon column correct?"
			print "    Problem variant: " ,
			self.printVariant(',')
		return [ "" , "" , "" ]

	def hasCodonPositionOnly( self , hgvsc ):
#		print "biomine::variant::mafvariant::hasCodonPositionOnly - " ,
		noncoding = self.hgvscIsNonCoding( hgvsc )
		if noncoding:
			return True
		else:
			pattern = re.compile( ".*\d+.*([a-zA-Z])" )
			pos = pattern.match( hgvsc )
			#print( hgvsc + "  " + str( pos ) )
			if pos and noncoding != "NULL":
				return False
			return True

	def splitSNVHGVSc( self , hgvsc , noncoding = False , override = False ):
		#print( "splitSNVHGVSc: " + hgvsc + ", " + str( noncoding ) + ", " + str( override ) )
		ref = ""
		pos = ""
		mut = ""
		pattern = None
		if noncoding:
			pattern = re.compile( "(.*)([a-zA-Z\*])>([a-zA-Z\*])" )
		else:
			pattern = re.compile( "(\d+?)([a-zA-Z\*])>([a-zA-Z\*])" )
		change = pattern.match( hgvsc )
		if change:
			ref = change.group( 2 )
			pos = change.group( 1 )
			mut = change.group( 3 )
			if override or self.reference == "" and self.alternate == "":
				self.reference = ref
				self.alternate = mut
			self.positionCodon = pos
		return [ ref , pos , mut ]

	def defaultNull( self , val , null = "." ):
		if val:
			return val
		else:
			return null

	def splitComplexHGVSc( self , hgvsc , multiple = False , override = False , null = "." ):
		#print "biomine::variant::mafvariant::splitComplexHGVSc - " ,
		matches = re.match( "(.*)del(.*)ins(.*)" , hgvsc )
		parts = matches.groups()
		#print( parts )
		ref = matches.group( 2 )
		pos = ""
		if multiple:
			positions = matches.group( 1 ).split( "_" )
			pos = positions[0]
		else:
			pos = matches.group( 1 )
		mut = matches.group( 3 )
		ref = self.defaultNull( ref , null = null )
		mut = self.defaultNull( mut , null = null )
		if override:
			self.reference = ref
			self.alternate = mut
		self.positionCodon = pos
		#print( ', '.join( [ str( ref ) , str( pos ) , str( mut ) ] ) )
		return [ ref , pos , mut ]

	def splitDeletionHGVSc( self , hgvsc , multiple = False , override = False , null = "." ):
		#print "biomine::variant::mafvariant::splitDeletionHGVSc - " ,
		matches = re.match( "(.*)del(.*)" , hgvsc )
		parts = matches.groups()
		#print( parts )
		ref = matches.group( 2 )
		pos = ""
		if multiple:
			positions = matches.group( 1 ).split( "_" )
			pos = positions[0]
		else:
			pos = matches.group( 1 )
		ref = self.defaultNull( ref , null = null )
		if override:
			self.reference = ref
			self.alternate = null
		self.positionCodon = pos
		#print( ', '.join( [ str( ref ) , str( pos ) , str( mut ) ] ) )
		return [ ref , pos , null ]

	def splitInsertionHGVSc( self , hgvsc , multiple = False , override = False , null = "." ):
		#print "biomine::variant::mafvariant::splitInsertionHGVSc - "
		matches = re.match( "(.*)ins(.*)" , hgvsc )
		parts = matches.groups()
		#print( parts )
		mut = matches.group( 2 )
		pos = ""
		if multiple:
			positions = matches.group( 1 ).split( "_" )
			pos = positions[0]
		else:
			pos = matches.group( 1 )
		mut = self.defaultNull( mut , null = null )
		if override:
			self.reference = null
			self.alternate = mut
		self.positionCodon = pos
		#print( ', '.join( [ str( ref ) , str( pos ) , str( mut ) ] ) )
		return [ null , pos , mut ]

	def splitInversionHGVSc( self , hgvsc , noncoding = False , multiple = False , override = False ):
		#print "biomine::variant::mafvariant::splitInversionHGVSc - "
		matches = re.match( "(.*)inv(.*)" , hgvsc )
		parts = matches.groups()
		#print( parts )
		mut = matches.group( 2 )
		pos = ""
		if multiple:
			positions = matches.group( 1 ).split( "_" )
			pos = positions[0]
		else:
			pos = matches.group( 1 )
		mut = self.defaultNull( mut , null = null )
		if override:
			self.reference = null
			self.alternate = mut
		self.positionCodon = pos
		#print( ', '.join( [ str( ref ) , str( pos ) , str( mut ) ] ) )
		if multiple:
			if noncoding:
				pattern = re.compile( "([\*|\-|\+]*\d+?[\-|\+]*[\d+?])([a-zA-Z\*]+)inv([a-zA-Z\*])" )
			else:
				pattern = re.compile( "(\d+?)([a-zA-Z\*]+)inv([a-zA-Z\*])" )
		else:
			if noncoding:
				pattern = re.compile( "([\*|\-|\+]*\d+?[\-|\+]*[\d+?])([a-zA-Z\*]+)inv([a-zA-Z\*])" )
			else:
				pattern = re.compile( "(\d+?)([a-zA-Z\*]+)inv([a-zA-Z\*])" )
		pattern = re.compile( "(\d+?)[\_\d.?]inv(\w+)" )
		matches = pattern.match( hgvsc )
		parts = matches.groups()
		mut = parts[1]
		pos = parts[0]
		self.reference = parts[1].reverse()
		self.positionCodon = parts[0]
		return [ pos , mut ]

	def splitDuplicationHGVSc( self , hgvsc , noncoding = False , multiple = False , override = False ):
		if multiple:
			if noncoding:
				pattern = re.compile( "([\*|\-|\+]*\d+?[\-|\+]*[\d+?])([a-zA-Z\*]+)dup([a-zA-Z\*])" )
			else:
				pattern = re.compile( "(\d+?)([a-zA-Z\*]+)dup([a-zA-Z\*])" )
		else:
			if noncoding:
				pattern = re.compile( "([\*|\-|\+]*\d+?[\-|\+]*[\d+?])([a-zA-Z\*]+)dup([a-zA-Z\*])" )
			else:
				pattern = re.compile( "(\d+?)([a-zA-Z\*]+)dup([a-zA-Z\*])" )
		pattern = re.compile( "([\*|-]\d+?[-\d.?|+\d.?).*dup(\w+)" )
		matches = pattern.match( hgvsc )
		parts = matches.groups()
		mut = parts[2]
		self.reference = parts[1]
		pos = parts[0]
		self.positionCodon = parts[0]
		return [ pos , mut ]

	def splitSimpleIndelHGVSc( self , hgvsc , noncoding = False , multiple = False , override = False ):
#		print "biomine::variant::mafvariant::splitSimpleIndelHGVSc - " ,
		pattern = re.compile( "(del)" )
		matches = pattern.match( hgvsc )
		if matches:
			self.splitDeletion( hgvsc )
		else:
			pattern = re.compile( "(dup)" )
			matches = pattern.match( hgvsc )
			if matches:
				self.splitDuplicationHGVSc( hgvsc )
			else:
				pattern = re.compile( "(inv)" )
				matches = pattern.match( hgvsc )
				if matches:
					self.splitInversionHGVSc( hgvsc )
				else:
					self.splitInsertion( hgvsc )

	def plausibleCodonFrame( self , otherVar ):
		start = int(self.start)
		#stop = self.stop
		ostart = int(otherVar.start)
		#ostop = otherVar.stop
		maxShift = 2
		if start+maxShift >= ostart and start-maxShift <= ostart:
			return True
		return False
'''
mu = variant(gene="BRAF",chromosome=7,start=12345,stop=123456,reference="AT",alternate="GC",referencePeptide="A123R")
mu.printVariant('\t')
nu = variant(gene="BRAF",referencePeptide="A123R")
nu.printVariant('\t')
ou = variant(gene="BRAF",chromosome=7,start=12345,stop=123456,reference="AT",alternate="GC")
ou.printVariant('\t')
pu = variant(chromosome=7,start=12345,stop=123456,reference="AT",alternate="GC")
pu.printVariant('\t')
qu = variant(chromosome=7,start=12345,reference="AT",alternate="GC")
qu.printVariant('\t')

print qu.attr()
'''
