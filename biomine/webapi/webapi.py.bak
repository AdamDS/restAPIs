#!/usr/bin/python
# author: Adam D Scott (amviot@gmail.com)
# first created: 2015*09*28

#http://docs.python-requests.org/en/latest/user/authentication/?highlight=authentication
import requests
from requests.auth import HTTPDigestAuth
import xml.etree.ElementTree as ET
import json
import AdvancedHTMLParser
import time

class webapi(object):
	'''Web api class, has 
		endpoint = the endpoint
		subset = the subset realm of the Web service
		action = the query, file upload, etc.'''
	nullXML = "<None></None>"
	def __init__(self,endpoint,subset):
		self.response = None
		self.endpoint = endpoint
		self.subset = subset
		self.action = ""
		self.headers = {}
		self.data = {}

	def printInfo( self ):
		print "Response = " + str( self.response ) + ", ok? " + str( self.response.ok )
		print "Endpoint = " + self.endpoint
		print "Subset = " + self.subset
		print "Action = " + self.action
		print "URL = " + self.buildURL()
		print "Headers = " ,
		print self.headers
		print "Data = " ,
		print self.data

	def setSubset(self,subset):
		self.subset = subset
		self.action = ""
	def resetURL(self):
		self.action = ""
	def fullReset( self ):
		self.resetAction()
		self.resetData()
		self.resetHeaders()
	def resetAction( self ):
		self.action = ""
	def resetHeaders( self ):
		self.headers = {}
	def resetData( self ):
		self.data = {}

	def beginQuery(self):
		self.action = ""

	def queryAll(self):
		return self.site + "?query=NOT%20asdf"

	def partial(self,term,value):
		self.action += term + ":" + str(value)
	def addPartial(self,condition,term,value):
		self.action += "%20" + condition + "%20" + term + ":" + str(value)

	def exact(self,term,value):
		self.action += term + "=" + str(value)
	def addExact(self,condition,term,value):
		self.action += "%20" + condition + "%20" + term + "=" + str(value)

	def multiOr(self,term,values):
		values = map(str,values)
		stripped = [ v.strip() for v in values ]
		self.action += term + ":OR(" + ','.join(stripped) + ")"
	def addMultiOr(self,condition,term,values):
		values = map(str,values)
		stripped = [ v.strip() for v in values ]
		self.action += "%20" + condition + "%20" + term + ":OR(" + ','.join(stripped) + ")"

	def buildURL(self):
		return self.endpoint + self.subset + self.action

	def buildURLJSON(self):
		return self.endpoint + self.subset + self.action + "&_type=json"

	def addData( self , field , value ):
		self.data[field] = value
	def buildData( self ):
#		print "biomine::webapi::webapi::buildData"
#		print self.data
		data = ""
		if self.data:
			data = json.dumps( self.data )
			self.headers["Accept"] = "application/json"
#		print self.data
		return data

	def addHeader( self , field , value ):
		self.headers[field] = value
	def buildHeader( self ):
		return json.dumps( self.headers )

	def submit(self,**kwargs):
		contentHeaders = kwargs.get("content",'')
		dataIn = kwargs.get("data",'')
		doPost = kwargs.get("post",False)
		#asSession = kwargs.get( "session" , False )
		timeout = kwargs.get( "timeout" , None )
		#if asSession:
		#	request = requests.Session()
		if contentHeaders:
			self.addHeader( "Content-Type" , contentHeaders )
		url = self.buildURL()
		headers = self.headers #buildHeader()
		data = self.buildData()
#		print url
#		print headers
#		print data
		try:
			if self.headers:
				if self.data:
					if doPost:
						self.response = requests.post( url , headers = headers , data = data , timeout = timeout )
#					print self.response
					else:
						self.response = requests.get( url , headers = headers , data = data , timeout = timeout )
				else:
					if doPost:
						self.response = requests.post( url , headers = headers , timeout = timeout )
					else:
						self.response = requests.get( url , headers = headers , timeout = timeout )
			else:
				if self.data:
					self.buildData()
					if doPost:
						self.response = requests.post( url , data = data , timeout = timeout )
					else:
						self.response = requests.get( url , data = data , timeout = timeout )
				else:
					if doPost:
						self.response = requests.post( url , timeout = timeout )
					else:
						self.response = requests.get( url , timeout = timeout )
			self.errorCheck()
		except:
			print "biomine::webapi::submit failed: " ,
			self.errorCheck()
			pass
		return self.response
	
	def errorCheck( self ):
		if not self.response.status_code:
			print "biomine::webapi::errorCheck Warning: no status code when trying " + self.buildURL()
			return
		code = self.response.status_code
		if code != 200:
			if code == 204:
				print "biomine::webapi::errorCheck Warning: no content from " + self.buildURL()
			else:
				print "biomine::webapi::errorCheck Warning: response from " + self.buildURL()
				print "received status code = " + str(code)
				self.printInfo()
	def testURL( self , **kwargs ):
		skip = kwargs.get( 'skip' , True )
		if skip:
			print( "biomine::webapi::testURL Warning: site test skipped - " + self.url )
			return True
		print( "biomine::webapi::testURL behavior: testing url - " + self.url )
		response = self.submit()
		if response.status_code < 300:
			return True
		return False

	def submitDigest(self,username,password):
		url = self.buildURL()
		self.response = requests.get( url , auth=HTTPDigestAuth( username , password ) )
		return self.response

	def submitJSON(self):
		url = self.buildURL()
		self.submit( content = "application/json")

	def submitXML(self):
		url = self.buildURL()
		self.submit( content = "text/xml")
	
	def parseHTMLResponse( self ):
		page = AdvancedHTMLParser.AdvancedHTMLParser()
		page.parseStr( self.response.text )
		return page
	
	def getXMLRoot( self ):
		try:
			return ET.fromstring( self.response.text )
		except:
			print "biomine::webapi::getXMLRoot Warning: no root xml"
			return ET.fromstring( webapi.nullXML )
	def getEntry( self , generator , text ):
		try:
			for entrygen in generator.iter( text ):
				return entrygen.text
		except:
			print "biomine::webapi::getEntry Warning: no xml entry for " + text + " from " + generator.text
			return webapi.nullXML
	def getElement( self , generator , text ):
		try:
			for entrygen in generator.iter( text ):
				return entrygen
		except:
			print "biomine::webapi::getElement Warning: no xml element for " + text + " from " + generator.text
			return ET.fromstring( webapi.nullXML )

	@staticmethod
	def runTime( initialTime ):
		return time.time() - initialTime
	@staticmethod
	def printRunTime( step , interval ):
		print "Running " + step + " took " + str( interval ) + "seconds"
