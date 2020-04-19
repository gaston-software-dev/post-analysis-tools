#!/usr/bin/python
# -*- coding: utf8 -*-

from .imports.readontology import output_str

class Error(Exception):
	"""Base class for exceptions in this module."""
	pass

class InputError(Error):
	"""Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

	def __init__(self, expression, message):
		self.expression = expression
		self.message = message
	
	@output_str	
	def __str__(self):
		st = "\nUnexpected value has occurred: Check %s\n"%(self.expression,)
		st += "Explanation: %s"%(self.message, )
		return st

