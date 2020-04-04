# coding: utf-8
"""
PySML-dev applications
======================

This module implements some common semantic 
similarity measure applications: Enrichment
analysis, Entity identification and cluste-
ring or classification.
"""

__all__ = ["conceptenrichment", "entityidentification"]

from .conceptenrichment import ConceptEnrichment
from .entityidentification import EntityIdentification
