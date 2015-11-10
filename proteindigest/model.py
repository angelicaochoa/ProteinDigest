# -*- coding: utf-8 -*-
"""This module contains the data model of the application."""


import pkg_resources
pkg_resources.require("SQLObject>=0.10.1")

from turbogears.database import PackageHub
# import some basic SQLObject classes for declaring the data model
# (see http://www.sqlobject.org/SQLObject.html#declaring-the-class)
from sqlobject import SQLObject, SQLObjectNotFound, RelatedJoin
from sqlobject.inheritance import InheritableSQLObject
# import some datatypes for table columns from SQLObject
# (see http://www.sqlobject.org/SQLObject.html#column-types for more)
from sqlobject import StringCol, UnicodeCol, IntCol, DateTimeCol

__connection__ = hub = PackageHub('proteindigest')


class EnzymeDigest:
    def __init__(self,name='',pepseq='',enzyme=''):
        self.name = name
        self.pepseq = pepseq
        self.enzyme = enzyme
        
    def read(self,filename):
        #Import the protein sequence
        self.pepseq = ''.join(open(filename).read().split())

    def setEnzyme(self,enzyme_name):
        #Set the enzyme attribute according to the selected enzyme name
        if enzyme_name == 'Trypsin':
            self.enzyme = 'Trypsin'
        elif enzyme_name == 'Proteinase K':
            self.enzyme = 'ProteinaseK'
        elif enzyme_name == 'Pepsin (pH=1.3)':
            self.enzyme = 'Pepsin13'
        elif enzyme_name == 'Pepsin (pH>2.0)':
            self.enzyme = 'Pepsin2'
    
    def peptidedigest(self,cleavage):
        #Digest the peptide according to which enzyme has been selected
        if self.enzyme == 'Trypsin':
            subMasses = self.Trypsin()
            missedCleaves = self.MCs(subMasses,cleavage)
        elif self.enzyme == 'ProteinaseK':
            subMasses = self.ProteinaseK()
            missedCleaves = self.MCs(subMasses,cleavage)
        elif self.enzyme == 'Pepsin13':
            subMasses = self.Pepsin13()
            missedCleaves = self.MCs(subMasses,cleavage)
        elif self.enzyme == 'Pepsin2':
            subMasses = self.Pepsin2()
            missedCleaves = self.MCs(subMasses,cleavage)                                   
        return missedCleaves
    
##    def MCs(self,subMasses,cleavage):
    def MCs(self,cleavageSites,cleavage):
        #Generate the peptide list with missed cleavages
        peptideList = {}
        sequence = ''.join(self.pepseq[:].strip())

        for mc in range(cleavage+1):
            masses = []
            for i in range(len(cleavageSites)-mc-1):
                i2 = i+1+mc
                ind1 = cleavageSites[i]
                ind2 = cleavageSites[i2]

                mass = sequence[ind1:ind2]
                masses.append([ind1,mass])
            peptideList[mc] = masses

        return peptideList
    
        
    def Trypsin(self):
        #Digest at C-terminal of K or R
        #Exceptions: if P is C-terminal to K or R
        peptide = self.pepseq

        siteList = 'KR'
        cleaveSites_ = [i+1 for i,p in enumerate(peptide) if p in siteList]

        #Remove any cleavage sites that do not agree with exception
        cleaveSites = [0]
        for i in range(len(cleaveSites_)):
            site = cleaveSites_[i]
            if site+1 < len(peptide):
                if peptide[site+1] != 'P':
                    cleaveSites.append(site)
        cleaveSites.append(len(peptide))

        return cleaveSites

    def ProteinaseK(self):
        #Digest at C-terminal of A, F, Y, W, L, I, V
        #Exceptions: N/A
        peptide = self.pepseq
        
        siteList = 'AFYWLIV'
        cleaveSites = [i+1 for i,p in enumerate(peptide) if p in siteList]
        cleaveSites.extend([0,len(peptide)])
        cleaveSites = sorted(cleaveSites)

        return cleaveSites

    def Pepsin13(self):
        #Digest at C-terminal of F, L
        #Exceptions: N/A
        peptide = self.pepseq
        
        siteList = 'FL'
        cleaveSites = [i+1 for i,p in enumerate(peptide) if p in siteList]
        cleaveSites.extend([0,len(peptide)])
        cleaveSites = sorted(cleaveSites)

        return cleaveSites   

    def Pepsin2(self):
        #Digest at C-terminal of F, L, W, Y, A, E, Q
        #Exceptions: N/A
        peptide = self.pepseq
        
        siteList = 'FLWYAEQ'
        cleaveSites = [i+1 for i,p in enumerate(peptide) if p in siteList]
        cleaveSites.extend([0,len(peptide)])
        cleaveSites = sorted(cleaveSites)

        return cleaveSites  



# functions for populating the database

def bootstrap_model(clean=False):
    """Create all database tables and fill them with default data.

    This function is run by the 'bootstrap' function from the command module.
    By default it creates all database tables for your model.

    You can add more functions as you like to add more boostrap data to the
    database or enhance the function below.

    If 'clean' is True, all tables defined by your model will be dropped before
    creating them again.

    """
    create_tables(clean)

def create_tables(drop_all=False):
    """Create all tables defined in the model in the database.

    Optionally drop existing tables before creating them.

    """
    from turbogears.util import get_model
    from inspect import isclass

    model = get_model()
    if not model:
        from proteindigest.command import ConfigurationError
        raise ConfigurationError(
            "Unable to create database tables without a model")

    try:
        so_classes = [model.__dict__[x] for x in model.soClasses]
    except AttributeError:
        so_classes = model.__dict__.values()

    if drop_all:
        print "Dropping all database tables defined in model."
        for item in reversed(so_classes):
            if isclass(item) and issubclass(item, SQLObject) and \
                    item is not SQLObject and item is not InheritableSQLObject:
                item.dropTable(ifExists=True, cascade=True)

    # list of constraints we will collect
    constraints = list()

    for item in so_classes:
        if isclass(item) and issubclass(item, SQLObject) and \
                item is not SQLObject and item is not InheritableSQLObject:
            # create table without applying constraints, collect
            # all the constaints for later creation.
            # see http://sqlobject.org/FAQ.html#mutually-referencing-tables
            # for more info
            collected_constraints = item.createTable(ifNotExists=True,
                applyConstraints=False)

            if collected_constraints:
                constraints.extend(collected_constraints)

    # now that all tables are created, add the constaints we collected
    for postponed_constraint in constraints:
        # item is the last processed item and we borrow its connection
        item._connection.query(postponed_constraint)

    print "All database tables defined in model created."

