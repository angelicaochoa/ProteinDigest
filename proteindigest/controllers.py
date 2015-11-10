# -*- coding: utf-8 -*-
"""This module contains the controller classes of the application."""

# symbols which are imported by "from proteindigest.controllers import *"
__all__ = ['Root']


# third-party imports
from turbogears import controllers, expose
import Bio.SeqIO
import Bio.Alphabet
#=============================SET UP FORM FOR WELCOME PAGE===========================================
from turbogears import validate, validators
from turbogears import widgets, error_handler

class SearchFields(widgets.WidgetsList):
    peptide_seq = widgets.TextArea(label="")
    enzyme_sel = widgets.SingleSelectField(label="Select an enzyme: ",
                                           options=["Trypsin",
                                                    "Proteinase K",
                                                    "Pepsin (pH=1.3)",
                                                    "Pepsin (pH>2.0)"],
                                           default="Trypsin")
    missed_cleaves = widgets.SingleSelectField(label="# of Missed Cleavages: ",
                                               options=["0","1","2","3","4","5"],
                                                default="0")
    min_mass = widgets.SingleSelectField(label="Min Peptide Mass: ",
                                         options=["0","500","750","1000","1250","1500","1750"],
                                         default="0")
    max_mass = widgets.SingleSelectField(label="Max Peptide Mass: ",
                                         options=["3000","4000","5000","6000","7000","8000","unlimited"],
                                         default="unlimited")
    min_length = widgets.TextField(label="Min Peptide Length: ",default="3")
    max_length = widgets.TextField(label="Max Peptide Length: ",default="100")
    
class SearchFieldsSchema(validators.Schema):
    peptide_seq = validators.String(min=1,strip=True)
    enzyme_sel = validators.OneOf(["Trypsin","Proteinase K","Pepsin (pH=1.3)","Pepsin (pH>2.0)"])
    missed_cleaves = validators.OneOf(["0","1","2","3","4","5"])
    min_mass = validators.OneOf(["0","500","750","1000","1250","1500","1750"])
    max_mass = validators.OneOf(["3000","4000","5000","6000","7000","8000","unlimited"])
    min_length = validators.All(validators.PlainText(),
                                validators.String(min=1,strip=True))
    max_length = validators.All(validators.PlainText(),
                                validators.String(min=1,strip=True))  

search_form = widgets.TableForm(
    fields = SearchFields(),
    validator = SearchFieldsSchema(),
    action = "digest",
    submit_text = "Perform"
    )
#===================================================================================================


from model import EnzymeDigest
from PDcalcs import *


#=========================READ THE PROTEIN SEQUENCE OR IMPORT FROM UNIPROT==========================
def readSequence(peptide_seq):
    if peptide_seq.find('>') >= 0:
        #FASTA format
        fasta_data = peptide_seq.split('\n')
        name = fasta_data[0][1:]
        sequence = ''.join(fasta_data[1:])
        sequence = ''.join(sequence.strip())
        return name,sequence
    elif not peptide_seq.isalpha() and len(peptide_seq) == 6:
        #Import UniProt FASTA file
        import urllib
        fasta_url = 'http://www.uniprot.org/uniprot/'+peptide_seq+'.fasta'
        fasta_data = urllib.urlopen(fasta_url).read().split('\n')
        name = fasta_data[0][1:]
        sequence = ''.join(fasta_data[1:])
        sequence = ''.join(sequence.strip())

        return name,sequence
    else:   
        name = 'Unnamed_Protein'
        sequence = ''.join(peptide_seq[:].splitlines())
      
        return name,sequence


#=========================STORE RELEVANT PROTEIN DATA INTO PEP_DICT=================================
def getPeptideData(digested_frags,sequence):
    pep_dict = {}

    for mc,frag in digested_frags.items():            
        for f in frag:        
            peptide = f[1]  #peptide sequence
            peptide_mass = aveMW(peptide)   #peptide mass based on average amino-acid values
            start_pos = f[0]+1  #start position of the peptide sequence
            end_pos = f[0]+len(peptide) #end position of the peptide sequence

            #peptide position withni the protein sequence
            pep_position = str(start_pos)+'-'+str(end_pos)
            
            #Get amino-acids to the left and right of the peptide
            if start_pos == 1:
                left_aa = ''
            else:
                left_aa = sequence[start_pos-2]

            if end_pos+2 == len(sequence):
                right_aa = ''
            else:
                try:
                    right_aa = sequence[end_pos+1]
                except IndexError:
                    right_aa=''

            #combine the left and right amino-acids into a string
            LnR_aa = left_aa+' | '+right_aa

            #put all of the values together in a list for the dictonairy
            flist = [peptide_mass,pep_position,len(peptide),mc,LnR_aa,peptide]

            if pep_dict.get(peptide,[]) == []:
                pep_dict[peptide] = flist        
            else:
                mcentry = pep_dict.get(peptide)[2]
                if mc > mcentry:
                    pep_dict[peptide] = flist
            
    return pep_dict.values()
#===================================================================================================



#================================CONTROLLERS.PY FUNCTIONS===========================================    
class Root(controllers.RootController):
    """The root controller of the application."""

    @expose(template="proteindigest.templates.welcome")
    def index(self):
        """"Show the welcome page."""
        # log.debug("Happy TurboGears Controller Responding For Duty"))
        return dict(form=search_form)

    @expose(template="proteindigest.templates.digest")
    @expose(template="proteindigest.templates.peptideresultsxml", as_format="xml",format="xml")
    @validate(form=search_form)
    @error_handler(index)
    def digest(self,peptide_seq,enzyme_sel,missed_cleaves,min_mass,max_mass,min_length,max_length):
        #Digest the protein according to the parameters specified by the user
        
        pepInput = EnzymeDigest()
        pepInput.name,pepInput.pepseq = readSequence(peptide_seq)
        pepInput.setEnzyme(enzyme_sel)
        
        #generate list of digested peptide fragements including the missed cleavages
        digested_frags = pepInput.peptidedigest(int(missed_cleaves))
        
        #calculate the peptide data (i.e., mass, length, etc.)
        unfiltered_pep_fragments = getPeptideData(digested_frags,pepInput.pepseq)

        #Filter out peptides that do not meet mass or length criteria
        mass_filtered_pep_fragments = []
        for peptide in unfiltered_pep_fragments:
            mass = peptide[0]
            if max_mass == 'unlimited':
                if mass >= int(min_mass):
                    mass_filtered_pep_fragments.append(peptide)
            else:
                if mass >= int(min_mass) and mass <= int(max_mass):
                    mass_filtered_pep_fragments.append(peptide)
        length_filtered_pep_fragments = []
        for peptide in mass_filtered_pep_fragments:
            length = len(peptide[-1])
            if length >= int(min_length) and length <= int(max_length):
                length_filtered_pep_fragments.append(peptide)
            
        #sort the filtered peptide fragments by their mass
        pep_fragments = sorted(length_filtered_pep_fragments,key=lambda mass: -mass[0])

        #Total peptides returned after filtering..
        pct=round(100*len(pep_fragments)/float(len(unfiltered_pep_fragments)),2)
        percent_covered = [len(pep_fragments),len(unfiltered_pep_fragments),pct]

        #generate the url that will take the user to the XML view of the table
        urlstring = '/digest/'+pepInput.pepseq+'/'+enzyme_sel+'/'+missed_cleaves+'/'+min_mass+'/'+max_mass+'/'+min_length+'/'+max_length+'?tg_format=xml'
                
        return dict(name=pepInput.name,sequence=pepInput.pepseq,enzyme=enzyme_sel,mc=missed_cleaves,minmass=min_mass,maxmass=max_mass,minlength=min_length,maxlength=max_length,pep_fragments=pep_fragments,urlstring=urlstring,percent_covered=percent_covered)


            
            
    
