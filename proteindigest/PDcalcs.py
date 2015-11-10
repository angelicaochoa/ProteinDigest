def getAAMasses():
    import os.path
    #Import the amino-acid masses
    __dir__=os.path.dirname(os.path.abspath(__file__))
    filepath = os.path.join(__dir__,'aaMasses.txt')
    filedata = open(filepath,'rU').readlines()
##    filedata = open('aaMasses.txt').readlines()
    aaMasses = {}
    for l in filedata:
        l = l.split()
        aaMasses[l[0]] = map(float,l[1:])
    return aaMasses    

def monoMW(peptide):
    #Calculate the monoisotopic mass of the peptide
    weight = 0
    aaMasses = getAAMasses()
    for aa in peptide:
        weight += aaMasses.get(aa,[0,0])[0]
    return weight

def aveMW(peptide):
    #Calculate the average mass of the peptide
    weight = 0
    aaMasses = getAAMasses()
    for aa in peptide:
        weight += aaMasses.get(aa,[0,0])[1]
    return weight




