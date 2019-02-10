import pandas as pd
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class Sequence:
    def __init__(self, handle):
        self.handle = handle

    def feature(self):
        for feature in self.handle.features:
            if feature.type == 'CDS':
                return feature.extract(self.handle.seq)

    def aminoacid(self):
        cds = self.feature()
        translatedseq = cds.translate()
        # print(translatedseq)
        protein = ProteinAnalysis(str(translatedseq))
        aminodic = protein.count_amino_acids()
        aminolist = list(aminodic.values())
        return aminodic, aminolist




Entrez.email = "hey@temple.edu" ##online requrests using Entrez  should include a valid email address
## search the nucleotide data base for sequences from Humans that have "type II diabetes"  somewhere in the text and have a CDS in their feature table
resulthandle = Entrez.esearch(db="nucleotide", retmax=10, term=""" ((("asthma")) AND "Homo sapiens"[Organism] AND CDS[Feature key]) """)
ereaddic = Entrez.read(resulthandle)  ## make a dictionary from the results.
resulthandle.close()

## get accession numbers for the Gene ID numbers in ereaddic["IdList"]
fetchhandle = Entrez.efetch(db="nucleotide", id=ereaddic["IdList"], rettype="acc")
accnums = fetchhandle.read().splitlines()
fetchhandle.close()


srlist = []
# srdic = {}
for accnum in accnums:
    try:
        fetchhandle = Entrez.efetch(db="nucleotide", id=accnum, rettype="gb", retmode="text")
        gbrecord = SeqIO.read(fetchhandle, "genbank")
        fetchhandle.close()
        # srdic[gbrecord.id] = gbrecord
        grs = Sequence(gbrecord)
        srlist.append(grs)
        # print(srdic)
        # print(' found')

    except:
        pass
        # print (" not found")



## write CDS sequences to a fasta file
cds_filepath = './codingsequences.fa'
with open(cds_filepath, 'w') as cds:
    for sr in srlist:
        idname = sr.handle.id
        seq = sr.feature()
        # print(idname)
        # print(seq)
        # print(sr.handle.seq)
        cds.write('>' + idname + '\n' + str(seq) + '\n')

## write amino acid counts to a file with a table
table_rows = []
for sr in srlist:
    aminodic, aminolist = sr.aminoacid()
    aminodic['id'] = sr.handle.id
    table_columns = pd.Series(aminodic)
    # print(table_columns)
    table_rows.append(table_columns)
# print(table_rows)
aminoacidfilepath = './aminoacids.csv'
aminoacids = pd.concat(table_rows, axis=1).T
aminoacids = aminoacids.reindex(columns=['id', 'A', 'C', 'D', 'E',
                                         'F', 'G', 'H', 'I', 'K',
                                         'L', 'M', 'N', 'P', 'Q',
                                         'R', 'S', 'T', 'V', 'W', 'Y'])
aminoacids.to_csv(aminoacidfilepath, sep='\t', index=False)



