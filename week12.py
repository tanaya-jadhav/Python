import pandas as pd
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Entrez
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

biomarthumanfile = '/Users/tanayajadhav/Desktop/CompGen/Tanaya_Jadhav_week12/biomarthuman.tsv'
bmrthum = pd.read_csv(biomarthumanfile, sep='\t', header=None).fillna('.')
bmrthum = bmrthum[1:]
filtered_bmrthum = bmrthum[bmrthum[1] != '.']

Entrez.email = "tanayaj@temple.edu"
rows = []
for index, ids in filtered_bmrthum.iterrows():
    # print(ids[0], ids[1])
    ensembl_id = ids[0]
    mrna_id = ids[1]
    resulthandle = Entrez.efetch(db="nucleotide", id=mrna_id, rettype="gb", retmode="text")
    genrecord = SeqIO.read(resulthandle, "genbank")
    seqrecord_id = genrecord.id
    columns = {'ensembl_id': ensembl_id, 'mrna_id': mrna_id, 'seqrecord_id': seqrecord_id}
    column_series = pd.Series(columns)
    rows.append(column_series)
    print(len(rows))
allrows = pd.concat(rows, axis=1).T
allrows = allrows.reindex(columns=['ensembl_id', 'mrna_id', 'seqrecord_id'])
bmrthum_filepath = './week12human.tsv'
allrows.to_csv(bmrthum_filepath, sep='\t', index=False, header=False)

# for id in idlist:
#     mrna = str(id)
#     resulthandle = Entrez.efetch(db="nucleotide", id=mrna, rettype="gb", retmode="text")
#     # print(resulthandle.readline().strip())
#     ereaddic = Entrez.read(resulthandle)
#
#     # resulthandle.close()