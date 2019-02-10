from Bio import SeqIO
import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

fasta_path = '/Users/tanayajadhav/Desktop/CompGen/Tanaya_Jadhav_week13/ch_ests_6.fasta'
rows = []
hum_genes = []
for fasta in SeqIO.parse(fasta_path, 'fasta'):
    sequence = str(fasta.seq)
    print(fasta.id)
    blast_result_handle = NCBIWWW.qblast(program="blastx", database="refseq_protein", sequence=sequence,
                                         entrez_query="txid9606[ORGN]")

    blast_records = NCBIXML.parse(blast_result_handle)
    blast_record = next(blast_records)
    if len(blast_record.alignments) != 0:
        # print("# alignments:", len(blast_record.alignments))
        tophit = blast_record.alignments[0]
        goat_id = str(fasta.id)
        human_ortholog_id = str(tophit.hit_id)
        human_gene = str(tophit.hit_def)
        columns = {'goat_id': goat_id, 'human_ortholog_id': human_ortholog_id, 'human_gene': human_gene}
        print(columns)
        column_series = pd.Series(columns)
        rows.append(column_series)
        hum_genes.append(human_ortholog_id)
    else:
        pass
blast_result_handle.close()


## print file with list of human orthologs
with open('./genelist.txt', 'w') as o:
    for gene in hum_genes:
        gene = gene.split('|')[1]
        o.write(gene + '\n')

allrows = pd.concat(rows, axis=1).T
allrows = allrows.reindex(columns=['goat_id', 'human_ortholog_id', 'human_gene'])
table_filepath = './week13table.tsv'
allrows.to_csv(table_filepath, sep='\t', index=False, header=False)


