import csv
import MySQLdb


db = MySQLdb.connect('localhost',user='tuf08739')
cur =  db.cursor()
sql = 'CREATE DATABASE tuf08739_week11'
cur.execute(sql)

database = 'USE tuf08739_week11;'

cur.execute(database)

maketablestr = """CREATE TABLE gene_sequences (id INT UNSIGNED AUTO_INCREMENT PRIMARY KEY UNIQUE,
                    seq_id VARCHAR(30) NOT NULL,sequence VARCHAR(30) NOT NULL,cds_sequence VARCHAR(30) NOT NULL);"""

cur.execute(maketablestr)
db.commit()

genefile = open('week11.tsv','r')
genereader =csv.reader(genefile,delimiter = '\t')

for row in genereader:
    id = int(row[0])
    seq_id = '\"' + row[1] +'\"'
    sequence = '\"' + row[2] +'\"'
    cds_sequence = '\"' + row[3] + '\"'
    data = (id,seq_id,sequence,cds_sequence)
    insertstatement = """INSERT INTO gene_sequences VALUES (%d,%s,%s,%s);"""
    cur.execute(insertstatement % data)
db.commit()
db.close()
genefile.close()

