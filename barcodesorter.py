#!/usr/bin/env python3
# Author:Tanaya Jadhav
# Uses a list of known barcodes to find them in a fastq file and
# create 10 fastq files
# all sequences in 1 fastq file contain the same barcode

from itertools import islice


def main():
    barcodes = ['ATGAGATCTT', 'AGCTCATTTC', 'TGAAAATCTT', 'TATCCAGCCA', 'AGGCAGGCAG', 'CTTGTTACTA', 'AAGGCACAAG',
                'TGCTCGCTGA', 'GTACCGCCGT', 'CCTCACCAGC']
    filelist = []
    for i in barcodes:
        file_name_string = i + '_seq.fq'
        filelist.append(file_name_string)
    filedict = dict(zip(barcodes, filelist))
    # print(filedict)
    with open('/Users/tanayajadhav/Desktop/CompGen/week3_seqs_bc.fq', 'r') as infile:
        x = 5
        while x == 5:
            lines = []
            x=6
            for line in infile:
                x=5
                lines.append(line)
                if len(lines) >= 4:
                    bc = lines[1][:10]
                    if bc in filedict:
                        with open(filedict[bc], 'a') as o:
                            o.write(lines[0] + lines[1] + lines[2] + lines[3])
                    break

        print('done')


if __name__ == '__main__':
    main()
