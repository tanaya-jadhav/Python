#!/usr/bin/env python3
# Author:Tanaya Jadhav

# from subprocess import call
import random
import sys


def get_seq_name(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_name = line.split()[0][1:]
                return seq_name

def run_samtools(filepath, start_base, stop_base):
    seq_name = get_seq_name(filepath)
    x = filepath.split('/')[-1].split('.')[0]
    outfile = './' + x + '_' + str(start_base) + '-' + str(stop_base) + '.fasta'
    cmd = ['samtools faidx', filepath, "\"" + seq_name + '\":' +
           str(start_base) + '-' + str(stop_base), '>', outfile]
    # print(outfile)
    # print(seq_name)
    # with open(outfile, 'w') as f:
    #     call(cmd, stdout=f)

##Takes a string and finds the reverse complement. Has been modified from previous version of this function which
##only found the complement but not the reverse complement. This now finds reverse complement.
def rev_comp(string):
    basedic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y',
               'Y': 'R', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
               'D': 'H', 'H': 'D'}
    newstring = ''
    for c in string:
        if c in basedic:
            newstring += basedic[c]
        else:
            newstring += c

    newstring = newstring[::-1]

    return newstring


def randread(seqlength, readlength):
    start = random.randint(0, seqlength-readlength-1)
    return  start, start + readlength


def main(filename, outfile, readlength, coverage):
# def main():
    # file_path = '/Users/tanayajadhav/Desktop/CompGen/Tanaya_Jadhav_Week2/chr22_206583718-500000-600000.fasta'
    file_path = filename
    seq = ""
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                pass
            else:
                seq += line.strip('\n')

    interval = randread(100000, 100)

    #obtaining reverse complement of the sequence
    rev_seq = rev_comp(seq)

    # ##for testing calculations for number of reads
    # coverage = 10
    # readlength = 100

    # calculating number of reads needed
    read_num = (coverage * int(len(seq)))//readlength
    loop_num = read_num//2

    seq_name = get_seq_name(filename)
    #generating file with reads
    with open(outfile, 'w') as o:
        contg_count = 0
        for reads in range(loop_num):
            contg_count += 1
            rand_interval = randread(100000, 100)
            read = seq[rand_interval[0]:rand_interval[1]]
            o.write('@' + seq_name + ' contg' + str(contg_count) + '\n' + read + '\n' + '+' + '\n' + 'I'*100 + '\n')
            contg_count += 1
            rev_rand_interval = randread(100000, 100)
            revcomp_read = rev_seq[rev_rand_interval[0]:rev_rand_interval[1]]
            o.write('@' + seq_name + ' contg' + str(contg_count) +
                    '\n' + revcomp_read + '\n' + '+' + '\n' + 'I'*100 + '\n')










if __name__ == '__main__':
    filename = sys.argv[1]
    outfile = sys.argv[2]
    readlength = int(sys.argv[3])
    coverage = int(sys.argv[4])
    main(filename, outfile, readlength, coverage)
    # main()