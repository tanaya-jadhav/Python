#!/usr/bin/env python3
## Authors: Tanaya Jadhav

def rev_comp(line):
    basedic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'R': 'Y',
               'Y': 'R', 'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
               'D': 'H', 'H': 'D'}
    newline = ''
    for c in line:
        if c in basedic:
            newline += basedic[c]
        else:
            newline += c

    return newline

def totals(list, dic):
    sum = 0
    for i in list:
        sum += dic[i]
    return sum

def main():
    fasta_path = './GRCH38p2_chr22.fasta'
    with open(fasta_path, 'r') as f:
        rev = ''
        countdic = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'R': 0, 'Y': 0,
                    'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0,
                    'H': 0, 'V': 0, 'N': 0, '-': 0}
        rev_countdic = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'R': 0, 'Y': 0,
                    'S': 0, 'W': 0, 'K': 0, 'M': 0, 'B': 0, 'D': 0,
                    'H': 0, 'V': 0, 'N': 0, '-': 0}
        sig_count = 0
        for line in f:
            if line.startswith('>') or line.startswith('@') or line.startswith('+'):
                pass
            else:
                ## counts bases in file
                count_sig = line.count('ATGTTG')
                sig_count += count_sig

                count_A = line.count('A')
                countdic['A'] += count_A

                count_T = line.count('T')
                countdic['T'] += count_T

                count_C = line.count('C')
                countdic['C'] += count_C

                count_G = line.count('G')
                countdic['G'] += count_G

                count_R = line.count('R')
                countdic['R'] += count_R

                count_Y = line.count('Y')
                countdic['Y'] += count_Y

                count_S = line.count('S')
                countdic['S'] += count_S

                count_W = line.count('W')
                countdic['W'] += count_W

                count_K = line.count('K')
                countdic['K'] += count_K

                count_M = line.count('M')
                countdic['M'] += count_M

                count_B = line.count('B')
                countdic['B'] += count_B

                count_D = line.count('D')
                countdic['D'] += count_D

                count_H = line.count('H')
                countdic['H'] += count_H

                count_V = line.count('V')
                countdic['V'] += count_V

                count_N = line.count('N')
                countdic['N'] += count_N

                count_gap = line.count('-')
                countdic['-'] += count_gap

                rev += rev_comp(line)

        ## calculates totals and other metrics
        baselist = ['A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'D', 'H', 'V']
        tot_bases = totals(baselist, countdic)
        uncalledlist = ['N', '-']
        tot_uncalled = totals(uncalledlist, countdic)
        twoamblist = ['R', 'Y', 'S', 'W', 'K', 'M']
        tot_twoamb = totals(twoamblist, countdic)
        threeamblist = ['B', 'D', 'H', 'V']
        tot_threeamb = totals(threeamblist, countdic)
        GCcontent = ((countdic['G'] + countdic['C'])/tot_bases)*100


        ## counts bases in reverse complement
        for line in rev:
            count_A = line.count('A')
            rev_countdic['A'] += count_A

            count_T = line.count('T')
            rev_countdic['T'] += count_T

            count_C = line.count('C')
            rev_countdic['C'] += count_C

            count_G = line.count('G')
            rev_countdic['G'] += count_G

            count_R = line.count('R')
            rev_countdic['R'] += count_R

            count_Y = line.count('Y')
            rev_countdic['Y'] += count_Y

            count_S = line.count('S')
            rev_countdic['S'] += count_S

            count_W = line.count('W')
            rev_countdic['W'] += count_W

            count_K = line.count('K')
            rev_countdic['K'] += count_K

            count_M = line.count('M')
            rev_countdic['M'] += count_M

            count_B = line.count('B')
            rev_countdic['B'] += count_B

            count_D = line.count('D')
            rev_countdic['D'] += count_D

            count_H = line.count('H')
            rev_countdic['H'] += count_H

            count_V = line.count('V')
            rev_countdic['V'] += count_V

            count_N = line.count('N')
            rev_countdic['N'] += count_N

            count_gap = line.count('-')
            rev_countdic['-'] += count_gap


        ## Print statements
        print('The given sequence has the following base counts: ',
              'A:', countdic['A'],
              'T:', countdic['T'],
              'G:', countdic['G'],
              'C:', countdic['C'])
        print('Total number of called bases in sequence: ', tot_bases)
        print('Number of uncalled bases in sequence: ', tot_uncalled)
        print('Number of two base ambiguities: ', tot_twoamb)
        print('Number of three base ambiguities: ', tot_threeamb)
        print('GC content: ', GCcontent, '%')
        print("Number of occurrences of 'ATGTTG' in the sequence: ", count_sig)

        print('The reverse complement of the given sequence has the following base counts: ',
              'A:', rev_countdic['A'],
              'T:', rev_countdic['T'],
              'G:', rev_countdic['G'],
              'C:', rev_countdic['C'])





if __name__ == '__main__':
    main()