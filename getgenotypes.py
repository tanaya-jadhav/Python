
import math

def phredprob(Qchar):
    """
        calculate the probability of an error using the character code from the read data
        assumes that character coding is as follows:
            phredscore = (ASCII value of the character) minus 33
            use the ord() function to get ASCII value
    """
    phredscore = ord(Qchar)-33
    return pow(10,-float(phredscore)/10)



##for i in range(33,127):
##    print i,chr(i)

def getrefseq(fname):
    """
        return a string containing the first sequence in a fasta file
    """
    lines = open(fname).readlines()
    s = ""
    for l in lines[1:]:
        if l=='' or l[0] == '\n' or l[0]=='>':
            break
        else:
            s += l.strip()
    return s

def getreadinfo(fname):
    """
        pull info on the base position, sequence and quality from a set of reads taken form a bam file
        excludes those reads that have a CIGAR score that is not equal to the length of the sequence
        i.e. avoid reads spanning indels
    """
    f = open(fname)
    firstbases=[]
    seqs = []
    quals = []
    for line in f:
        # print(line)
        v = line.split()
        # print(v)
        seq = v[9]
        seqlen = len(seq)
        noindelCIGARstr = str(seqlen) + "M"
        if v[5] == noindelCIGARstr:  ## check to see if CIGAR string indicates no indels
            firstbases.append(int(v[3]))
            seqs.append(seq)
            quals.append(v[10])
    return firstbases,seqs,quals


def matchreads(refseq,refbase1num,firstbases,seqs,quals):
    """
        make a list
            one item for each base in refseq in order
            each item is a list and contains, in order
                the base number
                the reference base
                a list of bases that occurred in reads
                a corresponding list of quality values that occured in reads
    """
    ## by python numbering the first base in refseq is at position 0
    ## need to renumber of firstbases[] values, so the base positions line up
    r = []
    numbases = len(refseq)
    # print(numbases)
    for i in range(numbases):
        r.append([i,refseq[i],[],[]])
    numreads = len(firstbases)
    for j in range(numreads):
        k = firstbases[j]
        for ci,c in enumerate(seqs[j]):
            # print(ci, c)
            renum1 = (k+ci) - refbase1num
            # print(renum1)
            if 0 <= renum1 < numbases:
                # print(r[renum1][2])
                r[renum1][2].append(c)
                r[renum1][3].append(quals[j][ci])
    return r


refname = "hs37d5_22.fa"
samsampname = "NA19247.mapped.chr22.sorted.sam"
output_filename = 'genotypecallspartb.txt'
with open(output_filename, 'w') as o:
    o.write('#Position' + '\t' + 'likelihood_homozygous_alternate' + '\t'
            + 'likelihood_heterozygous' + '\t' + 'likelihood_homozygous_reference'
            + '\t' + 'genotype call' + '\n')
refbase1num = 1
refseq = getrefseq(refname)
print(refseq)
firstbases,seqs,quals = getreadinfo(samsampname)
# print(firstbases, seqs, quals)
readinfo = matchreads(refseq,refbase1num,firstbases,seqs,quals)
print(readinfo)
varrr = []
for rr in readinfo:
    isvariable = rr[2].count(rr[1]) != len(rr[2])
    if isvariable:
        varrr.append(rr)
# print(varrr)
##    print refbase1num+rr[0],rr[1],rr[2],rr[3],isvariable
# print ("positions with variable reads:" )
for rr in varrr:
    errorprobs = []
    for e in rr[3]:
        errorprobs.append(phredprob(e))
    # print (refbase1num+rr[0],rr[1],rr[2],rr[3],errorprobs)


    matchreferrors = []
    alterrors = []
    i = 0
    for baseval in rr[2]:
        if baseval == rr[1]:
            matchreferrors.append(errorprobs[i])
        else:
            if baseval != 'N':
                alterrors.append(errorprobs[i])
        i = i+ 1

    prob_list = []
    line_list = [rr[0]]
    for g in [0, 1, 2]:
        x_tot = 1
        for r in matchreferrors:
            x = ((2 - g) * r) + (g * (1 - r))
            x_tot = x_tot * x
        # print('x_tot', x_tot)
        y_tot = 1
        for a in alterrors:
            y = ((2 - g) * (1 - a) + (g * a))
            y_tot = y_tot * y
        # print(y_tot)
        if y_tot == 1 and x_tot == 1:
            p = (0.5 ** (len(rr[2]))) * x_tot * y_tot
            # print('if')
        elif y_tot == 1:
            p = (0.5 ** (len(rr[2]))) * x_tot
            # print('elif1')
        elif x_tot == 1:
            p = (0.5 ** (len(rr[2]))) * y_tot
            # print('elif')
        else:
            p = (0.5 ** (len(rr[2]))) * x_tot * y_tot
            # print('else')
        prob_list.append(p)
    line_list.append(prob_list)
    if prob_list[0] > prob_list[1] and prob_list[0] > prob_list[2]:
        genotype = 'homozygous alternate'
        # print('genotype: homozygous alternate')
    elif prob_list[1] > prob_list[0] and prob_list[1] > prob_list[2]:
        genotype = 'heterozygous'
        # print('genotype: heterozygous')
    elif prob_list[2] > prob_list[0] and prob_list[2] > prob_list[1]:
        genotype = 'homozygous reference'
        # print('genotype: homozygous reference')
    else:
        genotype = 'uncalled base'
        # print('uncalled base')
    line_list.append(genotype)
    # print(line_list)
    with open('genotypecallspartb.txt', 'a') as o:
        o.write(str(line_list[0]) + '\t'
                + str(line_list[1][0]) + '\t' + str(line_list[1][1]) + '\t' + str(line_list[1][2])
                + '\t' + line_list[2] + '\n')
    # print (rr, matchreferrors, alterrors  )
    # print(varrr)







