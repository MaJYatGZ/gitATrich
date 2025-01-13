from Bio import SeqIO
import sys
import re
import numpy as np
import gzip
import random
import os


ATrepeatPattern = re.compile(r"(AT){3,}|(TA){3,}")

tmp = "ATATATAAAATATATATCCCCCATACTACCAAAAAAAAAAAAAAATATATATATATATATATATATCGACGTACGATCGATCGATCGAGATATATATATATATATATATATCGATAGCTAGCTAGCTATATATATATATATATATATATATATATATATATATAGCGCATAAT"

#mymatch = re.findall(ATrepeatPattern,tmp)

def FindATrepeat(myseq):
    print("The original seq is: "+myseq)
    ATrep = ["no",0,""]
    mymatch = re.findall(ATrepeatPattern,myseq)
    if mymatch != []:
        ATrep[0] = "yes"
        matches = [match.group(0) for match in re.finditer(ATrepeatPattern,myseq)]
        k = ""
        for i in matches:
            print("match: "+i+"; length: "+str(len(i)))
            k = k + str(i)
        ATrep[1] = len(k)
        #print(k)
    print("The length of total matched ATrich seqs is: "+str(ATrep[1]))
    return ATrep

ATnumber = FindATrepeat(tmp)

