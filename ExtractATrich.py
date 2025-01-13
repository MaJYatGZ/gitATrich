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
    ATrep = ["no",0,""]
    mymatch = re.findall(ATrepeatPattern,myseq)
    if mymatch != []:
        ATrep[0] = "yes"
        matches = [match.group(0) for match in re.finditer(ATrepeatPattern,myseq)]
        k = ""
        for i in matches:
            k = k + str(i)
        ATrep[1] = len(k)
        #print(k)
    return ATrep

FindATrepeat(tmp)


def ExtractATrichSeq(fastaFile, outputFile, ATpairnumber, ATpairnumber2):
    myATrich = {}
    myATrich_lens = {}  #ATrich seqs' lengths on each chromosome.
    myChrom_lens = {}  #chrom or contig? length?
    Fout = open(outputFile+".genome.ATrich.txt","w")
    FoutBed = open(outputFile+".genome.ATrich.bed","w")    
    myIDs = []
    print(fastaFile)
    with gzip.open(fastaFile, "rt") as handle:
        #print("OK")
        for seq_record in SeqIO.parse(handle,"fasta"):
            myID = seq_record.id
            mySeq = str(seq_record.seq)
            mySeq = mySeq.upper()
            myLen = len(seq_record)
            if myLen > 1000000:
                myIDs.append(myID)
                myChrom_lens[myID] = myLen
                myATrich_lens[myID] = []
                myStart = 100
                myEnd = myStart + 100
                myTile = 100
                while myEnd < myLen-1000:
                    frag = mySeq[myStart:myEnd]
                    frag = frag.upper()
                    countA = frag.count("A")
                    countT = frag.count("T")
                    countN = frag.count("N")
                    freq1 = 100*(countA+countT)/100.0
                    if countN < 10 and freq1 >= 90.0:
                        newfrag = mySeq[myStart-100:myEnd+1000]
                        newfrag = newfrag.upper()
                        ATstart = 0
                        ATend = 0
                        ATrateStart = []
                        ATrateEnd = []                
                        for i in range(0,100):
                            ATrate = (newfrag[i:].count("A")+newfrag[i:].count("T"))/len(newfrag[i:])
                            ATrateStart.append(ATrate)
                        max_idx = ATrateStart.index(max(ATrateStart))
                        ATstart = myStart - 100 + max_idx                
                        for i in range(1,1000):
                            ATrate = (newfrag[:-i].count("A")+newfrag[:-i].count("T"))/len(newfrag[:-i])
                            ATrateEnd.append(ATrate)
                        max_idx = ATrateEnd.index(max(ATrateEnd))
                        ATend = myEnd + 1000 - max_idx
                        #########################
                        #Biopython chrom position
                        myATseq = mySeq[(ATstart):(ATend-1)]
                        #######
                        ATrep = FindATrepeat(str(myATseq))
                        #######
                        AT_number = int(ATrep[1])
                        if ATrep[0] == "yes" and AT_number > 2*ATpairnumber:
                            if AT_number <= 2*ATpairnumber2:
                                myATrich_lens[myID].append(len(myATseq))
                                #####################
                                #UCSC chrom position
                                line2 = myID+"\t"+str(int(ATstart+1))+"\t"+str(int(ATend-1))+"\t"+"."+"\t"+""+"\t"+"."+"\n"
                                line = myID+":"+str(int(ATstart+1))+"-"+str(int(ATend-1))+"\t"+myATseq.upper()+"\n"
                                Fout.write(line)
                                FoutBed.write(line2)                        
                                ATseqID = myID+":"+str(int(ATstart+1))+"-"+str(int(ATend-1))
                                if myID in myATrich.keys():
                                    myATrich[myID].append([ATseqID,int(ATstart+1),int(ATend-1),myATseq.upper()])
                                else:
                                    myATrich[myID] = []
                                    myATrich[myID].append([ATseqID,int(ATstart+1),int(ATend-1),myATseq.upper()])
                        myStart = myEnd + 1000
                        myEnd = myStart + 100
                    else:
                        myStart = myStart + int(myTile/2)
                        myEnd = myEnd + int(myTile/2)
    Fout2 = open(outputFile+".ATrich.report.txt","w")
    for k in myIDs:
        chromLen = myChrom_lens[k]
        ATrichNumber = len(myATrich_lens[k])
        ATrichMean = np.mean(myATrich_lens[k])
        ATrichMedian = np.median(myATrich_lens[k])
        Fout2.write(k+"\t"+str(chromLen)+"\t"+str(ATrichNumber)+"\t"+str(ATrichMean)+"\t"+str(ATrichMedian)+"\n")
    Fout.close()
    Fout2.close()
    FoutBed.close()
    return myATrich


print("###################")
ExtractATrichSeq("./testfiles/anoGam3.fa.gz", "./testfiles/anoGam3", 25, 5000)  
print("###################")


