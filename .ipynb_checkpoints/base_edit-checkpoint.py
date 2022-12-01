# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:51:30 2019
This programm shall find all spacer designs to insert intendet mutation into your genome
@author: chdavo
"""

import Bio
from Bio import GenBank
from Bio import SeqIO
from Bio.Seq import Seq

pam="GG"

win_of_edi=6
dis_from_pam=18
end_of_win=dis_from_pam-win_of_edi
smallpam=3-len(pam)
codons_fw=["CAA","CAG","CGA"]
codons_rv=["TGG"]
genome=SeqIO.read("KT2440genome.fasta.txt", "fasta").seq #genome sequence to search against
cutoff=1 # proportion of coding sequence which is searched 1=100%; 0.5 only first half of gene
pam_rv_c=Seq(pam).reverse_complement()

count_targetedgenes=0
count_totalgenes=0
placestopcodon=[]

def Find_codon_position(full, codon):
    sub_index = 0
    position = -1
    for ch_i,ch_f in enumerate(full) :
        if ch_f.lower() != codon[sub_index].lower():
            position = -1
            sub_index = 0
        if ch_f.lower() == codon[sub_index].lower():
            if sub_index == 0 :
                position = ch_i

            if (len(codon) - 1) <= sub_index :
                break
            else:
                sub_index += 1

    return position

def Spacer_search():
     for n in range(int(len(seq_record.seq)*cutoff)):
        
        #forward strand
        if pam in seq_record[n+smallpam:3+n]: # search forward strand
            for codon in codons_fw:
                if codon in seq_record.seq[n-dis_from_pam:n-end_of_win+2]: # is editable codon in right distance to PAM ; +2 because the edited C is the first base of the codon, so they other bases can lay outside
                        codon_position=Find_codon_position(seq_record.seq[n-dis_from_pam:n-end_of_win],codon) #gives start index of search codon
                        if (codon_position+n-dis_from_pam)%3==0: #codon is in reading frame
                            if seq_record.seq[codon_position+n-dis_from_pam-1]!="G": #no GC
                                counts=genome.count(seq_record.seq[n-13:n])+genome.count(seq_record.seq[n-13:n].reverse_complement()) # spacer present in genome
                               # print(counts)
                              #  print(seq_record.seq[n-13:n])
                                if counts<2:
                              #      print ((codon_position+n-dis_from_pam)/len(seq_record.seq))
                                    placestopcodon.append(round(codon_position+n-dis_from_pam,2))
                    #                print (seq_record.seq[n-20:n])
                    #                print (seq_record.description)
                                    return(1)
                              
        
        if pam_rv_c in seq_record[n-3:n-smallpam]:
            for codon in codons_rv:
                if codon in seq_record.seq[n+end_of_win:n+dis_from_pam+1] or codon in seq_record.seq[n+dis_from_pam:n+dis_from_pam+2] and seq_record.seq[n+dis_from_pam-1]!="C":# +1 for the middle C and +2 for the last C, but than no G infront
                    codon_position=Find_codon_position(seq_record.seq[n+end_of_win:n+dis_from_pam],codon)
                    if (codon_position+n+end_of_win)%3==0:
                        counts=genome.count(seq_record.seq[n:n+13])+genome.count(seq_record.seq[n:n+13].reverse_complement())
                        if counts<2:
              #              print ((codon_position+n-dis_from_pam)/len(seq_record.seq))
                            placestopcodon.append(round((codon_position+n-dis_from_pam),2))
                         #   print (seq_record.seq[n:n+20].reverse_complement())
                  #          print (seq_record.description)
                            return(1)
     return(0)

for seq_record in SeqIO.parse("KT2440codingsequence.fasta", "fasta"):
    count_targetedgenes+=Spacer_search()
    count_totalgenes+=1
print(count_targetedgenes, count_totalgenes)
#print(placestopcodon)
with open ("results_PPngg.txt", "w+") as file:
     for item in placestopcodon:
        file.write("%s;" % item)
with open ("results_percentsuccessPPNGG.txt", "w+") as file:
    file.write(str(count_targetedgenes/count_totalgenes))
    
    
   
                        
          