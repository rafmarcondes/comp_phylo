# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 11:28:36 2015

@author: Rafael

This is code written to complete a sequence manipulation exercise for Jeremy Brown's CompPhylo class at LSU, Spring 2015.
Comments preceded by "INSTRUCTIONS" are Jeremy's literal instructions for the exercise; all other comments are mine.
"""


"""
INSTRUCTIONS:
- Create a new Python script (text file)
- At the beginning of the script, define a DNA sequence (taken from 
https://github.com/jembrown/CompPhylo_Spr2015/blob/master/CodingSeq.txt)
- Print the length of the sequence to the screen along with text explaining 
the value
"""
seq="aaaagctatcgggcccataccccaaacatgttggttaaaccccttcctttgctaattaatccttacgctatctccatcattatctccagcttagccctgggaactattactaccctatcaagctaccattgaatgttagcctgaatcggccttgaaattaacactctagcaattattcctctaataactaaaacacctcaccctcgagcaattgaagccgcaactaaatacttcttaacacaagcagcagcatctgccttaattctatttgcaagcacaatgaatgcttgactactaggagaatgagccattaatacccacattagttatattccatctatcctcctctccatcgccctagcgataaaactgggaattgccccctttcacttctgacttcctgaagtcctacaaggattaaccttacaaaccgggttaatcttatcaacatgacaaaaaatcgccccaatagttttacttattcaactatcccaatctgtagaccttaatctaatattattcctcggcttactttctacagttattggcggatgaggaggtattaaccaaacccaaattcgtaaagtcctagcattttcatcaatcgcccacctaggc"

print 'The length of this DNA sequence is', (len(seq)), 'basepairs'

# INSTRUCTIONS: Create and store the RNA equivalent of the sequence, then print to screen.
rna=seq.replace('t','u')
print rna

#INSTRUCTIONS: Create and store the reverse complement of your sequence, then print to screen.

"""
Complement: 
 I couldn't simply directly replace a for t, g for c, and so on, because it would cause confusion
 between the 'real' a's and the a's that had just replaced t's. Thus I do the replacements in two steps,
 first replacing each nucleotide for a placeholder letter, and after that replacing the placeholder
 for the actual complementary nucleotide
 """
   
comp=seq.replace('a', 'b')
comp=comp.replace('t','v')
comp=comp.replace('c', 'd')
comp=comp.replace('g', 'h')
comp=comp.replace('b', 't')
comp=comp.replace('v', 'a')
comp=comp.replace('h', 'c')
comp=comp.replace('d', 'g')

#Reverse:

revcomp=comp[::-1]

print(revcomp)


#INSTRUCTIONS: Extract the bases corresponding to the 13rd and 14th codons from the sequence, then print them to the screen.

"""
Codon 1 corresponds to bases 0, 1 and 2. Codon 2 corresponds to bases 3, 4, and 5.
Generalizing, codon x starts at base 3*(x-1).
Thus, codons 13 and 14 correspond to bases 36 to (including) 41
Therefore:
"""
print (comp[36:42])

"""
INSTRUCTIONS: - Create a function to translate the nucleotide sequence to amino acids 
using the vertebrate mitochondrial genetic code 
"""

def translate(seq):


#this bit of code breaks the input sequence into codons, and stores the codons in a list called "CodonSeq"
    CodonSeq=[]
    codon=""
    for z in seq:
        codon=codon+z
        if len(codon)==3 :
            CodonSeq.append(codon)
            codon=""

#define strings corresponding to the amino acids and to each codon position in the vertebrate mt genetic code:

    aa="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"
    pos0="TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
    pos1="TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
    pos2="TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"



#change the above strings to lower case, since the input sequence should be in lower case too. Ridiculously enough,
#it took me almost an entire day to realize the reason my code was not working was that the cases didn't match.
#So, remember: PYTHON IS CASE-SENSITIVE!!!!

#Note: i write the multi-line comments inside the function with hashtags instead of """ because
#using the """ was, for some reason, causing problems with the indentation.

    aa=aa.lower()
    pos0=pos0.lower()
    pos1=pos1.lower()
    pos2=pos2.lower()


#the code below iterates over each codon in CodonSeq (the codon list created from the input sequence) and
#over the aa's in the genetic code, testing if each base in a codon from CodonSeq corresponds to
#the respective position in the codon for each aa. If all three positions correspond, it appends
#the aa to a string called 'polypeptide', wich is the output of the function


    polypeptide=""
    for i in range (0,len(CodonSeq)) :
        for k in range (0,64) :
            if pos0[k]==CodonSeq[i][0] and pos1[k]==CodonSeq[i][1] and pos2[k]==CodonSeq[i][2] :
                polypeptide=polypeptide+aa[k]
    return polypeptide
    
#INSTRUCTIONS: - Translate the sequence and print it to the screen.
    
#run the funtion I just defined and print the output (a polypeptide):

print(translate(seq))
