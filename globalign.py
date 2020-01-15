"""
Ronald Randolph
CS594: Bioinformatics - HW2 Resubmission
This program computes the global alignment score of two sequences
using the following parameters: 

+2 for a match, -1 for a mismatch, and -2 for a gap. 

After, it performs the traceback and displays the alignment in a 
user-friendly format.
"""
import sys
import numpy as np

#function reads and stores sequences for alignment
def get_sequences():
    
    with open(sys.argv[1]) as f:
        data1 = f.read().splitlines()
        data1.pop(0)

    with open(sys.argv[2]) as f:
        data2 = f.read().splitlines()
        data2.pop(0)

    seq1 = "".join(data1)
    seq2 = "".join(data2)

    return seq1, seq2

#function creates matricices for sequence alignment
def initialize(n,m):
    
    #create substitution and traceback matrices
    T = np.zeros((n,m), dtype=int)
    A = np.zeros((n,m), dtype=int)

    T[0,0] = 0
    A[0,0] = -1

    #initialize matrix for global alignment 
    for j in range(1,m):
        T[0][j] = (T[0][j-1] - 2)
        A[0][j] = -1
    
    for i in range(1,n):
        T[i][0] = (T[i-1][0] - 2)
        A[i][0] = -1
    
    return T,A

#function fills substitution and traceback matrices
def fill(T,A,s,t):
    f = [0,0,0]

    #Fill substitution matrix
    for i in range(1,len(s)+1):
        for j in range(1,len(t)+1):
            f[0] = (T[i-1][j-1] + score(s[i-1],t[j-1]))
            f[1] = (T[i-1][j] - 2)
            f[2] = (T[i][j-1] - 2)
            
            #record traceback 
            T[i][j] = max(f)
            A[i][j] = np.argmax(f)
    
    #return score and filled matrices
    return T,A,T[len(s)][len(t)]

#function returns the score of matching two nucleotides
def score(a,b):
    if a == b:
        return 2
    else:
        return -1

#function prints the optimal alignment of sequences
def print_alignment(A,s,t):

    i = len(s)
    j = len(t)
    
    stmp = ""
    ttmp = ""
    ptmp = ""
    S = ""
    P = ""
    T = ""

    #perfrom traceback
    while(A[i][j] != -1):

        #diagonal
        if(A[i][j] == 0):
            stmp += s[i-1]
            ttmp += t[j-1]
            ptmp += "|"
            i = i-1
            j = j-1
        #up
        elif(A[i][j] == 1):
            stmp += s[i-1]
            ttmp += "-"
            ptmp += " "
            i = i-1
        #back
        elif(A[i][j] == 2):
            stmp += "-"
            ttmp += t[j-1]
            ptmp += " "
            j = j-1
   
    index = len(stmp)
    nl = 0
    
    #print out strings of aligned sequences 
    while index >= 0:
        if(index > 0):
            S += stmp[index-1] + " "
            T += ttmp[index-1] + " "
            P += ptmp[index-1] + " "
            nl += 1
        
        if((nl % 40 == 0 and index != 0) or (index == 0 and S != "")):
            print("\ns:",S)
            print("  ",P)
            print("t:",T)
            S = ""
            T = ""
            P = ""

        index = index-1

def main():
    
    #argument check
    if(len(sys.argv) != 3):
        print("Usage: python3 globalign.py seq1.fasta seq2.fasta")
        sys.exit()
    
    s,t = get_sequences()

    n = len(s)+1
    m = len(t)+1

    T,A = initialize(n,m)
    T,A,score = fill(T,A,s,t)
    print_alignment(A,s,t)
    
    print("\nOptimal Alignment Score:",score)
    print()

if __name__ == "__main__":
    main()