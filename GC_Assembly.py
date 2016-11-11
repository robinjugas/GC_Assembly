#Dependencies
from Bio import SeqIO #fasta loading
#import cmath #complex numbers
import numpy as np 
import matplotlib.pyplot as plt

#Open FASTA file new variable
lengths=list()
sequences=list()
handle=open("Ostreococcus tauri_kontigy.fasta","r")
for record in SeqIO.parse(handle, "fasta"):
    print(record.id)
    print(len(record.seq))
    lengths.append(len(record.seq))
    sequences.append(str(record.seq))

num_seq=list()
#Convert to genomic signals
for x in range(len(sequences)): 
    temp=np.zeros(lengths[x])
    for y in range(lengths[x]):
        if sequences[x][y]=="A":
            temp[y]=np.angle(1+1j)
        elif sequences[x][y]=="C":
            temp[y]=np.angle(-1-1j)
        elif sequences[x][y]=="G":
            temp[y]=np.angle(-1+1j)
        elif sequences[x][y]=="T":
            temp[y]=np.angle(1-1j)
        else:
            temp[y]=0
    num_seq.append(np.cumsum(temp))

for x in range(len(sequences)): 
    X=range(np.size(num_seq[x]))        
    plt.plot(X,num_seq[x])
plt.show()
#pokraovani
    




