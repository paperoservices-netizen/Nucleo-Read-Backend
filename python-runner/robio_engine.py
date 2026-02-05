import re
import matplotlib.pyplot as plt
import numpy as np
import os

START_CODON="ATG"
STOP_CODONS={"TAA","TAG","TGA"}

GENETIC_CODE={
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
 'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
 'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
 'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
 'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
 'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

RESTRICTION_ENZYMES={
 "EcoRI":"GAATTC","BamHI":"GGATCC","HindIII":"AAGCTT","NotI":"GCGGCCGC",
 "XhoI":"CTCGAG","KpnI":"GGTACC","AluI":"AGCT","HaeIII":"GGCC","TaqI":"TCGA"
}

def parse_fasta(text):
    lines=text.strip().splitlines()
    seq="".join(l for l in lines if not l.startswith(">")).upper()
    if not seq or not re.fullmatch("[ATUGCN]+",seq):
        raise ValueError("Invalid FASTA")
    return seq.replace("U","T")

def find_orfs(seq):
    orfs=[]
    for frame in range(3):
        for i in range(frame,len(seq)-2,3):
            if seq[i:i+3]==START_CODON:
                for j in range(i+3,len(seq)-2,3):
                    if seq[j:j+3] in STOP_CODONS:
                        orfs.append((i,j+3))
                        break
    return orfs

def translate(seq):
    return "".join(GENETIC_CODE.get(seq[i:i+3],"X") for i in range(0,len(seq)-2,3))

def plot_gc(seq,out):
    window=50
    vals=[(seq[i:i+window].count("G")+seq[i:i+window].count("C"))/window*100 for i in range(len(seq)-window)]
    plt.figure(figsize=(8,3))
    plt.plot(vals)
    plt.savefig(out,dpi=150)
    plt.close()

def draw_gel(length,cuts,out):
    fragments=[length]
    if cuts:
        cuts=sorted(cuts)
        fragments=[cuts[0]]+[cuts[i+1]-cuts[i] for i in range(len(cuts)-1)]+[length-cuts[-1]]
    plt.figure(figsize=(4,6))
    plt.bar([1]*len(fragments),fragments)
    plt.savefig(out,dpi=150)
    plt.close()

def analyze_fasta(fasta, job):
    os.makedirs("python-runner/images", exist_ok=True)
    os.makedirs("python-runner/results", exist_ok=True)

    dna=parse_fasta(fasta)
    gc=(dna.count("G")+dna.count("C"))/len(dna)*100
    rna=dna.replace("T","U")
    comp=dna.translate(str.maketrans("ATGC","TACG"))
    rev=comp[::-1]

    orfs=find_orfs(dna)
    start,end=max(orfs,key=lambda x:x[1]-x[0])
    protein=translate(dna[start:end])

    cuts=[]
    for enz,site in RESTRICTION_ENZYMES.items():
        for m in re.finditer(site,dna):
            cuts.append(m.start())

    plot_gc(dna,f"python-runner/images/{job}_gc.png")
    draw_gel(len(dna),cuts,f"python-runner/images/{job}_gel.png")

    return {
        "length":len(dna),
        "gc":round(gc,2),
        "rna":rna,
        "complement":comp,
        "reverse_complement":rev,
        "orf":[start,end],
        "protein":protein
    }
