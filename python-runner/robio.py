import re

def clean_fasta(fasta):
    lines = fasta.splitlines()
    seq = "".join([l.strip() for l in lines if not l.startswith(">")])
    return seq.upper()

def dna_to_rna(dna):
    return dna.replace("T","U")

def complement(dna):
    return dna.translate(str.maketrans("ATCG","TAGC"))

def reverse_complement(dna):
    return complement(dna)[::-1]

def gc_content(dna):
    return round((dna.count("G")+dna.count("C"))/len(dna)*100,2)

def find_orfs(seq):
    stops = ["TAA","TAG","TGA"]
    orfs=[]
    for frame in range(3):
        for i in range(frame, len(seq)-2,3):
            if seq[i:i+3]=="ATG":
                for j in range(i,len(seq)-2,3):
                    if seq[j:j+3] in stops:
                        orfs.append((i,j+3))
                        break
    return orfs

codon_table = {
 "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
 "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
 "ATG":"M","TAA":"*","TAG":"*","TGA":"*"
}

def translate(seq):
    prot=""
    for i in range(0,len(seq)-2,3):
        prot+=codon_table.get(seq[i:i+3],"X")
    return prot

def analyze(fasta, job_id):
    dna=clean_fasta(fasta)
    rna=dna_to_rna(dna)

    orfs=find_orfs(dna)
    best = max(orfs, key=lambda x:x[1]-x[0]) if orfs else None

    protein=""
    if best:
        protein = translate(dna[best[0]:best[1]])

    return {
        "job_id": job_id,
        "length": len(dna),
        "gc": gc_content(dna),
        "rna": rna,
        "complement": complement(dna),
        "reverse_complement": reverse_complement(dna),
        "orf": best,
        "protein": protein
    }
