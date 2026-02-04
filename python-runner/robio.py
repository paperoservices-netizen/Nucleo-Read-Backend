import re

START = "ATG"
STOPS = {"TAA","TAG","TGA"}

CODON = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

ENZYMES={"HaeIII":"GGCC","EcoRI":"GAATTC","BamHI":"GGATCC"}

def parse_fasta(text):
    lines=text.strip().splitlines()
    return "".join(lines[1:]).upper().replace("U","T")

def reverse_complement(d):
    return d.translate(str.maketrans("ATGC","TACG"))[::-1]

def find_orf(d):
    best=None
    for frame in range(3):
        for i in range(frame,len(d)-2,3):
            if d[i:i+3]==START:
                for j in range(i+3,len(d)-2,3):
                    if d[j:j+3] in STOPS:
                        l=j+3-i
                        if not best or l>best[2]:
                            best=(i,j+3,l)
                        break
    return best

def translate(seq):
    return "".join(CODON.get(seq[i:i+3],"X") for i in range(0,len(seq)-2,3))

def restriction_map(dna):
    sites={}
    for k,v in ENZYMES.items():
        sites[k]=[m.start()+1 for m in re.finditer(v,dna)]
    return sites

def analyze(fasta,job):
    dna=parse_fasta(fasta)
    if len(dna)==0: raise Exception("Empty FASTA")

    rc=reverse_complement(dna)
    gc=round((dna.count("G")+dna.count("C"))/len(dna)*100,2)

    orf=find_orf(dna)
    protein=""
    if orf:
        protein=translate(dna[orf[0]:orf[1]])

    return {
        "job_id":job,
        "length":len(dna),
        "gc":gc,
        "rna":dna.replace("T","U"),
        "complement":dna.translate(str.maketrans("ATGC","TACG")),
        "reverse_complement":rc,
        "orf":[orf[0],orf[1]] if orf else None,
        "protein":protein,
        "enzymes":restriction_map(dna),
        "dna":dna
    }
