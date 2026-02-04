import matplotlib.pyplot as plt

def gc_plot(seq,out):
    window=50
    vals=[(seq[i:i+window].count("G")+seq[i:i+window].count("C"))/window*100
          for i in range(len(seq)-window)]
    plt.figure(figsize=(10,4))
    plt.plot(vals)
    plt.savefig(out,dpi=150)
    plt.close()
