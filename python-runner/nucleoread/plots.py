import matplotlib.pyplot as plt

def save_gc_plot(seq, out_path):
    window = 50
    seq = seq.replace("U", "T")

    if len(seq) < window:
        return None

    values = [
        (seq[i:i+window].count("G") + seq[i:i+window].count("C")) / window * 100
        for i in range(len(seq) - window)
    ]

    plt.figure(figsize=(9, 3))
    plt.plot(values)
    plt.title("GC Content (Sliding Window)")
    plt.xlabel("Position")
    plt.ylabel("GC %")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()

    return out_path
