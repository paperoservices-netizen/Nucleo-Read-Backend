import matplotlib.pyplot as plt

def virtual_gel(cuts,length,out):
    frags=[]
    last=0
    for c in sorted(cuts):
        frags.append(c-last)
        last=c
    frags.append(length-last)

    y=range(len(frags))
    plt.scatter(frags,y)
    plt.yticks([])
    plt.xlabel("Fragment size (bp)")
    plt.title("Virtual Gel")
    plt.savefig(out)
    plt.close()
