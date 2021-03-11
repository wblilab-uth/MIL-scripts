import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import rcParams
import seaborn as sns
from decimal import Decimal

rcParams["font.size"] = 8
rcParams["font.family"] = "Arial"
rcParams['pdf.fonttype'] = 42


def super_mils(ml1df, title, plotname):
    y = np.interp(ml1df["meth"], (ml1df["meth"].min(), ml1df["meth"].max()), (0, +1))
    x = np.interp(ml1df["meth_rank"], (ml1df["meth_rank"].min(), ml1df["meth_rank"].max()), (0, +1))
    min_x = x.shape[0]
    min_idx = []
    min_k = 0

    for k in np.arange(0, 10000):
    
        f = x - k / float(10000)

        idx = np.argwhere(np.diff(np.sign(f - y))).flatten()

        if idx != []:
            if (np.max(idx) - np.min(idx[0])) < min_x:
                min_x = (np.max(idx) - np.min(idx[0]))
                min_idx = idx
                min_k = k

    fig, ax = plt.subplots(figsize=(3, 3.5))

    y = np.interp(ml1df["meth"], (ml1df["meth"].min(), ml1df["meth"].max()), (0, +1))
    x = np.interp(ml1df["meth_rank"], (ml1df["meth_rank"].min(), ml1df["meth_rank"].max()), (0, +1))
    f = x - min_k / float(10000)
    ax.plot(x, f, '-', linewidth=2)
    ax.plot(x, y, '-', color="#E55050", linewidth=2)
    ax.set_ylim([-0.1, 1.1])
    idx = np.argwhere(np.diff(np.sign(f - y))).flatten()
    ax.plot(x[idx], f[idx], 'go')
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(1)
    old_x = ml1df["meth_rank"]
    old_y = ml1df["meth"]
    ax.set_xticks(([x.min(), np.quantile(x, 0.5), x.max()]))
    ax.set_yticks(([y.min(), 0.5 * y.max(), y.max()]))

    ax.set_xticklabels(np.round([old_x.min() - 1, np.quantile(old_x, 0.5) - 1, old_x.max() - 1]).astype(int))
    ax.set_yticklabels((np.round([old_y.min(), 0.5 * old_y.max(), old_y.max()])))
    ax.set_xlabel("methylation ranking")
    ax.set_ylabel("Relative m6A level")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(plotname + ".pdf", dpi=300, transparent=True)
    return min_idx, min_k



###Peak based MILs identification and characteriztion:
###Input:1). peak overlapped intron L1 bed, 2). MINT-Seq,TT-Seq,RNA-Seq quantification 3). L1 annotation file; 4). RefSeq gene homer annotation file
pdf = pd.read_csv("/path/to/data/peak.merged.overlapped.L1.bed",sep = "\t",names = ["chr","start","end","name","length","strand"])
anno = pd.read_csv("/path/to/data/L1.intron.anno.txt",sep="\t",header=None)
pdf = pd.merge(pdf,anno,left_on = ["chr","start","end"],right_on = [1,2,3])
df = pd.read_csv("/path/to/data/MINT-TT.L1.rpkm.txt",sep = "\t")


###TT-Seq cutoff, and L1 length cutoff
df = df[df.iloc[:,5]>200]
df = df[df.iloc[:,0].isin(pdf.iloc[:,3])]

l1df = df[df.iloc[:,8]>0.1]
l1df = l1df[l1df.iloc[:,0].isin(anno.iloc[:,0])]
psedo_count = 1
l1df = l1df[l1df.iloc[:,0].isin(agdf.iloc[:,0])]
l1df["meth"] = (l1df.iloc[:,9]+psedo_count)/(l1df.iloc[:,8]+psedo_count)
ml1df["meth_rank"] = ml1df['meth'].rank(ascending = True)
ml1df = ml1df.sort_values(["meth"]).reset_index()

###Idenfity super-MILs
min_x,min_k =super_mils(ml1df,"Super-MILs","FigX-Super-MILs-ranking")
s_bool = []
for i in ml1df.index:
    if ml1df.loc[i,"meth_rank"] > np.mean(min_x):
        s_bool.append("True")
    else:
        s_bool.append("False")

ml1df["super"] = s_bool
ml1df[ml1df["super"] == "False"].to_csv('./Typical-MILs.txt',sep = "\t")
ml1df[ml1df["super"] == "True"].to_csv('./Super-MILs.txt',sep = "\t")

###Annotating host genes of MILs
adf = pd.merge(l1df,anno,left_on = l1df.columns[0],right_on = anno.columns[0])
adf["pos"] = adf.apply( lambda i: i[7].split("(")[0],axis =1)
adf = adf[adf["pos"] == "intron "]
adf["gene"] = adf.apply( lambda i: i[7].split("(")[1].split(",")[0],axis =1)
gdf = pd.read_csv("/path/to/data/gene-annotation.txt",sep ="\t").iloc[:,[0,1,2,3,4,7,8,9,10,11]]
agdf = pd.merge(adf,gdf,left_on= "gene", right_on = gdf.columns[0])
agdf["gene"].to_csv("./mil.host.gene.txt",sep= "\t",index = False)