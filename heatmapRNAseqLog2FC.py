#!/usr/bin/python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

tempgenesIAC = []
IAC = {}
with open("../DESeq2/filtered_analysis_new/resShrinkIAC.csv","r") as set1:
    for i in set1:
        i = i.rstrip().replace("\"","")
        if "baseMean" not in i:
            i = i.split("\t")
            if i[5] != "NA":
                padj = float(i[5])
                log2 = round(float(i[2]),2)
                gen = i[0]
                IAC[gen] = [log2,padj]
                tempgenesIAC.append(gen)

tempgenesSP = []
SP = {}
with open("../DESeq2/filtered_analysis_new/resShrinkSP.csv","r") as set1:
    for i in set1:
        i = i.rstrip().replace("\"","")
        if "baseMean" not in i:
            i = i.split("\t")
            if i[5] != "NA":
                padj = float(i[5])
                log2 = round(float(i[2]),2)
                gen = i[0]
                SP[gen] = [log2,padj]
                tempgenesSP.append(gen)


tempgenesRB = []
RB = {}
with open("../DESeq2/filtered_analysis_new/resShrinkRB.csv","r") as set1:
    for i in set1:
        i = i.rstrip().replace("\"","")
        if "baseMean" not in i:
            i = i.split("\t")
            if i[5] != "NA":
                padj = float(i[5])
                log2 = round(float(i[2]),2)
                gen = i[0]
                RB[gen] = [log2,padj]
                tempgenesRB.append(gen)


treatments = ["IAC", "SP", "RB"] #treatments

genesIACSPRB = set(tempgenesIAC).intersection(tempgenesSP)
genesIACSPRB = genesIACSPRB.intersection(tempgenesRB)
genesIACSPRB = list(genesIACSPRB)

pcut = 0.05
DEGS = []
for i in genesIACSPRB:
    if float(IAC[i][1]) <= pcut or float(SP[i][1]) <= pcut or RB[i][1] <= pcut:
        DEGS.append(i)

print (len(DEGS),"nplots:",len(DEGS)/300)

arrays = []
geneis = []
for i in range(0,len(DEGS),50):
    addgeneis = []
    addarrays = []

    if i+50 < len(DEGS):
        maxi = i+50
    else:
        maxi = len(DEGS)

    for j in range(i,maxi):
        addgeneis.append(DEGS[j])
        addarrays.append([IAC[DEGS[j]][0],SP[DEGS[j]][0],RB[DEGS[j]][0]])
    arrays.append(addarrays)
    geneis.append(addgeneis)

################################################################################
countplot = 0
for reg in range(0,len(geneis),6):
    countplot += 1
    if reg + 6 < len(geneis):
        maxj = reg + 6
    else:
        maxj = len(geneis)

    addfig = 0
    fig = plt.figure(figsize=(18,10)) # width vs height in inches
    fig.subplots_adjust(hspace=0.28, wspace=0.39, left=0.08, bottom=0.13, right=0.98, top=0.97)
    myaxes = []
    for i in range(reg,maxj):
        addfig += 1
        ax = fig.add_subplot(1, 6, addfig) #two first numbers: rows vs columns
        myaxes.append(ax)
        im = ax.imshow(np.array(arrays[i]))
        # We want to show all ticks...
        ax.set_xticks(np.arange(len(treatments)))
        ax.set_yticks(np.arange(len(geneis[i])))
        # ... and label them with the respective list entries
        ax.set_xticklabels(treatments,fontsize=10)
        ax.set_yticklabels(geneis[i],fontsize=8)
        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=60, ha="right",rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        #for i in range(len(genes)):
        #    for j in range(len(treatments)):
        #        text = ax.text(j, i, harvest[i, j],
        #                       ha="center", va="center", color="w")

        ax.set_title("")
        fig.tight_layout()

    plt.colorbar(im, ax=myaxes, orientation="horizontal", fraction=0.02, pad=0.07,label='Log2 Fold Change\n(inoculated/control)')
    #plt.show()
    plt.savefig('plot_'+str(countplot)+'.png',dpi=100,bbox_inches='tight')
    plt.close(fig)
