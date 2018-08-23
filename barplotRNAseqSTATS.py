#!/usr/bin/python3

import math
import string
import itertools
import numpy as np
import pandas as pd
import itertools as it
from scipy import stats
import matplotlib.pyplot as plt
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison

tukeyout = open("tukey_output.txt","w")
myletters = list(string.ascii_lowercase)
cutoff = 0.05
barletters = {} # {'Sample1': {1: ['a'], 2: ['a', 'b'], 3: ['b'], 4: ['d'],n...},'Sample2':...}
#Filter susceptible resistant
#function to check if any element in a list a is also in a list b
def common_member(a, b):
    a_set = set(a)
    b_set = set(b)
    if len(a_set.intersection(b_set)) > 0:
        return(True)
    return(False)

################################################################################
#Create dict DEGS with DEG as keys
DE = pd.read_csv("/home/hugo/Dropbox/Esalq/Flowering_RNAseq_CLC/Flowering_RNA_SEQ_filtered_Table1_test.csv")
DEGS = {}
for i in range(len(DE)):
    if float(DE.loc[i]['padj']) <= cutoff:
        DEGS[DE.loc[i]['Reference']] = ''
del DE #release memory from pandas dataframe
print ("Total number of DEGs:",len(DEGS))

################################################################################
#Create dataframe for Log2 Fold Change file
df1 = pd.read_csv("log2foldchange.csv", sep='\t') #read the csv file
print ("Size of Log2 Fold Change dataframe:",len(df1))

#remove rows from dataframe df1 if are not identified as DE
listdrop = []
for i in range(len(df1)):
    identifier = df1.loc[i][0]
    if identifier not in DEGS:
        listdrop.append(i)
print ("Number of indexes not DE dropped:",len(listdrop))

#Reorganize dataframe df1
df1.drop(df1.index[listdrop], inplace=True) #drop from dataframe
df1.reset_index(inplace=True,drop=True) #reset index of rows

#release memory
del DEGS
del listdrop

#lists to feed dataframe df1
listctrmeans = []
listctrstdmeans = []
listtreatmeans = []
listtreatstdmeans = []
listttestpairs = []

#STATISTICS
fvalues = [] #fvalue for each index of dataframe
pvalues = [] #pvalue for each index of dataframe
for i in range(len(df1)): #for each index in dataframe
    identifier = df1.loc[i][0]
    count = 0
    listanova = [] #List with the values of N listanova [[val1,val2,val3],[val1,val2,val3]]
    counttreats = {} #Dict {Treat number:[val1,val2,val3]}
    tukeydata = []
    tukeygroups = []
    for j in range(1,len(df1.columns),3): #read from 3 to 3 columns, considering treats are 3 replicates
        listanova.append([df1.loc[i][j],df1.loc[i][j+1],df1.loc[i][j+2]])
        counttreats[count] = [df1.loc[i][j],df1.loc[i][j+1],df1.loc[i][j+2]]
        count += 1
        tukeydata.append(df1.loc[i][j])
        tukeydata.append(df1.loc[i][j+1])
        tukeydata.append(df1.loc[i][j+2])
        for group in range(3):
            tukeygroups.append(count)

    f_value, p_value = stats.f_oneway(*listanova) #compute ANOVA among all listanova
    fvalues.append(f_value) #append ANOVA fvalue for the index
    pvalues.append(p_value) #append ANOVA pvalue for the index

    ###################################################################### TUKEY
    mc = MultiComparison(tukeydata, tukeygroups)
    result = mc.tukeyhsd()
    #print (identifier,"\n",result)
    ############################################################################
    tukeyout.write(identifier + "\n" + str(result) + "\n")
    ############################################################################
    countt = 0
    treatments = [] #list of treatments [1,2,3,4,n...]
    for ix in listanova:
        countt += 1
        treatments.append(countt)

    groupPAIR = list(it.combinations(treatments, 2)) #all possible treatments pairs
    groupRELATION = list(result.reject) # True, the two group's means ARE significantly different
                                        # False, the two group's means ARE NOT significantly different
    treat_letter = {} #letters in each of treatments {1:['a','b'],2:['a'],3:['b'],n...}
    countl = 0 #count letters
    combination_blacklisted = [] #black list treatments pairs

    elemletters = {}
    cl = 0
    for u in range(len(groupPAIR)):
        ele1 = groupPAIR[u][0]
        ele2 = groupPAIR[u][1]
        if str(groupRELATION[u]) == "False":
            if ele1 in elemletters and ele2 in elemletters:
                equal_letter = any(elem in elemletters[ele1] for elem in elemletters[ele2])
                #print (ele1,elemletters[ele1],ele2,elemletters[ele2],groupRELATION[u],equal_letter)
                if equal_letter == False:
                    cl += 1
                    ele12l = myletters[cl]
                    add1 = elemletters[ele1]
                    add2 = elemletters[ele2]
                    add1.append(ele12l)
                    add2.append(ele12l)
                    elemletters[ele1] = add1
                    elemletters[ele2] = add2
            elif ele1 not in elemletters and ele2 not in elemletters:
                elemletters[ele1] = [myletters[cl]]
                elemletters[ele2] = [myletters[cl]]
            elif ele1 in elemletters and ele2 not in elemletters:
                ele2l = elemletters[ele1][0]
                elemletters[ele2] = [ele2l]
        else: #if relation is True
            if ele1 in elemletters and ele2 in elemletters:
                if elemletters[ele1] == elemletters[ele2]: #se forem diferentes nao precisa add
                    cl += 1
                    elemletters[ele1] = [myletters[cl]]
                    for uu in range(len(groupPAIR)):
                        if str(groupRELATION[uu]) == "False" and ele1 in groupPAIR[uu]:
                            for elem in groupPAIR[uu]:
                                if elem != ele1:
                                    add1 = elemletters[elem]
                                    add1.append(myletters[cl])
                                    elemletters[elem] = add1
            if ele1 in elemletters and ele2 not in elemletters:
                cl += 1
                elemletters[ele2] = [myletters[cl]]
            if ele1 not in elemletters and ele2 in elemletters:
                cl += 1
                elemletters[ele1] = [myletters[cl]]
            if ele1 not in elemletters and ele2 not in elemletters:
                elemletters[ele1] = [myletters[cl]]
                cl += 1
                elemletters[ele2] = [myletters[cl]]
    #print (elemletters)
    barletters[identifier] = elemletters
    ############################################################################
    # Conduct t-test on each pair
    pairs = list(itertools.combinations(list(counttreats.keys()), 2)) #compute all possible pairs for N listanova
    bonferroni = cutoff / len(pairs) #adjust for this multiple comparison problem by dividing the statistical significance level by the number of comparisons made. In this case, if we were looking for a significance level of 5%, we'd be looking for p-values of 0.05/10 = 0.005 or less. This simple adjustment for multiple comparisons is known as the Bonferroni correction.
    letters = {}
    count = 0

    for trat1, trat2 in pairs:
        ttest = stats.ttest_ind(counttreats[trat1],counttreats[trat2]) #compute t-test
        pvaluettest = ttest[1]
        if float(p_value) <= cutoff and float(pvaluettest) <= cutoff: #significance level
            #print (df1.loc[i]['genes'],trat1,trat2)

            if trat1 not in letters:
                letters[trat1] = [myletters[count]]
            else:
                addlet = letters[trat1]
                addlet.append(myletters[count])
                letters[trat1] = addlet

            if trat2 not in letters:
                letters[trat2] = [myletters[count]]
            else:
                addlet2 = letters[trat2]
                addlet2.append(myletters[count])
                letters[trat2] = addlet2
            count += 1

    finalletters = []
    for x in range(len(counttreats)):
        if x not in letters:
            finalletters.append([])
        else:
            finalletters.append(letters[x])
    listttestpairs.append(finalletters)

    ############################################################################
    controls = []
    stdctr = []
    treats = []
    stdtre = []
    #compute controls
    for jj in range(0,len(listanova),2):
        controls.append(np.mean(listanova[jj]))
        stdctr.append(np.std(listanova[jj]))
    for jjj in range(1,len(listanova),2):
        treats.append(np.mean(listanova[jjj]))
        stdtre.append(np.std(listanova[jjj]))

    listctrmeans.append(tuple(controls))
    listctrstdmeans.append(tuple(stdctr))
    listtreatmeans.append(tuple(treats))
    listtreatstdmeans.append(tuple(stdtre))

#feed dataframe
df1.insert(len(df1.columns),'fvalue', fvalues )
df1.insert(len(df1.columns),'pvalue', pvalues )
df1.insert(len(df1.columns),'ctrmeans', listctrmeans )
df1.insert(len(df1.columns),'ctrstdmeans', listctrstdmeans )
df1.insert(len(df1.columns),'treatmeans', listtreatmeans )
df1.insert(len(df1.columns),'treatstdmeans', listtreatstdmeans )
df1.insert(len(df1.columns),'ttestpairs', listttestpairs )

################################################################################
#FILTERS

#remove from dataframe if ANOVA pvalue not significant
anovadroplist = []
for i in range(len(df1)):
    if float(df1.loc[i]['pvalue']) > cutoff:
        anovadroplist.append(i)

#Reorganize dataframe df1
df1.drop(df1.index[anovadroplist], inplace=True) #drop from dataframe
df1.reset_index(inplace=True,drop=True) #reset index of rows
del anovadroplist

#FILTER 1
susresistdroplist = []
for i in range(len(df1)):
    a = df1.loc[i]['ttestpairs'][1] #treatment1 letters
    b = df1.loc[i]['ttestpairs'][3] #treatment2 letters
    c = df1.loc[i]['ttestpairs'][2] #controle2 letters
    comonab = common_member(a,b)
    trat_a = float(df1.loc[i]['treatmeans'][0]) #treatment1 means
    trat_b = float(df1.loc[i]['treatmeans'][1]) #treatment2 means
    trat_c = float(df1.loc[i]['treatmeans'][2]) #treatment2 means
    ctr_a = float(df1.loc[i]['ctrmeans'][0])
    ctr_b = float(df1.loc[i]['ctrmeans'][1])
    ctr_c = float(df1.loc[i]['ctrmeans'][2])
    if float(trat_a) < float(ctr_a) and float(trat_b) < float(ctr_b) and float(trat_c) < float(ctr_c) or float(ctr_a) > 0.0 or float(ctr_b) > 0.0 or float(ctr_c) > 0.0: #conditions to remove from dataframe
        susresistdroplist.append(i)

        print (trat_a,ctr_a,trat_b,ctr_b,trat_c,ctr_c)

#Reorganize dataframe df1
#DESLIGA O FILTRO COMENTANDO A LINHA ABAIXO
df1.drop(df1.index[susresistdroplist], inplace=True) #drop from dataframe
df1.reset_index(inplace=True,drop=True) #reset index of rows
del susresistdroplist

################################################################################
#Write a dataframe file
wdf = open("dataframe.txt","w")
wdf.write(df1.to_string())
wdf.close()
print ("Final dataframe size:",len(df1))
################################################################################
#PARAMETERS FOR THE PLOT
N = 3
ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
################################################################################
#PLOT
print ("Number of barplots:",round(len(df1)/6.0))
countplot = 0 #number of the plot
for nplot in range(0,len(df1),6): #read df1 from row six to six records
    addn = 0
    fig = plt.figure(figsize=(18,10)) # width vs height in inches
    fig.subplots_adjust(hspace=0.4, wspace=0.4, left=0.07, bottom=0.04, right=0.95, top=0.96)
    labels = ('IAC', 'SP', 'RB')
    countplot += 1
    finalrange = nplot + 6
    if finalrange > len(df1):
        finalrange = len(df1)

    for row in range(nplot,finalrange):
        addn += 1
        ax = fig.add_subplot(2, 3, addn) #two first numbers: rows vs columns
        xval1 = df1.loc[row]['ctrmeans']
        xyerr1 = df1.loc[row]['ctrstdmeans']
        xval2 = df1.loc[row]['treatmeans']
        xyerr2 = df1.loc[row]['treatstdmeans']
        ylim = max(xval1 + xval2) + max(xyerr1 + xyerr2) + 2.0

        if float(df1.loc[row]['pvalue']) <= 0.05: #if ANOVA is significant, print * on title
            anoval = "*"
        else:
            anoval = ""

        rects1 = ax.bar(ind + .65, df1.loc[row]['ctrmeans'], width,
                color='#999999', yerr=df1.loc[row]['ctrstdmeans'], ecolor='black') #control

        rects2 = ax.bar(ind + .65 + width, df1.loc[row]['treatmeans'], width,
                color='#ef8a62', yerr=df1.loc[row]['treatstdmeans'], ecolor='black') #treatments

        ax.set_alpha(0.8)
        ax.set_xticks(ind + .65 + width)
        ax.set_xticklabels(labels,rotation=20,fontsize=12) #config labels
        ax.set_ylabel('Log2 Fold Change',fontsize=12)
        ax.set_ylim([0,ylim])
        gene = df1.loc[row]['genes']
        ax.set_title(gene + anoval, horizontalalignment='center', fontname='Verdana',
                            fontsize=14, fontstyle='normal', fontweight='bold')

        count = 0
        countlet = 1
        for i in rects1:
            height = i.get_height()
            x = i.get_x() #get position in x axis
            ax.text(x + width/2, (xyerr1[count] + height + (0.02 * (ylim - xyerr1[count] + height))),
                str(round(height, 2)) + "\n" + "".join(barletters[gene][countlet]),
                horizontalalignment ='center',fontsize=12)
            count += 1
            countlet += 2

        count = 0
        countlet = 2
        for i in rects2:
            height = i.get_height()
            x = i.get_x()
            ax.text(x + width/2, (xyerr2[count] + height + (0.02 * (ylim - xyerr2[count] + height))),
                str(round(height, 2)) + "\n" + "".join(barletters[gene][countlet]),
                horizontalalignment ='center',fontsize=12)
            count += 1
            countlet += 2

    fig.legend((rects1[0], rects2[0]), ('Control', 'Treatment'),'lower right')
    plt.savefig('plot_'+str(countplot)+'.png',dpi=100,bbox_inches='tight')
    #plt.show()
    plt.close(fig) #release pyplot memory after each figure plot is saved
