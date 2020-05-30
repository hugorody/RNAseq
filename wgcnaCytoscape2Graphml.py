#!/usr/bin/python3

import networkx as nx
import glob, os
import re

################################################################################
# INPUTS

dir1 = "/PATH/MODULES/" #directory of cytoscape input graphs
diroutput = "/PATH/OUTPUTS/"
functionalannotation = "compgg_functionalannotation_HR28may20.tsv"
targetsfile = "TARGETcandidateGENES.txt" #target genes separated by breakline

################################################################################
# PARSE GRAPH INPUTS
os.chdir(dir1) #list all files in the input directory 1

listofnodefiles = []
for file in glob.glob("*nodes*"):
    listofnodefiles.append(file)

################################################################################
# PARSE FUNCTIONAL ANNOTATION
compggfunc = {}
bh_ath = {}
with open(functionalannotation,"r") as set1:
    for i in set1:
        i = i.rstrip().split("\t")
        if "#" not in i[0]:
            transcriptid = i[0]
            orfid = i[1]
            Orthogroup = i[2]
            MAPKs = i[3]
            TFs = i[4]
            RGA = i[5]
            IPR = i[6]
            Pfam = i[7]
            TMHMM2 = i[8]
            TargetP2 = i[9]
            SignalP5 = i[10]
            BesthitAth = re.sub(".[0-9]$","",i[11])
            BesthitAthannotation = i[12]
            #BesthitZmays = i[13]
            SP_UPDOWN = i[14]
            SP_LFC = i[15]
            IAC_UPDOWN = i[16]
            IAC_LFC = i[17]
            RB05_UPDOWN = i[18]
            RB05_LFC = i[19]
            RB200_UPDOWN = i[20]
            RB200_LFC = i[21]
            Blast2GO = i[22]
            GO = i[23]
            KOG = i[24]
            #Transcriptseq = i[25]
            #ORFseq = i[26]

            #FEED DICTIONARIES
            compggfunc[transcriptid] = i
            bh_ath[transcriptid] = [BesthitAth,BesthitAthannotation]


################################################################################
# PARSE TARGET GENES
targetgenes = {}
with open(targetsfile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        targetgenes[i] = ''

################################################################################
# PARSE WGCNA GRAPHS

targetmodules = {}
Gwgcna = {}
for file in listofnodefiles:
    with open(file,"r") as set1:
        for i in set1:
            i = i.rstrip().split("\t")
            if "nodeName" not in i[0]:
                node = i[0]
                color = i[2]
                Gwgcna[node] = color

                if node in targetgenes: #FILTER FOR ONLY TARGETED GENES
                    targetmodules[color] = re.sub("-nodes-","-edges-",file) #select module to parse graph
                    targetgenes[node] = color #update targetgenes dict

################################################################################
print ("Targeted modules:",", ".join([m for m in targetmodules.keys()]))
print ("Genes and modules:")
for i in targetgenes.items():
    print (i[0],i[1])
################################################################################
# GRAPH MODELING

#CREATE ONE SEPARATED GRAPH FOR EACH OF TARGETED MODULES
for m in targetmodules.items():
    os.chdir(dir1)
    targetmodulecolor = m[0]
    targetedgefile = m[1]

    print ("Modeling:",targetmodulecolor)

    g = nx.Graph() # Define the undirected graph and type
    for i in Gwgcna.items():
        node = i[0]
        modulecolor = i[1]
        if modulecolor == m[0]: #filter only add nodes from module being targeted
            addnode = g.add_node(node,
                    db="wgcna",
                    color=modulecolor,
                    athbesthit=bh_ath[node][0],
                    athbesthitannot=bh_ath[node][1],
                    orfid = compggfunc[node][1],
                    Orthogroup = compggfunc[node][2],
                    MAPKs = compggfunc[node][3],
                    TFs = compggfunc[node][4],
                    RGA = compggfunc[node][5],
                    IPR = compggfunc[node][6],
                    Pfam = compggfunc[node][7],
                    TMHMM2 = compggfunc[node][8],
                    TargetP2 = compggfunc[node][9],
                    SignalP5 = compggfunc[node][10],
                    SP_UPDOWN = compggfunc[node][14],
                    SP_LFC = compggfunc[node][15],
                    IAC_UPDOWN = compggfunc[node][16],
                    IAC_LFC = compggfunc[node][17],
                    RB05_UPDOWN = compggfunc[node][18],
                    RB05_LFC = compggfunc[node][19],
                    RB200_UPDOWN = compggfunc[node][20],
                    RB200_LFC = compggfunc[node][21],
                    Blast2GO = compggfunc[node][22],
                    GO = compggfunc[node][23],
                    KOG = compggfunc[node][24])

    print ("Parsing edges from file:",targetedgefile)
    with open(targetedgefile,"r") as set1:
        for i in set1:
            i = i.rstrip().split("\t")
            if "fromNode" not in i[0]:
                fromnode = i[0]
                tonode = i[1]
                weight = float(i[2])
                if weight >= cutoff_weight:
                    addedge = g.add_edge(fromnode,tonode,
                                    source="wgcna co-expression",
                                    modulecolor=targetmodulecolor,
                                    weight=weight)


    print ("Saving targeted module in",diroutput)
    os.chdir(diroutput)
    export_graphml = nx.write_graphml(g, "WGCNA_targeted_module_"+targetmodulecolor+".xml")
