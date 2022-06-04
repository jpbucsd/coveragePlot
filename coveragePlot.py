import sys
import os
import matplotlib.pyplot as plt
import numpy as np

td = False #program uses a tag directory location to plot
tdf = False
tdname = '' #tag directory directory

nf = False #the program uses a file comparing tag directory name to the names. otherwise it will name the files by tag directory
nff = False
name = ''

op = False #the program has an output directory to use, otherwise default
opf = False
out = '' #output directory

sf = False
ef = False
start = 0;
end = 10000;

lf = False

cont = True

helpS = """Options\n\n
\t-t <tag directory>    - tag directory produced by HOMER containing bedgraph files, required\n
\t-d <.txt>             - txt file with each line containing the tag directory name and the name for the plot produced, separated by a tab\n
\t-o <output directory> - output directory for images\n
\t-s <integer>          - start index in genome to begin graph\n
\t-t <integer>          - end index in genome to begin graph\n
\t-l                    - an option to include the start and stop location in file name
"""

for arg in sys.argv:
    if(tdf):
        td = True
        tdname = arg
        tdf = False
    elif(nff):
        nf = True
        name = arg
        nff = False
    elif(opf):
        op = True
        out = arg
        opf = False
    elif(sf):
        start = int(arg)
        sf = False
    elif(ef):
        end = int(arg)
        ef = False
    elif (arg == "-t") or (arg == "-T") or (arg == "--tag"):
        tdf = True
    elif (arg == "-d") or (arg == "-D") or (arg == "--dict"):
        nff = True
    elif (arg == "-o") or (arg == "-O") or (arg == "--out"):
        opf = True
    elif (arg == "-h") or (arg == "-H") or (arg == "--help"):
        cont = False
        print(helpS)
    elif (arg == "-s") or (arg == "-S") or (arg == "--Start"):
        sf = True
    elif (arg == "-l") or (arg == "-L") or (arg == "--Location"):
        lf = True
    elif (arg == "-e") or (arg == "-E") or (arg == "--End"):
        ef = True
    
#print("names: " + name)
#print("tag dir: " + tdname)
#print("output dir: " + out)
#print("start loc: " + str(start))
#print("end loc: " + str(end))
#SRR3401655 -> Strand specific RNA-seq on WT G27 strain nickel treated, replica B
#SRR3401654 -> Strand specific RNA-seq on WT G27 strain nickel treated, replica A
#SRR3401651 -> Strand specific RNA-seq on WT G27 strain untreated, replica B
#SRR3401650 -> Strand specific RNA-seq on WT G27 strain untreated, replica A

# SRR3401647	 -> Strand specific RNA-seq on deltaNikR strain nickel treated, replica B
#12. SRR3401643 -> Strand specific RNA-seq on deltaNikR strain nickel treated, replica A
#13. SRR3401641 -> Strand specific RNA-seq on deltaNikR strain untreated, replica B
#14. SRR3401620 -> Strand specific RNA-seq on deltaNikR strain untreated, replica A
if end < start:
    cont = False
    print("end must be greater than or equal to start")
accessions = []
names = {}
if(cont): 
    #create names dict
    
    nl = []
    
    if(nf):
        nfile = open(name,"r")
        for line in nfile.readlines():
            n = False
            w = ''
            key = ''
            for char in line:
                if(n):
                    if char != '\n':
                        w += char
                elif char == '\t':
                    n = True
                    key = w
                    w = ''
                else:
                    w += char
            names[key] = w
            nl.append(key)
    if(td):
        for subdir in os.listdir(tdname):
            accessions.append(str(subdir))
    if(op):
        found = False
        for subdir in os.listdir():
            if(subdir == out):
                found = True
        if(not found):
            os.system("mkdir " + out)
            
    if(nf and td):
        for acc in accessions:
            found = False
            for nn in nl:
                if acc == nn:
                    found = True
            if(found == False):
                print("tag directory " + acc + " missing in dictionary file")
            cont = found
        
           
if(cont):
    total = end - start + 1;
    
    #names = {
    #    "SRR3401655": "Strand specific RNA-seq on WT G27 strain nickel treated, replica B",
    #    "SRR3401654": "Strand specific RNA-seq on WT G27 strain nickel treated, replica A",
    #    "SRR3401651": "Strand specific RNA-seq on WT G27 strain untreated, replica B",
    #    "SRR3401650": "Strand specific RNA-seq on WT G27 strain untreated, replica A",
    #    "SRR3401647": "Strand specific RNA-seq on deltaNikR strain nickel treated, replica B",
    #    "SRR3401643": "Strand specific RNA-seq on deltaNikR strain nickel treated, replica A",
    #    "SRR3401641": "Strand specific RNA-seq on deltaNikR strain untreated, replica B",
    #    "SRR3401620": "Strand specific RNA-seq on deltaNikR strain untreated, replica A"
    #}

    #accessions = ["SRR3401655","SRR3401654","SRR3401651","SRR3401650","SRR3401647","SRR3401643","SRR3401641","SRR3401620"]
    pArray = {
        "loci" : [0 for i in range(total)]
    }
    nArray = {
        "loci" : [0 for i in range(total)]
    }   
    for i in range(total):
        pArray["loci"][i] = i + start
        nArray["loci"][i] = i + start


    SR = 0
    for SRA in accessions:

        #initialize storage
        pArray[SRA] = [0 for i in range(total)]
        nArray[SRA] = [0 for i in range(total)]
        #prep files
        os.system("gzip -d tagDirsRNA/" + SRA + "/" + SRA + ".ucsc.bedGraph.gz")
        file = open("tagDirsRNA/" + SRA + "/" + SRA + ".ucsc.bedGraph","r")
        ln = 0
        base = 0
        pos = True
        for line in file.readlines(): 
            if ln > 0:
                place = 0
                sN = 0 #base start for line
                eN = 0 #base end for line
                score = 0
                col = ''
                for char in line:
                    if (char == '\t'):
                        if(place == 1):
                            sN = int(col)
                        elif(place == 2):
                            eN = int(col)
                        elif(place == 3):
                            score = float(col)
                        place += 1
                        col = ''
                    else:
                        col += char
                score = float(col)
                if(eN < base):
                    #we are now on the negative strand
                    base = 0
                    pos = False
                if(base >= start and base < end - 1):
                    for i in range(eN - sN + 1):
                        if(base - start + i >= total):
                            doNothing = 0
                        else:
                            if(pos):
                                #print(base - start + i)
                                pArray[SRA][base - start + i] = score
                            else:
                                nArray[SRA][base - start + i] = score
                base = base + (eN - sN)
            ln += 1

        SR += 1
        #positive strand
        fig = plt.figure()
        plot = fig.add_subplot(111)
        plot.plot(pArray["loci"],pArray[SRA],linewidth=1.0,color = "blue")
        plot.fill_between(pArray["loci"],pArray[SRA],alpha=0.2, color='blue')
        plot.set_xlabel("loci")
        plot.set_ylabel("Coverage")
        title = ''
        if(nf):
            title += names[SRA]
        else:
            title += SRA
        title += "(+strand)"
        plot.set_title(title)
        if(op):
            title = out + "/" + title
        if(lf):
            title += "_" + str(start) + "-" + str(end)
        title += ".png"
        plt.savefig(title)
        print("created plot: " + title)
        #negative strand
        fig2 = plt.figure()
        plot2 = fig2.add_subplot(111)
        plot2.plot(nArray["loci"],pArray[SRA],linewidth=1.0,color = "red")
        plot2.fill_between(pArray["loci"],pArray[SRA],alpha=0.2, color='red')
        plot2.set_xlabel("loci")
        plot2.set_ylabel("Coverage")
        title = ''
        if(nf):
            title += names[SRA]
        else:
            title += SRA
        title += "(-strand)"
        plot2.set_title(title)
        if(op):
            title = out + "/" + title
        if(lf):
            title += "_" + str(start) + "-" + str(end)
        title += ".png"
        plt.savefig(title)
        print("created plot: " + title)