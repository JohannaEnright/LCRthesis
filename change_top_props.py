#!/usr/bin/python3

import subprocess
import tempfile
from re import search
from sys import argv as a

def findindexedletters(tfName, start, end):

    codstochange = []
    with open(tfName, "r") as fh:
        f = fh.readlines()
        for num in range(start, end+1):
            #print(f[num])
            cod = search(r"^(\w)*\S", f[num])
            codstochange.append(cod.group())
    return codstochange


def alterproportion(prop, percent):

    percent = percent/100  #convert percent to decimal
    newprop = prop*(percent)
    changedby = prop - newprop    #will be positive if the newproportions were decreased and negative if it was increased
    #print([round(newprop,5), round(changedby,5)])
    return [round(newprop,5), round(changedby,5)]


#will take a nt_aa_cod file and will change the letters from start, end, by the given proportion
def changegeneralprops(fileR, start, end, percent):

    start = int(start); end = int(end); percent = float(percent)/100  #convert percent to decimal
    if end < start:
        raise ValueError("end must be greater than start")

    dcod = {"ttt":1, "ttc":2, "tta":3, "ttg":4, "tct":5, "tcc":6, "tca":7, "tcg":8, "tat":9, "tac":10,"tgt":11,\
             "tgc":12, "tgg":13, "ctt":14, "ctc":15, "cta":16, "ctg":17, "cct":18, "ccc":19, "cca":20, "ccg":21, "cat":22,\
             "cac":23, "caa":24, "cag":25, "cgt":26, "cgc":27, "cga":28, "cgg":29, "att":30, "atc":31, "ata":32, "atg":33,\
             "act":34, "acc":35, "aca":36, "acg":37, "aat":38, "aac":39, "aaa":40, "aag":41, "agt":42, "agc":43, "aga":44,\
             "agg":45, "gtt":46, "gtc":47, "gta":48, "gtg":49, "gct":50, "gcc":51, "gca":52, "gcg":53, "gat":54, "gac":55,\
             "gaa":56, "gag":57, "ggt":58, "ggc":59, "gga":60, "ggg":61}    #use this to index at the proper location for the corresponding proportion
    
    daa = {"a":1, "c":2, "d":3, "e":4, "f":5, "g":6, "h":7, "i":8, "k":9,"l":10,\
            "m":11, "n":12, "p":13, "q":14, "r":15, "s":16, "t":17, "v":18, "w":19, "y":20}  #this list is in the same order as the proportions from the nt_aa_codon bias file
    
    dnt = {"a":1, "c":2, "g":3, "t":4}


    with open(fileR, "r") as fh:
        f = fh.readlines()
        h1 = f[0]
        h2 = f[1][:-1]+". Frequencies altered by "+str(percent*100)+"%: "+"letters from "+str(start)+" to "+str(end)+" when arranged most used to least used.\n"
        letters = f[2]   #line containing the letters
        #print(let)
        freq = f[3].split("\t")[:-1]   #remove "\n" at the end
        props = f[4].split("\t")[:-1]  #remove weird "" at the end

        #print(freq)
        #print(len(freq))

    if len(freq) == 62:   #codonfile
        if start  > 61 or end > 61:
            raise ValueError("must input start and end values less than 62")
        ordered = subprocess.check_output(["ordercodons", fileR, "codon"])
        d = dcod  #for indexing the dictionary
    elif len(freq) == 21:  #aafile
        if start > 20 or end > 20:
            raise ValueError("must input start and end values less than 21")
        ordered = subprocess.check_output(["ordercodons", fileR, "aa"])
        d = daa  #for indexing the dictionary
    elif len(freq) == 5:  #ntfile
        if start > 4 or end > 4:
            raise ValueError("must input start and end values less than 5")
        ordered = subprocess.check_output(["ordercodons", fileR, "nt"])
        d = dnt  #for indexing the dictionary
    tfName = tempfile.NamedTemporaryFile(delete = False).name
    with open(tfName, "w") as fh:
        fh.write(ordered.decode())

    list_letterstochange = findindexedletters(tfName, start, end)
    subprocess.call(["rm", tfName])
    #print(list_letterstochange)


    for let in list_letterstochange:
        letfreq = float(freq[d[let]])
        #print(d[let])
        newfreq = round(letfreq*percent, 5)
        #print(newfreq)
        freq[d[let]] = newfreq

    new_total = 0
    for f in freq[1:]:   #remove word 'float'
        #print(f)
        new_total += float(f)
    newprops = ["prop"]
    for f in freq[1:]: 
        newprop = float(f)/new_total
        newprops.append(str(round(newprop,5)))
    #print(newprops)

    #now convert everything back to a string to write to a file
    freq_line = ""
    for i in freq:
        freq_line += str(i)+"\t"
    freq_line += "\n"
    prop_line = ""
    for i in newprops:
        prop_line += str(i)+"\t"

    to_file = h1+h2+letters+freq_line+prop_line
    if "/" in fileR:
        i = fileR.rindex("/")
        fileW = fileR[i+1:]+"alt_"+str(start)+"_"+str(end)+"_"+str(int(percent*100))
    else:
        fileW = fileR+"alt_"+str(start)+"_"+str(end)+"_"+str(int(percent*100))

    #print(to_file)
    with open(fileW, "w") as fh:
        fh.write(to_file)
    #print("\n\n\n")

#changegeneralprops("/home/JosLaptop/proj1/nt_cod_aa_proportions/LCRcounts/ScereLCRcodon_count", "1", "5", "200")
#changegeneralprops("/home/JosLaptop/proj1/nt_cod_aa_proportions/LCRcounts/ScereLCRnt_count", "1", "4", "80")
#changegeneralprops("/home/JosLaptop/proj1/nt_cod_aa_proportions/LCRcounts/ScereLCRaa_count", "1", "5", "80")
try:
    changegeneralprops(a[1], a[2], a[3], a[4])
except:
    print("IndexError: arguments are filename (nt_cod_aa format, start (int), end (int), and percent)")
