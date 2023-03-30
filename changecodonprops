#!/usr/bin/python3

from re import search
from sys import argv as a
import subprocess
import tempfile

#this will create a temp file of ordered codons from most used to least used (allows you to determine what proportions to change later)
def OrderCodons(fileR, Type):
    to_file = "#ordered letter propotions from "+fileR+"\n"
    with open(fileR, "r") as fh:
        f = fh.readlines()

        if Type == "codon":
            codons = f[2].split("\t")
            prop = f[4].split("\t")
            new_codons = []
            for c in codons:
                cod = search(r"\w\w\w", c)  #just get rid of the wierd spaces after the codons and such
                if cod:
                    new_codons.append(cod.group())
            new_prop = []
            for p in prop:
                pr = search(r"\S\..*\S", p)  #get rid of the weird spaces after the digits
                if pr:
                    new_prop.append(pr.group())
            sorted_new_prop = sorted(new_prop, reverse=True)  #sort codon proportions from highest to lowst
            d = {}
            for cod, pro in zip(new_codons, new_prop):
                d[cod] = pro
            sorted_new_codons = []
            for p in sorted_new_prop:  #now put the codons in order from most frequently used to least frequently used
                for k,v in d.items():
                    if p == v:
                        sorted_new_codons.append(k)
                        d[k]= None  #change it to 0, so codons aren't used more than once if values are the same

            ordered_cod_prop = zip(sorted_new_codons, sorted_new_prop)

            for tup in ordered_cod_prop:
                to_file += tup[0] + "\t" + tup[1] + "\n"

    tfName1 = tempfile.NamedTemporaryFile(delete = False).name
    with open(tfName1, "w") as fh:
        fh.write(to_file)
    return tfName1

#will return the top n codons of type codtype whose proportions you should change (as a list)
def findtopcodons(tfName1, codtype, n):

    count = 0
    codonprops_tochange = []
    with open (tfName1, "r") as fh:
        f = fh.readlines()
        for line in f:
            if count < int(n):
                if "#" not in line:
                    cod_prop = line[:-1].split("\t")
                    a = cod_prop[0][0]; b = cod_prop[0][1]; c = cod_prop[0][2]

                    if codtype == "1":  #class1 (all the same nts)
                        if  a == b and a == c:  #indexes the individual letters within the codon
                            codonprops_tochange.append(cod_prop[0])
                            count += 1
                    elif codtype == "2":  #class2
                        if (a == b and a != c) or (a == c and a != b) or (b == c and a != b):
                            codonprops_tochange.append(cod_prop[0])
                            count += 1
                    elif codtype == "3":  #class3
                        if a != b and a != c and b != c:
                            codonprops_tochange.append(cod_prop[0])
                            count += 1
    #print(codonprops_tochange)
    return codonprops_tochange

#takes propotion as float type and percent to increase or decrease the proportion by that amount.
#exp 90% (will change the proportion to 90% of what it was (decrease)
#120 will change the proportions to 120% of what they were
def alterproportion(prop, percent):

    percent = percent/100  #convert percent to decimal
    newprop = prop*(percent)
    changedby = prop - newprop    #will be positive if the newproportions were decreased and negative if it was increased
    return [round(newprop,5), round(changedby,5)]
    #return [newprop, changedby]

#fileR form nt_cod_aa_proportions, codon type (1, 2, 3), n as in how many (from the highest to lowest to change),
# percent to alter the codon proportion by, and if you are increasing the proportion or decreasing the proportion
def changecodonprops(fileR, codtype, n, percent):

    d = {"ttt":1, "ttc":2, "tta":3, "ttg":4, "tct":5, "tcc":6, "tca":7, "tcg":8, "tat":9, "tac":10,"tgt":11,\
         "tgc":12, "tgg":13, "ctt":14, "ctc":15, "cta":16, "ctg":17, "cct":18, "ccc":19, "cca":20, "ccg":21, "cat":22,\
         "cac":23, "caa":24, "cag":25, "cgt":26, "cgc":27, "cga":28, "cgg":29, "att":30, "atc":31, "ata":32, "atg":33,\
         "act":34, "acc":35, "aca":36, "acg":37, "aat":38, "aac":39, "aaa":40, "aag":41, "agt":42, "agc":43, "aga":44,\
         "agg":45, "gtt":46, "gtc":47, "gta":48, "gtg":49, "gct":50, "gcc":51, "gca":52, "gcg":53, "gat":54, "gac":55,\
         "gaa":56, "gag":57, "ggt":58, "ggc":59, "gga":60, "ggg":61}    #use this to index at the proper location for the corresponding proportion

    with open(fileR, "r") as fh:
        f = fh.readlines()
        h1 = f[0]
        h2 = f[1][:-1]+". props altered: class" + codtype + ", num changed: "+n+"\n"
        cod = f[2]
        freq = f[3]
        props = f[4].split("\t")
      #  print(props)


    if len(props) == 63:  #order the file according to most used/ least used codon
        tfName1 = OrderCodons(fileR, "codon")
        list_codstochange = findtopcodons(tfName1, codtype, n)
        subprocess.call(["rm", tfName1])   #delete tempfile1

        totalprop_changedby = 0
        for c in list_codstochange:
            codonprop = props[d[c]]  #this gives the proportion for that codon (can now change it... maybe decrease it by 10%)
            newprop_info = alterproportion(float(codonprop), float(percent))   #result is a list containing the [newprop, amount it has been changed by]
            props[d[c]] = newprop_info[0]  #changes the proportion in the list called 'props'
            totalprop_changedby += newprop_info[1]  #save this to decrease the other proportions comparatively so total equals 1

    #now we need to change the other proportions accordingly 
    addby = totalprop_changedby / (61-len(list_codstochange))  #will be negative if percent is over 100 (then the other proportions are subtracted from)
    #print(addby)
    for k, v in d.items():
        if k not in list_codstochange:
            props[v] = round(float(props[v]) + addby, 5)  #add other proportions by the corresponding fraction
    #now convert it back into a string to write to the file with the newly changed proportions
    prop_line = ""
    for i in props[:-1]:  #get rid of the weird ""
        prop_line += str(i) + "\t"


   # is1 = 0    #   THIS FOR LOOP IS JUST FOR CHECKING
   # for i in range(1, len(props)-1):  #ignore "prop" at start and " " at end for now
   #     is1 += float(props[i])
   # print(is1)

    #write new info to file
    to_file = h1 + h2 + cod + freq + prop_line
    #print(to_file)
    if "/" in fileR:
        i = fileR.rindex("/")
        fileW = fileR[i+1:] + "alt_"+codtype+"_"+n+"_"+percent
    else:
        fileW = fileR + "alt_"+codtype+"_"+n+"_"+percent
    with open(fileW, "w") as fh:
        fh.write(to_file)
   # print("file is written")

#changecodonprops("proj1/nt_cod_aa_proportions/LCRcounts/alteredprops/ScereLCRcodon_count_changedprops", "3", "2", "50")
try:
    changecodonprops(a[1], a[2], a[3], a[4])
except:
    print("IndexError: arguments are 1:filename(in nt_cod_aa_bias format), 2:codon type (1, 2, or 3), 3:number of codons (from top) that you want to change), 4: percent to change them by (exp: '90' will result in 10% decrease for those codons)")





