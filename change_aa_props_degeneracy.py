#!/usr/bin/python3

import subprocess
import tempfile
import re
from sys import argv as a

def findindexedletters(tfName, degtype, start, end):

    aadeg = {"1":["m", "w"], "2":["n", "d", "c", "q", "e", "h", "k", "f", "y"],
             "3":["i"], "4":["a", "g", "p", "t", "v"], "6":["r", "l", "s"]}
    letstochange = []

    if start != None and end != None:
        pass
        #with open(tfName, "r") as fh:
        #    f = fh.readlines()
        #    for line in f:
        #        let = re.search(r"^([a-z])*\S", line)
        #        if let:
        #            pass
    else:
        for num in degtype:
            for aa in aadeg[num]:
                letstochange.append(aa)
    return letstochange


def alterproportion(prop, percent):

    percent = percent/100  #convert percent to decimal
    newprop = prop*(percent)
    changedby = prop - newprop    #will be positive if the newproportions were decreased and negative if it was increased
    return [round(newprop,5), round(changedby,5)]


#will take a nt_aa_cod file and will change the frequencies of amino acids of specified degeneracies by the given proportion
def changegeneralprops(fileR, aadegeneracy, percent, start=None, end=None):

    daa = {"a":1, "c":2, "d":3, "e":4, "f":5, "g":6, "h":7, "i":8, "k":9,"l":10,\
            "m":11, "n":12, "p":13, "q":14, "r":15, "s":16, "t":17, "v":18, "w":19, "y":20}  #this list is in the same order as the proportions from the nt_aa_codon bias file
    
    with open(fileR, "r") as fh:
        f = fh.readlines()
        h1 = f[0]
        h2 = f[1][:-1]+". Frequencies of amino acids with degeneracies "+aadegeneracy + " altered by "+percent+"%\n"
        letters = f[2]   #line containing the letters
        freq = f[3].split("\t")[:-1]   #remove "\n" at the end
        props = f[4].split("\t")[:-1]  #remove weird "" at the end

    aadegeneracy = re.split("\s", aadegeneracy)   #make a list of the amino acid degeneracies to change
    for num in aadegeneracy:
        if num != "1" and num != "2" and num != "3" and num != "4" and num != "6":
            raise ValueError("aa degeneracy can only be 1, 2, 3, 4, or 6")
    percent = float(percent)/100  #convert percent to decimal

    if len(freq) == 21:  #aafile
        ordered = subprocess.check_output(["ordercodons", fileR, "aa"])   #orderthem just in case in the future I only want to change certain ones
    else:  #ntfile
        raise ValueError("must input file of amino acid proportions")

    tfName = tempfile.NamedTemporaryFile(delete = False).name
    with open(tfName, "w") as fh:
        fh.write(ordered.decode())

    list_letterstochange = findindexedletters(tfName, aadegeneracy, start, end)
    subprocess.call(["rm", tfName])


    for let in list_letterstochange:
        letfreq = float(freq[daa[let]])
        newfreq = round(letfreq*percent, 5)
        freq[daa[let]] = newfreq

    new_total = 0
    for f in freq[1:]:   #remove word 'float'
        new_total += float(f)
    newprops = ["prop"]
    for f in freq[1:]: 
        newprop = float(f)/new_total
        newprops.append(str(round(newprop,5)))

    #now convert everything back to a string to write to a file
    freq_line = ""
    for i in freq:
        freq_line += str(i)+"\t"
    freq_line += "\n"
    prop_line = ""
    for i in newprops:
        prop_line += str(i)+"\t"

    aadeg_note = ""
    for n in aadegeneracy:
        aadeg_note += n+"_"

    to_file = h1+h2+letters+freq_line+prop_line
    if "/" in fileR:
        i = fileR.rindex("/")
        fileW = fileR[i+1:]+"altAA_"+aadeg_note+str(int(percent*100))
    else:
        fileW = fileR+"altAAdeg_"+aadeg_note+str(int(percent*100))

#    print(fileW)
#    print(to_file)
    with open(fileW, "w") as fh:
        fh.write(to_file)

#changegeneralprops("/home/JosLaptop/proj1/nt_cod_aa_proportions/LCRcounts/ScereLCRaa_count", "1 2 3", "80")
try:
    changegeneralprops(a[1], a[2], a[3])
except:
    print("IndexError: arguments are filename (nt_cod_aa format), 'string of aadeg (integers) separated by whitespace', and percent")
