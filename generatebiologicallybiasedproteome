#!/usr/bin/python3

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
from sys import argv

#takes probability of selecting a,g,c,t and strings it together as a list of weights []
#makes 100 000, 100 codon long coding sequences (need 303 nts)
#seed defaults to 10
#num_iterations defaults to 100 000
#organism name defaults to 'test'
def biologicalbiasnts(weight, seed, num_iterations, organism_name):
    print("nt bias")

    random.seed(seed)
    nts = ["a", "c", "g", "t"]
    
    h = 0
    allcds = []
    while h < num_iterations:
        i = 0
        coding_sequence = ""
        while i < 100:
            codon=""            #add 3 nts to make sure the codon is not a stop codon
            j = 0
            while j < 3:
                nt = random.choices(nts, weights=weight, k=1)[0]
                codon += nt
                j+=1
            if codon == "taa" or codon == "tag" or codon == "tga":  #keep chossing nts until you don't get stop codon
                while codon == "taa" or codon == "tag" or codon == "tga":
                    codon=""  #reset to blank list
                    j = 0
                    while j < 3:
                        nt = random.choices(nts, weights=weight, k=1)[0]
                        codon += nt
                        j+=1
            coding_sequence += codon
            i += 1
        coding_sequence += "taa"
        allcds.append(coding_sequence)
        print(coding_sequence)
        print(len(coding_sequence))
        print("")
        h += 1
    to_file_DNA = ""
    ID = 0
    for sequence in allcds:
        to_file_DNA += ">" + str(ID) + "\n"\
                + sequence + "\n\n"
        ID += 1

    fastafilename = organism_name + "_nt_s" + str(seed) + "_" + str(weight[0]) + "_" + str(weight[1]) +\
            "_" + str(weight[2]) + "_" + str(weight[3]) + ".fasta"
    with open(fastafilename, "w") as fh:
        fh.write(to_file_DNA)
    return fastafilename


def biologicalbiascodons(weight, seed, num_iterations, organism_name):
    print("codon bias")
    codons = ["ttt", "ttc", "tta", "ttg", "tct", "tcc", "tca", "tcg", "tat", "tac", "tgt",\
              "tgc", "tgg", "ctt", "ctc", "cta", "ctg", "cct", "ccc", "cca", "ccg", "cat",\
              "cac", "caa", "cag", "cgt", "cgc", "cga", "cgg", "att", "atc", "ata", "atg",\
              "act", "acc", "aca", "acg", "aat", "aac", "aaa", "aag", "agt", "agc", "aga",\
              "agg", "gtt", "gtc", "gta", "gtg", "gct", "gcc", "gca", "gcg", "gat", "gac",\
              "gaa", "gag", "ggt", "ggc", "gga", "ggg"]   #this list is in the same order as the proportions from the nt_aa_codon bias file
    random.seed(seed)

    allcds = []
    h=0
    while h < num_iterations:  #generate 100 000 sequences with 101 codons each
        coding_sequence = ""
        for i in range(100):
            codon = random.choices(codons, weights=weight, k=1)[0]
            coding_sequence += codon
        coding_sequence += "taa"
      #  print(coding_sequence)
      #  print(len(coding_sequence))
        allcds.append(coding_sequence)
        h+=1

    ID=0
    to_file_DNA = ""
    for sequence in allcds:
        to_file_DNA += ">" + str(ID) + "\n"\
                + sequence + "\n\n"
        ID += 1
    fastafilename = organism_name + "_codon_s" + str(seed) + ".fasta"
    with open(fastafilename, "w") as fh:
        fh.write(to_file_DNA)
    return fastafilename


def biologicalbiasaas(weight, seed, num_iterations, organism_name): 
    print("aa bias")
    aas = ["a", "c", "d", "e", "f", "g", "h", "i", "k","l",\
           "m", "n", "p", "q", "r", "s", "t", "v", "w", "y"]  #this list is in the same order as the proportions from the nt_aa_codon bias file

    aatocod = {"a": ["gct", "gcc", "gca", "gcg"], "c":["tgt", "tgc"], "d":["gat", "gac"],
               "e":["gaa", "gag"], "f":["ttt", "ttc"], "g": ["ggt", "ggc", "gga", "ggg"], 
               "h":["cat", "cac"], "i":["att", "atc", "ata"], "k":["aaa", "aag"],"l":["tta", "ttg", "ctt", "ctc", "cta", "ctg"],
               "m":["atg"], "n":["aat", "aac"], "p":["cct", "ccc", "cca", "ccg"], "q":["caa", "cag"],
               "r":["aga", "agg", "cgt", "cgc", "cga", "cgg"], "s":["tct", "tcc", "tca", "tcg", "agt", "agc"],
               "t":["act", "acc", "aca", "acg"],"v":["gtt", "gtc", "gta", "gtg"], "w":["tgg"], "y":["tat", "tac"]}

    random.seed(seed)
    allcds = []
    h=0
    while h < num_iterations:
        coding_sequence = ""
        for i in range(100):  #make aa sequence according to proportions
            aa = random.choices(aas, weights=weight, k=1)[0]
            randomcod = random.choices(aatocod[aa], k=1)[0]  #pick a codon for that aa
            coding_sequence += randomcod  #choose a codon for this aa (with equal probability)
        coding_sequence += "taa"
        allcds.append(coding_sequence)
        h+=1

    #now make the fasta file
    ID=0
    to_file_DNA = ""
    for sequence in allcds:
        to_file_DNA += ">" + str(ID) + "\n"\
                + sequence + "\n\n"
        ID += 1
    fastafilename = organism_name + "_cdsaabias_s" + str(seed) + ".fasta"
    with open(fastafilename, "w") as fh:
        fh.write(to_file_DNA)
    return fastafilename



#takes a fasta file with DNA sequences and converts it to a fasta file
#translated protein sequences with the same ID
def MakeTranslatedFasta(CDSfasta):
    to_file_protein = ""
    for record in SeqIO.parse(CDSfasta, "fasta"):
        ID = record.id
        seq = record.seq
        to_file_protein += ">"+ID+"\n"+ str(Seq(seq).translate()[:-1])+"\n\n"
    filename = CDSfasta[:-6] + ".translate.fasta"
    with open(filename, "w") as fh:
        fh.write(to_file_protein)
    print("Fasta file is written")

##fasta = biologicalbiasnts(0.1, 0.1, 0.4, 0.4, 10, 100000, "test")  #for indv testing
##MakeTranslatedFasta(fasta)

#takes a fasta file of coding sequences and puts those sequnces into a
#genbankfile, with corresponding protein sequences
def MakeGenbankFile(CDSfasta):
    print("making genbank file")
  
    #concatenate all DNA seqs and make it the chromosome sequence
    chromosome_string = ""
    for rec in SeqIO.parse(CDSfasta, "fasta"):
        chromosome_string += rec.seq
    chromosome_object = Seq(chromosome_string)

    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

    name = "biologically_proportional_simulation"
    description = "biologically proportional genbank file containing 100 codon long\
            DNA sequences or 100 aa long protein sequences with a designated\
            proportion of letters"
    organism = "biologically proportional"
    accessions = "NC_"+CDSfasta
    genbankname = CDSfasta[:-6]+".gb" 
    print("genbankname", genbankname)

    record = SeqRecord(chromosome_object, id= dt_string, name=name,\
            description= description, annotations={"molecule_type":"DNA", "accessions":\
            [accessions], "organism": organism})

    for rec in SeqIO.parse(CDSfasta, "fasta"):
        print("adding features")
        ID = rec.id
        DNAseq = rec.seq
        start = int(ID)*303
        end = start + 303
        translation = Seq(DNAseq).translate()[:-1]
       
        #create feature annotations
        feature = SeqFeature(FeatureLocation(start=start, end=end), type = "CDS", strand=1,\
                  qualifiers= {"protein_id":[ID], "translation":[translation], "db_xref":["GeneID:"+ID]})
        record.features.append(feature)

    #save as genbank file
    with open(genbankname, "w") as fh:
        SeqIO.write(record, fh, "genbank")
    print("genbank file is written")

#takes a file produced from 'nt_pro_cod_genomebias' script (becuase it has an ugly output that needs to be parsed 
# into a nice list of proportions)
#returns a list of proportions (floats) (given in the same order as the file) 
def getproplist(filename):

    with open(filename) as fh:
        file = fh.readlines()
        for line in file:
            if "prop" in line:
                list_props = line.split("\t")
                weights = []  #this will have the nt or codon proportions appended to it
                for i in range(1,len(list_props)-1):  #avoid 'props' and '' at start and end of list_props
                    weights.append(float(list_props[i]))
    print(weights)
    return weights

#"proj1/nt_cod_aa_proportions/Homo_sapiens_codon_count"

#can specify ntcodpropfile OR weight, not both
def executefunctions(ntcodpropfile=None, weight=None, seed=10, num_iterations=100000, organismname="test"):

    print("started running")
    if weight == None: #if a proportion file is specified and not a list of weights, then parse proportion file
        weight = getproplist(ntcodpropfile)

    #else:
    #    weight = weight.split(" ") #make string of propotions into a list
    #    for i in range(len(weight)):  #convert all nums to float
    #        weight[i] = float(weight[i])

    #now call on the proper function: are you generating proteome based off nt, codon, or aa proportions
    if len(weight) == 4:
        CDSfasta = biologicalbiasnts(weight, seed, num_iterations, organismname) #generates fasta based on given nucleotide proportions
    elif len(weight) == 20:
        CDSfasta = biologicalbiasaas(weight, seed, num_iterations, organismname)
    elif len(weight) == 61:
        CDSfasta = biologicalbiascodons(weight, seed, num_iterations, organismname)
    else:
        raise ValueError("Did not enter correct number of letter proportions must be 4, 20, or 61")

    #weight argument is a list of floats (nt or aa or codon proportions)
    MakeTranslatedFasta(CDSfasta)  #a fasta of the proteins just in case you want to quickly run it through seg
    MakeGenbankFile(CDSfasta)      #the genbankfile to run through CombinedUsingSysCall program


##COMMAND LINE STUFF##
args = argv[1:]
executefunctions(ntcodpropfile=args[0], seed=10, num_iterations=100000,organismname=args[1] )
#if args[0] == "-h" or args == "--help":
#    raise TypeError("put '-f' or '-s' (file=-f, string of weights=-s), 'nt_pro_cod_genomebias_filename'\
#            OR string of weights 'a c g t' separated by a space, seed(d=10), number_of_sequences_generated(d=100 000), organism_name(d=test)(for file name)")
#elif len(args) < 2:
#    raise TypeError("must specify '-f'(filename) OR '-s'(string of weights), filename containing proportions OR string of nt/cod/aa weights in order")
#if args[0] != "-f" and args[0] != "-s":
#    raise ValueError("first argument must be a flag: '-f' or '-s'")
#if len(args) > 2:
#    seed = int(args[2])
#else:
#    seed = 10  #default
#if len(args) > 3:
#    num_iterations = int(args[3])
#else:
#    num_iterations = 100000  #default
#if len(args) > 4:
#    organismname = args[4]
#else:
#    organismname = "test"  #default
#
#if args[0]=="-f": #then a filename was provided
#    executefunctions(ntcodpropfile=args[1], seed=seed, num_iterations=num_iterations, organismname=organismname)
#elif args[0]=="-s": #then a string of weights was provided
#    executefunctions(weight=args[1], seed=seed, num_iterations=num_iterations, organismname=organismname)

#args[0] will eighter be file of nt/aa/codon proportions (specifically from 'nt_pro_cod_genomebias') output
#or it will be a string of weights entered manually on command line "0.2 0.3 0.2 0.4" in between "" and separated 
#by a space

