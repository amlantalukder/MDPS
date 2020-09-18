import sys, os, argparse, getopt, ast

MATCH = 0
WOBBLE = 1
MISMATCH = 2
BULGE_MIRNA = 3
BULGE_MRNA = 4

NUMBER_STATES = 5

#-----------------------------------------------
def readFile(filename):

    fl = open(filename, "r")

    data = fl.readlines()

    fl.close()

    return data

#-----------------------------------------------
def writeFile(filename, data, mode="w"):

    d = os.path.dirname(filename)

    if d != "":
        if not os.path.exists(d):
            os.makedirs(d)

    fl = open(filename, mode)
    fl.write(data)
    fl.close()

#-----------------------------------------------
def readFileInTable(filename, delim='\t'):

    fl = open(filename, "r")

    data = fl.readlines()

    fl.close()

    data = [item.strip().split(delim) for item in data]

    return data

#-----------------------------------------------
def writeDataTableAsText(data, filename, mode="w"):

    text = formatDataTable(data, "\t", "\n")
        
    writeFile(filename, text, mode)

#-----------------------------------------------
def formatDataTable(data, col_sep="\t", row_sep="\n"):

    return row_sep.join([col_sep.join([str(item1) for item1 in item]) for item in data])

# ------------------------------------
def fastaToDict(file_name, header_del=''):

    fasta_data = readFile(file_name)
    fasta_data = "".join(fasta_data)
    fasta_data = fasta_data.split(">")[1:]
    fasta_data = [item.strip().split("\n") for item in fasta_data]

    if header_del != '':
        fasta_dict = dict([[item[0].split(header_del)[0], [item[0].split(header_del)[1:], "".join(item[1:])]] for item in fasta_data])
    else:
        fasta_dict = dict([[item[0], "".join(item[1:])] for item in fasta_data])
    
    return fasta_dict

# ------------------------------------
def showPercBar(counter, size, perc, perc_inc=10):
    num_prints = 100/perc_inc
    progress = int(counter*10/size)*10
    if progress >= perc:
        sys.stdout.write('='*(int((progress-perc)/num_prints)+1))
        sys.stdout.flush()
        if progress >= 100:
            print('100%')
        perc = progress + perc_inc
    return perc
