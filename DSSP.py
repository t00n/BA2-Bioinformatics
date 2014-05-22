import os
import Score

ROOT_DIR = "dataset/"
DSSP_DIR = ROOT_DIR + "dssp/"

class DSSPChain:
    def __init__(self):
        self.id = ""
        self.protein = ""
        self.organism = ""
        self.sequence = ""
        self.structure = ""
        
    def __repr__(self):
        ret = ">"
        ret += self.id + "|" + self.protein + "|" + self.organism + "\n"
        ret += self.sequence + "\n"
        ret += self.structure
        return ret

class DSSPData(list):
    def __init__(self, *args):
        list.__init__(self, *args)

    def loadFromFile(self, filename, chains):
        file = open(DSSP_DIR + filename + ".dssp", 'r')
        start=False
        data = dict.fromkeys(chains, DSSPChain())
        for line in file:
            if(line.find('#') != -1):
                start=True
            elif(not start):
                if (line[0:6] == "COMPND"):
                    for c in data:
                        data[c].protein = line[line.find(": ")+2:line.find(";")]
                elif (line[0:6] == "SOURCE"):
                    for c in data:
                        data[c].organism = line[line.find(": ")+2:line.find(";")]
            elif(start):
                chain = line[11:12].upper()
                aa = line[13:14].upper()
                if (chain in data and aa not in ["B", "Z", "X", "J"]):
                    struct = line[16:17].upper()
                    data[chain].sequence += aa
                    if (struct in ["H", "G", "I"]):
                        data[chain].structure += "H"
                    elif (struct in ["C", "S", "B", " "]):
                        data[chain].structure += "C"
                    else:
                        data[chain].structure += struct
        for c in data:
            data[c].id = filename+c
            self.append(data[c])
        file.close()

class InfoData:
    def load(self, filename):
        CATH_info = open(filename)
        files = dict()
        for line in CATH_info:
            filename = line[:4]
            seqnum = line[4:5]
            if (filename not in files):
                files[filename] = []
            if (seqnum not in files[filename]):
                files[filename].append(seqnum)
        CATH_info.close()
        return files

if __name__ == '__main__':
    info = InfoData()
    files = info.load(ROOT_DIR + "CATH_info.txt")
    dssp = DSSPData()
    for filename in sorted(files.keys()):
        dssp.loadFromFile(filename, files[filename])
    cpt = 0
    for line in dssp:
        print(line)
        cpt += len(line.sequence)

    print(len(dssp))
    print(cpt)