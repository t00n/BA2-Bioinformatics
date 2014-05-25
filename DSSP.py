import os

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
    def loadMany(self, filename, dssp_dir):
        info = open(filename)
        files = dict()
        for line in info:
            filename = line[:4]
            seqnum = line[4:5]
            if (filename not in files):
                files[filename] = ""
            if (seqnum not in files[filename]):
                files[filename] += seqnum
        info.close()
        for filename in files:
            self.loadFromFile(dssp_dir + filename + ".dssp", files[filename])

    def loadFromFile(self, filename, chains):
        file = open(filename, 'r')
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
                chain = line[11:12]
                aa = line[13:14].upper()
                if (chain in data and aa not in "BZXJOU"):
                    struct = line[16:17].upper()
                    data[chain].sequence += aa
                    if (struct in "HGI"):
                        data[chain].structure += "H"
                    elif (struct in "CSB "):
                        data[chain].structure += "C"
                    elif (struct in "E"):
                        data[chain].structure += "E"
                    elif (struct in "T"):
                        data[chain].structure += "T"
                    else:
                        print("Structure Not Found : ", filename, line[0:5].strip())
        for c in data:
            data[c].id = os.path.basename(filename).split(".")[0]+c
            self.append(data[c])
        file.close()

    def saveResult(self, cache):
        with open(cache, "w") as output:
            for line in self:
                print(line, file=output)
            output.close()

if __name__ == '__main__':
    dssp = DSSPData()
    dssp.loadMany("dataset/CATH_info.txt", "dataset/dssp/")
    dssp.saveResult("dataset/summary.txt")
