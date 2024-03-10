######################33
#该模块的到的文件为fasta，下一步将其利用CD-hit软件去冗余

class Data_change():

    def __init__(self,file_path,data_name):
        # 文件路径
        self.file_path = r'data/'+file_path
        #将数据转换成fasta格式去冗余
        self.out_fasta_path = r'fasta/' + data_name

    def read_data(self):
        file = self.file_path
        with open(file,'r') as f:
            lines = f.readlines()
        NAME = []
        SEQ = []
        for line in lines:
            line = line.strip().split("\t")#0:CLMID 1:UniprotID 2:位点 3：修饰类型 4：基因 5：物种 6：序列 7&8：参考文献
            name = '>'+'\t'+line[1]+'\t'+line[5]+'\t'+line[2]+'\t'+line[3]#0:> 1:UniprotID 2:物种 3：位点 4：修饰类型
            seq = line[6]
            NAME.append(name)
            SEQ.append(seq)

        return NAME,SEQ


    def out_data(self):
        NAME, SEQ = self.read_data()
        Name_dic = {}

        for i in range(len(NAME)):
            name = NAME[i].strip().split("\t")
            ID = name[1]
            mod = name[4]
            site = int(name[3])
            if ID in Name_dic:
                if mod in Name_dic[ID]:
                    Name_dic[ID][mod].append(site)
                else:
                    Name_dic[ID][mod] = [site]
            else:
                mod_dic = {}
                mod_dic[mod] = [site]
                Name_dic[ID] = mod_dic

        new_NAME = []
        new_SEQ = []
        for i in range(len(NAME)):
            name = NAME[i].strip().split("\t")
            ID = name[1]
            new_ID = '>'+'\t'+ID
            if new_ID not in new_NAME:
                new_NAME.append(new_ID)
                new_SEQ.append(SEQ[i])

        return new_NAME,new_SEQ,Name_dic

    def out_fasta(self):
        NAME, SEQ ,Name_dic = self.out_data()
        out_fasta_path = self.out_fasta_path
        with open(out_fasta_path,'w') as f:
            for i in range(len(NAME)):
                print(NAME[i],file=f)
                print(SEQ[i],file=f)
        print("OK!")
        return Name_dic


if __name__ == '__main__':
    # Camellia sinensis.txt,Oryza sativa.txt,Triticum aestivum.txt
    data = Data_change(file_path='Camellia sinensis.txt', data_name = 'CS.txt')
    data.out_data()
    Name_dic = data.out_fasta()