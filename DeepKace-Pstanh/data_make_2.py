#读取txt数据
#筛选修饰类型
#筛选蛋白ID
#按照蛋白ID统计修饰位点集合mod
#按照蛋白ID统计赖氨酸位置Q
#定义扩展长度L
#(1)扩展序列需要先进行判断，防止数组越界(2)序列以fasta格式命名，并以名字中是否存在*作为正负样本判据
#按照修饰位点扩展序列获得正样本pos[*]
#按照[Q-pos]扩展序列获得负样本neg

from data_change_1 import Data_change
import random
import numpy as np

class Data_process():

    def __init__(self,file_path,mod_type,L,out_file_name,Dict):
        # 文件路径
        self.file_path = r'fasta_acd/'+file_path
        #修饰类型
        self.mod_type = mod_type
        #窗口大小
        self.L = L
        #蛋白质的修饰类型和修饰位点
        self.Dict = Dict
        #输出fasta文件的路径和文件名
        self.out_file_path = r'data_fasta/' + out_file_name

    def read_fasta(self):
        Name_dic = self.Dict
        file_path = self.file_path
        with open(file_path,'r') as f:
            lines = f.readlines()
        NAME = []
        SEQ = []
        i = 0
        while i < len(lines):
            line = lines[i].strip().split("\t")  # 0:> 1:UniprotID
            if '>' in line:
                ID = line[1]
                seq = lines[i+1].strip()
                for key in Name_dic[ID].keys():
                    for j in Name_dic[ID][key]:
                        name = '>'+'\t'+ID+'\t'+str(j)+'\t'+key# 0:> 1:UniprotID 2:位点 3:修饰类型
                        NAME.append(name)
                        SEQ.append(seq)
            i = i+2

        return NAME,SEQ



    def select_mod(self):
        mod_type = self.mod_type
        NAME, SEQ = self.read_fasta()
        pos_names = []
        pos_seqs = []
        neg_names = []
        neg_seqs = []
        short = []
        for i in range(len(NAME)):
            name = NAME[i].strip().split("\t")# 0:> 1:UniprotID 2:位点 3:修饰类型
            if name[-1] == mod_type:
                NAME_label = NAME[i] + '\t' + 'pos' # # 0:> 1:UniprotID 2:位点 3:修饰类型 4：标签
                pos_names.append(NAME_label)
                pos_seqs.append(SEQ[i])
                short_name_pos = name[0]+'\t'+name[1]+'\t'+name[2]
                short.append(short_name_pos)

        for i in range(len(NAME)):
            name = NAME[i].strip().split("\t")# 0:> 1:UniprotID 2:位点 3:修饰类型
            if name[-1] != mod_type:
                short_name_neg = name[0]+'\t'+name[1]+'\t'+name[2]
                if short_name_neg not in short:
                    NAME_label = NAME[i]+'\t'+'neg' # 0:> 1:UniprotID 2:位点 3:修饰类型 4：标签
                    neg_names.append(NAME_label)
                    neg_seqs.append(SEQ[i])

        return pos_names,pos_seqs,neg_names,neg_seqs


    def make_pos_neg(self):
        L = self.L
        pos_names, pos_seqs, neg_names, neg_seqs = self.select_mod()

        pos = []
        for i in range(len(pos_names)):
            name = pos_names[i].strip().split("\t")# 0:> 1:UniprotID 2:位点 3:修饰类型 4：标签
            seq = 'O'*L+pos_seqs[i]+'O'*L
            site = int(name[2])+L
            pos_short = seq[site-L-1:site+L]
            pos.append(pos_short)
        print('****正样本：',len(pos),'****')

        neg = []
        for i in range(len(neg_names)):
            name = neg_names[i].strip().split("\t")# 0:> 1:UniprotID 2:位点 3:修饰类型 4：标签
            seq = 'O' * L + neg_seqs[i] + 'O' * L
            site = int(name[2]) + L
            neg_short = seq[site - L - 1:site + L]
            neg.append(neg_short)
        print('****负样本：',len(neg),'****')

        return pos,neg
    def balance(self):
        pos_names, _, neg_names, _ = self.select_mod()
        pos, neg = self.make_pos_neg()

        pos_num = len(pos)
        neg_num = len(neg)

        if pos_num < neg_num:
            arr = np.arange(neg_num)
            selected_arr = random.sample(list(arr), pos_num)
            neg_names_new = []
            neg_new = []
            for i in selected_arr:
                neg_names_new.append(neg_names[i])
                neg_new.append(neg[i])
            neg_names = neg_names_new
            neg = neg_new
        else:
            arr = np.arange(pos_num)
            selected_arr = random.sample(list(arr), neg_num)
            pos_names_new = []
            pos_new = []
            for i in selected_arr:
                pos_names_new.append(pos_names[i])
                pos_new.append(pos[i])
            pos_names = pos_names_new
            pos = pos_new

        print('****正样本：',len(pos),'****')
        print('****负样本：',len(neg),'****')

        return pos_names,pos,neg_names,neg

    def out(self):
        # pos_names, _, neg_names, _ = self.select_mod()
        # pos, neg = self.make_pos_neg()
        pos_names, pos, neg_names, neg = self.balance()

        out_file_path = self.out_file_path
        with open(out_file_path,'w') as f:
            for i in range(len(pos_names)):
                print(pos_names[i],file=f)
                print(pos[i], file=f)
            for i in range(len(neg_names)):
                print(neg_names[i],file=f)
                print(neg[i], file=f)

        print("OK!")

if __name__ == '__main__':

    #Camellia sinensis.txt,Oryza sativa.txt
    data1 = Data_change(file_path='Camellia sinensis.txt', data_name = 'CS.txt')
    Name_dic = data1.out_fasta()
    #Acetylation,Crotonylation,Ubiquitination
    data2 = Data_process(file_path='CS.txt',mod_type='Acetylation',L=20,out_file_name='OS_ub.txt',Dict=Name_dic)
    #NAME,SEQ = data2.read_fasta()
    data2.out()
    

