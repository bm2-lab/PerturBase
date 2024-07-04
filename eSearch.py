import subprocess, os, sys
from collections import defaultdict
import numpy as np, pandas as pd
import glob
from multiprocessing import Pool
from sklearn.utils import shuffle
from tqdm import tqdm
import re

### 存在一个问题就是 有些文章有数据，但是没有projectID，原因是projectID的标题和原文标题没有进行关联，比如
### Pooled Knockin Targeting for Genome Engineering of Cellular Immunotherapies 有数据，但是并没有projectID，其数据存放在 POKI-Seq in Primary Human T Cells文章中



def download(x):
    ### cmd = r'esearch -db bioproject  -query "{}[Title]" |  efetch -format   medline'.format(x) ## 会漏掉一些信息
    ### cmd = r'esearch -db BioProject  -query "{}[Title]" |  efetch -format   medline'.format(x) ## 后续增加该命令判断有没有 geo data
    cmd = r'esearch -db BioProject  -query "{}"" |  efetch -format   medline'.format(x)   ### 试试这个命令
    result = subprocess.getoutput(cmd)
    if result:
        try:
            results = result.strip().split('\n\n\n')
            ID = [re.findall('PRJ\S+\d+', i)[0] for i in results]
            return (x, ID)
        except Exception as e:
            return (x, '')
            print (e, result, cmd)
    else:
        return (x, '')
            
    
        



def f_download(filein1, filein2):
    dat1 = pd.read_csv(filein1, sep=',')
    dat2 = pd.read_excel(filein2)
    mylist = set(list(dat1['Title']) + list(dat2['Article Title']))
    with Pool(3) as pool:  ### 设置太大会报错
        results = list(tqdm(pool.imap(download, mylist), total=len(mylist)))
    with open('result.tsv', 'w') as fout:
        fout.write('Title\tProjectID\n')
        for title, ID in results:
            fout.write('{}\t{}\n'.format(title, ','.join(ID)))

###第一个文件是PubMed, 第二个数据是WoS
if __name__ == '__main__':
    if len(sys.argv[:]) != 3:
        print ('input PubMed file and WoS file\n')
    else:
        f_download(sys.argv[1], sys.argv[2])
    