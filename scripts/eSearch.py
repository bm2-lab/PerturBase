import subprocess, os, sys
from collections import defaultdict
import numpy as np, pandas as pd
import glob
from multiprocessing import Pool
from sklearn.utils import shuffle
from tqdm import tqdm
import re

#### get bioproject id with study title

def download(x):
    ### cmd = r'esearch -db bioproject  -query "{}[Title]" |  efetch -format   medline'.format(x) 
    ### cmd = r'esearch -db BioProject  -query "{}[Title]" |  efetch -format   medline'.format(x)
    cmd = r'esearch -db BioProject  -query "{}"" |  efetch -format   medline'.format(x)
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