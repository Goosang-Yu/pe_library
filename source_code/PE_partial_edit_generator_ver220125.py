import sys
import pandas as pd
from itertools import product, permutations, combinations
from tqdm import tqdm

## Partially edited sequence generator ##
## Version.1.0 by Goosang Yu (2022.01.11.)
## Contact: gsyu93@gmail.com

base_DIR = r'C:\Users\home\Dropbox\Paired lib_PE2\작업 중\PE_off_group7'
sAnalysis_Tag = 'group7_DNMT1'

df = pd.read_csv(r'%s\input\%s.csv' % (base_DIR, sAnalysis_Tag))

list_all = []

for idx in tqdm(df.index, desc='Processing', ncols=100):
    data = df.loc[idx]
    seq_on = data['On-target']
    seq_MM = data['WT context']
    target = data['Target #']
    len_tgt = len(seq_on)

    ## Make mismatch info dictionary ##
    ## Dictionary / Key = mismatched position / value = nt of on-target at that position ##
    if   target == 'DNMT1': dict_edit = {26: 'c'}  # including 25 G to C intended edit as a default
    elif target == 'EMX1':  dict_edit = {25: 't'}  # including 25 G to C intended edit as a default
    elif target == 'FANCF': dict_edit = {26: 'c'}  # including 25 G to C intended edit as a default
    elif target == 'HEK3':  dict_edit = {21: 'cttT'}  # including 25 G to C intended edit as a default
    elif target == 'HEK4':  dict_edit = {22: 't'}  # including 25 G to C intended edit as a default
    elif target == 'RNF2':  dict_edit = {26: 'a'}  # including 25 G to C intended edit as a default
    elif target == 'RUNX1': dict_edit = {26: 'c'}  # including 25 G to C intended edit as a default
    elif target == 'VEGFA': dict_edit = {25: 't'}  # including 25 G to C intended edit as a default
    else: print('No tharget matched!')

    for pos in range(4, len_tgt - 3):
        if seq_on[pos] != seq_MM[pos]:
            dict_edit[pos] = seq_on[pos]

    ## Make all combinations for partial edit ##

    list_comb = []

    for var in range(0, len(dict_edit.keys())):  # from nC1 to nCn
        comb = list(combinations(dict_edit.keys(), var + 1))
        list_comb += comb

    ## Make all possible sequences ##
    list_partial_edited = []
    for sComKey in list_comb:
        sNewSeq = list(seq_MM)
        for editpos in sComKey:
            sNewSeq[editpos] = dict_edit[editpos]
        sNewSeq = ''.join(sNewSeq)
        list_partial_edited.append(sNewSeq)

    list_all.append(list_partial_edited)
# idx loop end

df_partial = pd.DataFrame(list_all)

df_with_partial_edit = pd.concat([df, df_partial], axis=1)
df_with_partial_edit.to_csv(r'%s/output/%s_wPartial.csv' % (base_DIR, sAnalysis_Tag), sep=',', index=False)
