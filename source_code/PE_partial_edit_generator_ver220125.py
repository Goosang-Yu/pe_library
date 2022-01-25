import sys
import pandas as pd
from itertools import product, permutations, combinations
from tqdm import tqdm

## Partially edited sequence generator ##
## Version.1.0 by Goosang Yu (2022.01.11.)
## Contact: gsyu93@gmail.com

df = pd.read_csv(r'E:\Dropbox\Paired lib_PE2\작업 중\PE_off_group5_6\PE_off_group56_list.csv')

list_all = []

for idx in tqdm(df.index, desc='Processing', ncols=100):
    seq_on = df['On-target'].iloc[idx]
    seq_MM = df['WT context'].iloc[idx]
    len_tgt = len(seq_on)

    ## Make mismatch info dictionary ##
    ## Dictionary / Key = mismatched position / value = nt of on-target at that position ##
    dict_edit = {25: 'c'}  # including 25 G to C intended edit as a default
    for pos in range(0, len_tgt):
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
df_with_partial_edit.to_csv(r'E:\Dropbox\Paired lib_PE2\작업 중\PE_off_group5_6\PE_off_group56_list_wPartial.csv', sep=',',
                            index=False)
