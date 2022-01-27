import sys
import pandas as pd

InputDIR = r'E:\Dropbox\Paired lib_PE2\Fig source data\PE_off-target\현종 정리 중\220118_group56_partial_edit\3_Analysis'

df = pd.read_csv('%s/python_code_for_analysis/MM6_filtered_output.csv' % InputDIR)

## output dictionary - subdevided by RTT lentgh (12, 20, 30)
dict_out = {'All': {}, 12: {}, 20: {}, 30: {}}
for i in range(1,48): dict_out['All'][i] = []
for i in range(1,30): dict_out[12][i] = []
for i in range(1,38): dict_out[20][i] = []
for i in range(1,48): dict_out[30][i] = []

for idx in df.index:

    ## Determine component
    data = df.loc[idx]

    mismatch_info = data['MMpos']
    mismatch_pos = mismatch_info.split('_', )[:-1]
    list_mismatch_pos = list(map(int, mismatch_pos))
    mismatch_num = len(list_mismatch_pos)

    on_target = data['On-trgt']
    WT_seq = data['WTSeq']
    efficiency = data['Relative_adj_PE_efficiency']
    RT = 0

    ## Determine RTT length
    if   len(on_target) == 36: RT = 12
    elif len(on_target) == 44: RT = 20
    elif len(on_target) == 54: RT = 30
    else: print('Target length error!', data['Brcd'])


    for n in range(0, mismatch_num):
        mm_pos = list_mismatch_pos[n] - 3
        if mm_pos not in dict_out[RT]: dict_out[RT][mm_pos] = []
        dict_out[RT][mm_pos].append(efficiency)
        dict_out['All'][mm_pos].append(efficiency)

df_All = pd.DataFrame(dict([(key, pd.Series(val)) for key, val in dict_out['All'].items()]))
df_R12 = pd.DataFrame(dict([(key, pd.Series(val)) for key, val in dict_out[12].items()]))
df_R20 = pd.DataFrame(dict([(key, pd.Series(val)) for key, val in dict_out[20].items()]))
df_R30 = pd.DataFrame(dict([(key, pd.Series(val)) for key, val in dict_out[30].items()]))

df_All.to_csv('%s/python_code_for_analysis/Sort_by_pos_All_MM%s.csv' % (InputDIR, mismatch_num))
df_R12.to_csv('%s/python_code_for_analysis/Sort_by_pos_R12_MM%s.csv' % (InputDIR, mismatch_num))
df_R20.to_csv('%s/python_code_for_analysis/Sort_by_pos_R20_MM%s.csv' % (InputDIR, mismatch_num))
df_R30.to_csv('%s/python_code_for_analysis/Sort_by_pos_R30_MM%s.csv' % (InputDIR, mismatch_num))

## 각 dictionary에 RTT 길이별로 efficiency들 모아서
## output file로 만들기
## 각각의 regione에 상관 없이 position들마다 efficiency를 넣으므로
## 실제 output이 만들어지면, 전체 input의 n배의 efficieny 개수가 나와야 정상이다. n = mismatch 개수

