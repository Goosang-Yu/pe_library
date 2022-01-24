import sys
import time
import pandas as pd
import numpy as np
from tqdm import tqdm

def main():
    print('\nPE off-target partial edit classification\n')

    ## Directory ##
    BaseDIR   = r'E:\Dropbox\Paired lib_PE2\Fig source data\PE_off-target\현종 정리 중\220118_group56_partial_edit\2_group56_partial_edit_sorting'
    InputDIR  = '%s/input' % BaseDIR
    OutputDIR =  '%s/output' % BaseDIR
    Input_file = '%s/group56_partial_count.csv' % InputDIR

    df_raw = pd.read_csv(Input_file)
    df_raw.columns = ['Brcd', 'Trgt', 'PBS', 'RTT', 'MM', 'posEdit', 'RefSeq', 'On-trgt', 'WTSeq', 'cPE', 'cBG']

    list_mismatch_num = list(df_raw['MM'].drop_duplicates())

    print(df_raw)
    print(list_mismatch_num)

    for MM_num in list_mismatch_num:
        print('Filter: Mismatch number', MM_num)
        df_by_mm_num = df_raw.loc[df_raw['MM'] == MM_num]
        df_by_mm_num.reset_index(inplace=True, drop=True)
        # df_by_mm_num.to_csv('%s/MM%s_df_by_MM.csv' % (InputDIR, MM_num)) ## for test / delete later
        print(df_by_mm_num)

        dict_cPE = {}
        dict_cBG = {}
        dict_out = {}

        for idx, brcd in tqdm(enumerate(df_by_mm_num['Brcd']),
                              desc='Make dictionaries', ncols=100, total=len(df_by_mm_num)):

            data = df_by_mm_num.loc[idx]

            if brcd not in dict_out:
                dict_cPE[brcd] = {}
                dict_cBG[brcd] = {}
                dict_out[brcd] = {'Trgt': brcd, 'PBS': data['PBS'], 'RTT': data['RTT'], 'MM': data['MM'],
                                  'posEdit': data['posEdit'], 'On-trgt': data['On-trgt'], 'WTSeq': data['WTSeq']}
                d = dict_out[brcd]
                temp_CltSeq = list(d['WTSeq'][:21] + d['On-trgt'][21:21 + d['RTT']] + d['WTSeq'][21 + d['RTT']:])
                temp_CltSeq[25] = 'c' ## pos +5 G to C transversion
                d['CompltEditSeq'] = ''.join(temp_CltSeq)

            dict_cPE[brcd][data['RefSeq']] = data['cPE']
            dict_cBG[brcd][data['RefSeq']] = data['cBG']

        for brcd in dict_out:
            d = dict_out[brcd]
            ## read count for PE-treated
            tot_cPE = sum(list(dict_cPE[brcd].values()))
            otr_cPE = dict_cPE[brcd]['Other']
            WT_cPE  = dict_cPE[brcd][d['WTSeq']]
            clt_cPE = dict_cPE[brcd][d['CompltEditSeq']]
            par_cPE = tot_cPE - otr_cPE - WT_cPE - clt_cPE

            ## read count for BG
            tot_cBG = sum(list(dict_cBG[brcd].values()))
            otr_cBG = dict_cBG[brcd]['Other']
            WT_cBG  = dict_cBG[brcd][d['WTSeq']]
            clt_cBG = dict_cBG[brcd][d['CompltEditSeq']]
            par_cBG = tot_cBG - otr_cBG - WT_cBG - clt_cBG

            ## Calculate editing efficiency














if __name__ == '__main__':
    if len(sys.argv) == 1:
        start_time = time.perf_counter()
        main()
        print("\n::::::::::: %.2f seconds ::::::::::::::" % (time.perf_counter() - start_time))
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    # if END: len(sys.argv)
# if END: __name__
