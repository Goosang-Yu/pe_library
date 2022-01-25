import sys
import time
import pandas as pd
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

    for MM_num in list_mismatch_num:
        print('-' * 50, '\nFilter: Mismatch number', MM_num, '\n', '-' * 50)
        df_by_mm_num = df_raw.loc[df_raw['MM'] == MM_num]
        df_by_mm_num.reset_index(inplace=True, drop=True)
        # df_by_mm_num.to_csv('%s/MM%s_df_by_MM.csv' % (InputDIR, MM_num)) ## for test / delete later

        dict_cPE = {}
        dict_cBG = {}
        dict_out = {}
        list_err = []

        for idx, brcd in tqdm(enumerate(df_by_mm_num['Brcd']),
                              desc='Make dictionaries', ncols=100, total=len(df_by_mm_num)):

            df = df_by_mm_num.loc[idx]

            if brcd not in dict_out:
                dict_cPE[brcd] = {}
                dict_cBG[brcd] = {}
                dict_out[brcd] = {'Brcd': brcd, 'PBS': df['PBS'], 'RTT': df['RTT'], 'MM': df['MM'],
                                  'posEdit': df['posEdit'], 'On-trgt': df['On-trgt'], 'WTSeq': df['WTSeq'].upper()}
                d = dict_out[brcd]
                temp_CltSeq = list(d['WTSeq'][:21] + d['On-trgt'][21:21 + d['RTT']] + d['WTSeq'][21 + d['RTT']:])
                temp_CltSeq[25] = 'c' ## pos +5 G to C transversion
                d['CompltEditSeq'] = ''.join(temp_CltSeq).upper()

            dict_cPE[brcd][df['RefSeq']] = df['cPE']
            dict_cBG[brcd][df['RefSeq']] = df['cBG']

        # for loop end: df_by_mm_num

        for brcd in dict_out:
            d = dict_out[brcd]
            try:
                ## read count for PE-treated
                tot_cPE = sum(list(dict_cPE[brcd].values()))
                otr_cPE = dict_cPE[brcd]['Other']
                WT_cPE  = dict_cPE[brcd][d['WTSeq']]
                clt_cPE = dict_cPE[brcd][d['CompltEditSeq']] ### key error for complt edit seq -> 아마도 cltSeq 만드는 code가 잘못된듯
                par_cPE = tot_cPE - otr_cPE - WT_cPE - clt_cPE
            except:
                error = ['PE', brcd, d['WTSeq'], d['CompltEditSeq'], d['On-trgt']]
                list_err.append(error)

            ## read count for BG
            try:
                tot_cBG = sum(list(dict_cBG[brcd].values()))
                otr_cBG = dict_cBG[brcd]['Other']
                WT_cBG  = dict_cBG[brcd][d['WTSeq']]
                clt_cBG = dict_cBG[brcd][d['CompltEditSeq']]
                par_cBG = tot_cBG - otr_cBG - WT_cBG - clt_cBG
            except:
                error = ['BG', brcd, d['WTSeq'], d['CompltEditSeq'], d['On-trgt']]
                list_err.append(error)

            ## Calculate editing efficiency

            # init value
            woOther_BG = tot_cBG - otr_cBG
            woOther_PE = tot_cPE - otr_cPE

            if woOther_BG:
                BG_All_Edit_Effi = 100 * ((clt_cBG + par_cBG) / woOther_BG)
            else:
                BG_All_Edit_Effi = 0
            PE_clt_only_effi = 0
            PE_Par_only_effi = 0
            PE_All_edit_effi = 0

            # read count BG normalization

            if woOther_BG:
                norm_Clt = (woOther_PE * (clt_cBG / woOther_BG))
                norm_Par = (woOther_PE * (par_cBG / woOther_BG))
                norm_All = (woOther_PE * ((clt_cBG + par_cBG) / woOther_BG))
            else:
                norm_Clt = 0
                norm_Par = 0
                norm_All = 0

            if woOther_PE - norm_Clt:
                PE_clt_only_effi = ((clt_cPE - norm_Clt) / (woOther_PE - norm_Clt)) * 100
            if woOther_PE - norm_Par:
                PE_Par_only_effi = ((par_cPE - norm_Par) / (woOther_PE - norm_Par)) * 100
            if woOther_PE - norm_All:
                PE_All_edit_effi = (((clt_cPE + par_cPE) - norm_All) / (woOther_PE - norm_All)) * 100

            # make output

            d['Total_cPE'] = tot_cPE
            d['Other_cPE'] = otr_cPE
            d['woOther_cPE'] = woOther_PE
            d['WTSeq_cPE'] = WT_cPE
            d['Complete_cPE'] = clt_cPE
            d['Partial_cPE'] = par_cPE

            d['Total_cBG'] = tot_cBG
            d['Other_cBG'] = otr_cBG
            d['woOther_cBG'] = woOther_BG
            d['WTSeq_cBG'] = WT_cBG
            d['Complete_cBG'] = clt_cBG
            d['Partial_cBG'] = par_cBG

            d['BG_AllEdit']   = BG_All_Edit_Effi # complete + partial edit in BG sample
            d['PE_CmpltEdit'] = PE_clt_only_effi # complete edit (expected product) fraction
            d['PE_PartlEdit'] = PE_Par_only_effi # partial edit (unexpected product; Not WT or complete edit)
            d['PE_AllEdit']   = PE_All_edit_effi # complete + partial edit

            list_mismatch_pos = []

            for pos in range(0, len(d['On-trgt'])):
                on_base = d['On-trgt'][pos]
                mm_base = d['WTSeq'][pos]
                if on_base != mm_base: list_mismatch_pos.append(pos)

            counted_mm = len(list_mismatch_pos)
            mm_pos_info = ''
            for p in list_mismatch_pos:
                mm_pos_info = '%s%s_' % (mm_pos_info, p)

            d['Counted MM number'] = counted_mm
            d['MMpos'] = mm_pos_info

        df_error = pd.DataFrame(list_err, columns=['Sample Type', 'Barcode', 'WTSeq', 'CompleteEditSeq', 'On-target'])
        df_error.to_csv('%s/error_list_MM%s.csv' % (OutputDIR, MM_num), index=False)
        df_result = pd.DataFrame(dict_out).T
        df_result.to_csv('%s/output_MM%s.csv' % (OutputDIR, MM_num), index=False)

        # for loop end: dict_out
    # for loop end: list_mismatch_num

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
