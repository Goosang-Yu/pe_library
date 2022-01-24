import sys
import time
import pandas as pd
import numpy as np
from tqdm import tqdm

def main():
    print('''\
    \    /\\
     )  ( ')
    (  /  )
     \(__)|
    '''
          )
    print('\n')
    print('-' * 30, '\n Edit efficiency calculator')
    print('-' * 30)

    ## Directory ##
    BaseDIR   = 'E:/Dropbox/Paired lib_PE2/Fig source data/PE_off-target/현종 정리 중/220118_new_error-free/2_Efficiency_caluration'
    InputDIR  = '%s/input' % BaseDIR
    OutputDIR =  '%s/output' % BaseDIR
    PE_Input_file = 'PE_efficiency_input_293T_PE_off_220118.csv'
    BG_Input_file = 'PE_efficiency_input_BG_PE_off_220118.csv'

    print('PE sample input file:', PE_Input_file)
    print('BG sample input file:', BG_Input_file)

    ## Parameters
    read_threshold = 200
    BG_effi_threshold = 5
    Total_read_type = 'woOther' ## woOther or Total

    ## Original version ##
    print('\nStart processing - On-target X3 version\n')
    df_PE = load_input_and_preprocessing(PE_Input_file, InputDIR)
    df_BG = load_input_and_preprocessing(BG_Input_file, InputDIR)

    df_temp_output = calc_norm_effi(df_PE, df_BG, Total_read_type)
    dict_on_target = make_dict_on_target(df_temp_output, Total_read_type, read_threshold, BG_effi_threshold)

    df_result = calc_relative_effi(df_temp_output, dict_on_target)
    df_result.to_csv('%s/result_output.csv' % OutputDIR)


    ## Combined On target combined read count version ##
    print('\nStart processing - On-target combined read count version\n')
    df_comb_on_PE = combined_on_read_count(df_PE)
    df_comb_on_BG = combined_on_read_count(df_BG)

    df_comb_temp_output = calc_norm_effi(df_comb_on_PE, df_comb_on_BG, Total_read_type)
    dict_comb_on_target = make_dict_on_target(df_comb_temp_output, Total_read_type, read_threshold, BG_effi_threshold)

    df_comb_result = calc_relative_effi(df_comb_temp_output, dict_comb_on_target)
    df_comb_result.to_csv('%s/comb_result_output.csv' % OutputDIR)

# def END: main

def load_input_and_preprocessing(InputFile, InputDIR):

    col_name = ['Barcode', 'WT', 'Alt', 'Intended', 'Mismatch', 'Other', 'Total', 'RT-PBS', 'On_Off']
    df = pd.read_csv('%s/%s' % (InputDIR, InputFile))
    df.columns = col_name
    df['On_Off'] = df['On_Off'].str.upper()

    list_woOther = []
    list_any_ED  = []

    for idx in df.index:
        woOther_cnt = df['Total'].iloc[idx] - df['Other'].iloc[idx]
        any_ED_cnt  = woOther_cnt - df['WT'].iloc[idx]

        list_woOther.append(woOther_cnt)
        list_any_ED.append(any_ED_cnt)

    df.insert(6, 'Any_edit', list_any_ED)
    df.insert(7, 'woOther', list_woOther)

    df.set_index('Barcode', inplace=True)

    return df

# def END: load_input_and_preprocessing

def calc_norm_effi(df_PE, df_BG, Total_read_type):

    list_output = []

    for barcode in tqdm(df_PE.index, desc='Make normalized efficiency table', ncols=100):
        PE = df_PE.loc[barcode]
        BG = df_BG.loc[barcode]

        raw_BG_effi = 0
        raw_PE_effi = 0
        compl_only_effi = 0
        Inten_only_effi = 0
        MisMt_only_effi = 0
        All_edited_effi = 0

        if BG[Total_read_type]: raw_BG_effi = 100 * (BG['Any_edit'] / BG[Total_read_type])
        if PE[Total_read_type]: raw_PE_effi = 100 * (PE['Any_edit'] / PE[Total_read_type])

        if BG[Total_read_type]:
            norm_ALt = (PE[Total_read_type] * (BG['Alt'] / BG[Total_read_type]))
            norm_Int = (PE[Total_read_type] * (BG['Intended'] / BG[Total_read_type]))
            norm_MM = (PE[Total_read_type] * (BG['Mismatch'] / BG[Total_read_type]))
            norm_Any = (PE[Total_read_type] * (BG['Any_edit'] / BG[Total_read_type]))
        else:
            norm_ALt = 0
            norm_Int = 0
            norm_MM = 0
            norm_Any = 0

        if PE[Total_read_type] - norm_ALt:
            compl_only_effi = ((PE['Alt'] - norm_ALt) / (PE[Total_read_type] - norm_ALt)) * 100
        if PE[Total_read_type] - norm_Int:
            Inten_only_effi = ((PE['Intended'] - norm_Int) / (PE[Total_read_type] - norm_Int)) * 100
        if PE[Total_read_type] - norm_MM:
            MisMt_only_effi = ((PE['Mismatch'] - norm_MM) / (PE[Total_read_type] - norm_MM)) * 100
        if PE[Total_read_type] - norm_Any:
            All_edited_effi = ((PE['Any_edit'] - norm_Any) / (PE[Total_read_type] - norm_Any)) * 100

        list_output += [[barcode, PE['RT-PBS'].upper(), PE['On_Off'].upper(),
                         BG['Total'], BG['woOther'], PE['Total'], PE['woOther'],
                         raw_BG_effi, raw_PE_effi, compl_only_effi, Inten_only_effi, MisMt_only_effi, All_edited_effi]]

    df_output = pd.DataFrame(np.array(list_output), columns=['Barcode', 'RT-PBS', 'On_Off', 'BG_Total', 'BG_woOther',
                                                             'PE_Total', 'PE_woOther', 'raw_BG_effi', 'raw_PE_effi',
                                                             'complete_only_effi', 'Int_only_effi',
                                                             'MM_only_effi', 'All_edit_effi'])
    df_output.set_index('Barcode', inplace=True)

    return df_output

# def END: calculate_normalized_efficiency

def make_dict_on_target(df_temp_output, Total_read_type, read_threshold, BG_effi_threshold):

    ## df_on-target filter
    df_on_target = df_temp_output.loc[df_temp_output['On_Off'] == 'ON']

    dict_on_target = {}

    BG_total, PE_total = 'BG_woOther', 'PE_woOther'
    if Total_read_type == 'Total': BG_total, PE_total = 'BG_Total', 'PE_Total'

    for barcode in tqdm(df_on_target.index, desc='Making On-target dictionary', ncols=100):
        data = df_on_target.loc[barcode]
        BG_read, PE_read = float(data[BG_total]), float(data[PE_total])
        sRTPBS = data['RT-PBS']

        if sRTPBS not in dict_on_target:
            dict_on_target[sRTPBS] = {'raw_BG': [], 'raw_PE': [],'Complt': [],'IntOly': [],'MsMOly': [],'AllEdt': []}

        if BG_read >= read_threshold: dict_on_target[sRTPBS]['raw_BG'] += [data['raw_BG_effi']]
        if PE_read >= read_threshold: dict_on_target[sRTPBS]['raw_PE'] += [data['raw_PE_effi']]

        if BG_read >= read_threshold and PE_read >= read_threshold and float(data['raw_BG_effi']) < BG_effi_threshold:
            dict_on_target[sRTPBS]['Complt'] += [data['complete_only_effi']]
            dict_on_target[sRTPBS]['IntOly'] += [data['Int_only_effi']]
            dict_on_target[sRTPBS]['MsMOly'] += [data['MM_only_effi']]
            dict_on_target[sRTPBS]['AllEdt'] += [data['All_edit_effi']]

    return dict_on_target


# def END: make_on_target_dictionary

def calc_relative_effi(df_temp_output, dict_on_target):

    list_brcd = []
    list_effi = []
    list_numb = []
    list_rltv = []

    edit_type = ['raw_BG', 'raw_PE', 'Complt', 'IntOly', 'MsMOly', 'AllEdt']
    effi_col_name = ['AveOnEffi_%s' % c for c in edit_type]
    numb_col_name = ['NumOn_%s' % c for c in edit_type]
    rltv_col_name = ['relative_%s' % c for c in edit_type[2:]]

    for Barcode, sRTPBS in zip(df_temp_output.index, df_temp_output['RT-PBS']):
        on_effi = dict_on_target[sRTPBS]
        list_effi_temp = []
        list_numb_temp = []

        for on_type in edit_type:
            num = int(len(on_effi[on_type]))
            effi = 0
            if num:
                efficiencies = list(map(float, on_effi[on_type]))
                effi = np.mean(efficiencies)

            list_effi_temp.append(effi)
            list_numb_temp.append(num)

        list_brcd.append(Barcode)
        list_effi.append(list_effi_temp)
        list_numb.append(list_numb_temp)

    on_target_array = np.concatenate((np.array(list_effi), np.array(list_numb)), axis=1)
    df_on_target_array = pd.DataFrame(on_target_array, columns=effi_col_name+numb_col_name)
    df_on_target_array['Barcode'] = list_brcd
    df_on_target_array.set_index('Barcode', inplace=True)

    df_result = pd.concat([df_temp_output, df_on_target_array], axis=1).reindex(df_temp_output.index)

    for Barcode in tqdm(df_result.index, desc='Calculate relative efficiencies', ncols=100):
        data = df_result.loc[Barcode]
        Ave_ontgt_effi = float(data['AveOnEffi_AllEdt'])
        rltv_Complt_effi = rltv_Intend_effi = rltv_Mismch_effi = rltv_Any_ed_effi = 0

        if Ave_ontgt_effi:
            rltv_Complt_effi =  (float(data['complete_only_effi']) / Ave_ontgt_effi)
            rltv_Intend_effi =  (float(data['Int_only_effi']) / Ave_ontgt_effi)
            rltv_Mismch_effi =  (float(data['MM_only_effi']) / Ave_ontgt_effi)
            rltv_Any_ed_effi =  (float(data['All_edit_effi']) / Ave_ontgt_effi)

        list_rltv.append([rltv_Complt_effi, rltv_Intend_effi, rltv_Mismch_effi, rltv_Any_ed_effi])

    df_rltv_effi = pd.DataFrame(np.array(list_rltv), columns=rltv_col_name)
    df_rltv_effi['Barcode'] = list_brcd
    df_rltv_effi.set_index('Barcode', inplace=True)
    df_final = pd.concat([df_result, df_rltv_effi], axis=1).reindex(df_result.index)

    return df_final

# def END: calculate_relative_efficiency


def combined_on_read_count(df):

    df_on_target_temp = df.loc[df['On_Off'] == 'ON']
    df_off_target_temp = df.loc[df['On_Off'] == 'OFF']

    dict_on_combined_read = {}

    for barcode in df_on_target_temp.index:
        data = df_on_target_temp.loc[barcode]
        sRTPBS = data['RT-PBS'].upper()
        if sRTPBS not in dict_on_combined_read:
            dict_on_combined_read[sRTPBS] = {'WT': 0, 'Alt': 0, 'Intended': 0, 'Mismatch': 0,'Other': 0,
                                             'Any_edit': 0, 'woOther': 0, 'Total': 0, 'RT-PBS': '', 'On_Off': ''}

        dict_on_combined_read[sRTPBS]['WT'] += int(data['WT'])
        dict_on_combined_read[sRTPBS]['Alt'] += int(data['Alt'])
        dict_on_combined_read[sRTPBS]['Intended'] += int(data['Intended'])
        dict_on_combined_read[sRTPBS]['Mismatch'] += int(data['Mismatch'])
        dict_on_combined_read[sRTPBS]['Other'] += int(data['Other'])
        dict_on_combined_read[sRTPBS]['Total'] += int(data['Total'])

    for sRTPBS in dict_on_combined_read:
        data_on = dict_on_combined_read[sRTPBS]
        data_on['Any_edit'] = int(data_on['Alt']) + int(data_on['Intended']) + int(data_on['Mismatch'])
        data_on['woOther'] = int(data_on['Total']) - int(data_on['Other'])
        data_on['RT-PBS'] = sRTPBS
        data_on['On_Off'] = 'ON'

    df_combined_on_read = pd.DataFrame(dict_on_combined_read).T
    df_result = pd.concat([df_combined_on_read, df_off_target_temp])

    return df_result

# def END: combined_on_read_count


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
