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
        print(df_by_mm_num)

        list_brcd = list(df_by_mm_num['Brcd'].drop_duplicates())

        for brcd in list_brcd:
            data = df_by_mm_num.loc[brcd]







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
