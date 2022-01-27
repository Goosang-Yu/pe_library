import sys, time
import pandas as pd

def main():
    ## Load input data
    InputDIR = r'E:\Dropbox\Paired lib_PE2\Fig source data\PE_off-target\현종 정리 중\220118_group56_partial_edit\3_Analysis'
    OutputDIR = r'%s/python_code_for_analysis/output' % InputDIR
    Mismatch_Tag = 2 ## set the number of mismatch which needs the output data

    df = pd.read_csv('%s/python_code_for_analysis/MM%s_filtered_output.csv' % (InputDIR, Mismatch_Tag))
    dict_out, dict_out_region = make_dict_sorted_by_mismatched_position_and_region(df, Mismatch_Tag)

    for key in dict_out.keys():
        dict_to_csv(dict_out[key], OutputDIR, 'Sort_by_pos_%s_MM%s' % (key, Mismatch_Tag))
    for key in dict_out_region.keys():
        dict_to_csv(dict_out_region[key], OutputDIR, 'Sort_by_region_%s_MM%s' % (key, Mismatch_Tag))

# def main(): End

def make_dict_sorted_by_mismatched_position_and_region(dataframe, Mismatch_Tag):

    ## Make initial dictionary sorted by position - subdevided by RTT lentgh (12, 20, 30)
    dict_out = {'All': {}, 12: {}, 20: {}, 30: {}}
    for i in range(1,48): dict_out['All'][i] = []
    for i in range(1,30): dict_out[12][i] = []
    for i in range(1,38): dict_out[20][i] = []
    for i in range(1,48): dict_out[30][i] = []

    ## Make initial dictionary sorted by region - subdevided by RTT lentgh (12, 20, 30)
    region1_range = list(range(1, 7))
    region2_range = list(range(7, 18))
    region3_range = list(range(18, 21))
    region4_range = list(range(21, 30))

    dict_out_region = {'All': {}, 12: {}, 20: {}, 30: {}}

    for idx in dataframe.index: ## Data sorting

        ## Determine component
        data = dataframe.loc[idx]
        mismatch_info = data['MMpos']
        list_mismatch_pos = list(map(int, mismatch_info.split('_', )[:-1]))
        mismatch_num = len(list_mismatch_pos)
        if mismatch_num != Mismatch_Tag:
            print('Mismatch number error!')
            sys.exit()

        efficiency = data['Relative_adj_PE_efficiency']
        RT = data['RTT']
        list_region = []

        for n in range(0, mismatch_num):
            mm_pos = list_mismatch_pos[n] - 3
            if mm_pos not in dict_out[RT]: dict_out[RT][mm_pos] = []
            dict_out[RT][mm_pos].append(efficiency)
            dict_out['All'][mm_pos].append(efficiency)

            if   mm_pos in region1_range: list_region.append('R1')
            elif mm_pos in region2_range: list_region.append('R2')
            elif mm_pos in region3_range: list_region.append('R3')
            elif mm_pos in region4_range: list_region.append('R4')
            else: print('Can not find exact region?!')

        ## Remove duplicated region
        list_region = remove_duplicated_component_in_list(list_region)

        label_region = ''
        for n in range(0, len(list_region)):
            label_region = '%s%s' % (label_region, list_region[n])

        if label_region not in dict_out_region['All']: dict_out_region['All'][label_region] = []
        if label_region not in dict_out_region[RT]: dict_out_region[RT][label_region] = []
        dict_out_region['All'][label_region].append(efficiency)
        dict_out_region[RT][label_region].append(efficiency)

    return dict_out, dict_out_region



def remove_duplicated_component_in_list(list_for_processing):
    list_new = []
    for component in list_for_processing:
        if component not in list_new: list_new.append(component)

    return list_new

def dict_to_csv(dictionary, OutDIR, Output_Tag):
    df = pd.DataFrame(dict([(key, pd.Series(val)) for key, val in dictionary.items()])).sort_index(axis=1)
    df.to_csv('%s/%s.csv' % (OutDIR, Output_Tag))

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
