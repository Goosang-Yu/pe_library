import Bio.SeqIO, os
import pandas as pd
import sys, time, regex
from tqdm import tqdm

start = time.time()

def main():
    sAnalysis_Tag = '63_GS_PE off-target_283T_2_1rxn_220118'
    BaseDIR = r'C:\Users\home\Desktop\220128_miniseq'
    FASTQ_file   = r'%s\%s\%s.fastq' % (BaseDIR, sAnalysis_Tag, sAnalysis_Tag.split('_')[0])
    Barcode_file = r'%s\C_PECV6K_with_refseq_211112.csv' % BaseDIR
    OutDIR = r'%s\%s\output_PECV6K' % (BaseDIR, sAnalysis_Tag)
    sRE = '[T]{4}[ACGT]{16}' ## for PECV6K = 16 / for off-target = 20
    sError = 'ErrorFree'

    dict_brcd = make_bc_list_dictionary(Barcode_file)
    os.makedirs(OutDIR, exist_ok=True)

    dict_brcd_count, dict_read_type_count = find_barcode_in_NGSread(FASTQ_file, sRE, dict_brcd, sError)

    dict_to_csv(dict_brcd_count, OutDIR, 'bc_count_%s' % sAnalysis_Tag, 1)
    dict_to_csv(dict_read_type_count, OutDIR, 'read_count_%s' % sAnalysis_Tag)



def find_barcode_in_NGSread(FASTQ_file, sRE, dict_brcd, sError):
    fastq_info = Bio.SeqIO.parse(FASTQ_file, 'fastq')

    dict_sOutput = {brcd: [] for brcd in dict_brcd.keys()}
    dict_sOutput2 = {brcd: 0 for brcd in dict_brcd.keys()}
    dict_sOutput3 = {'conv': {'WT': 0, 'ED': 0, 'Other': 0},
                     'opti': {'WT': 0, 'ED': 0, 'Other': 0},
                     'Error_prone': {'WT': 0, 'ED': 0, 'Other': 0}}

    for sSeqData in tqdm(fastq_info, desc='Sorting from FASTQ data', ncols=100, total=len(FASTQ_file)/4):

        sReadID = str(sSeqData.id)
        sNGSSeq = str(sSeqData.seq)

        for sReIndex in regex.finditer(sRE, sNGSSeq, overlapped=True):
            nIndexStart = sReIndex.start()
            nIndexEnd = sReIndex.end()
            sBarcodeMatch = sNGSSeq[nIndexStart:nIndexEnd]
            sRefSeqCheck = sNGSSeq[:nIndexStart+24]
            sTargetSeq = sNGSSeq[nIndexEnd-2:-40]

            ### Skip Non-barcodes ###
            try:
                dict_refSeq = dict_brcd[sBarcodeMatch]
            except KeyError:
                continue
            #########################

            ## Skip error in Refseq ##
            if sError == 'ErrorFree':
                if   dict_refSeq['convRef'] in sRefSeqCheck: read_type = 'conv'
                elif dict_refSeq['optiRef'] in sRefSeqCheck: read_type = 'opti'
                else: read_type = 'Error_prone'
            ##########################

            if   dict_brcd[sBarcodeMatch]['WTSeq'] in reverse_complement(sTargetSeq): product_type = 'WT'
            elif dict_brcd[sBarcodeMatch]['EDSeq'] in reverse_complement(sTargetSeq): product_type = 'ED'
            else: product_type = 'Other'

            dict_sOutput2[sBarcodeMatch] += 1
            dict_sOutput3[read_type][product_type] += 1

        # loop END: i, sReadLine
    # loop END: sSeqData

    return dict_sOutput2, dict_sOutput3


def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '.': '.', '*': '*',
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]


def make_bc_list_dictionary(Barcode_file):

    df_bc_list = pd.read_csv(Barcode_file)
    df_bc_list.columns = ['Barcode', 'convRef', 'optiRef', 'WTSeq', 'EDSeq']
    dict_brcd = {}

    for idx in df_bc_list.index:
        data = df_bc_list.loc[idx]
        barcode = data['Barcode'].upper()
        convRef = data['convRef'].upper()
        optiRef = data['optiRef'].upper()
        WT_Seq  = data['WTSeq'].upper()
        ED_Seq  = data['EDSeq'].upper()

        dict_brcd[barcode] = {'convRef': convRef, 'optiRef': optiRef, 'WTSeq': WT_Seq, 'EDSeq': ED_Seq}

    return dict_brcd

def dict_to_csv(dictionary, OutDIR, Output_Tag, T=0):
    df = pd.DataFrame(dict([(key, pd.Series(val)) for key, val in dictionary.items()])).sort_index(axis=1)
    if T == 1: df.T.to_csv('%s/%s.csv' % (OutDIR, Output_Tag))
    else: df.to_csv('%s/%s.csv' % (OutDIR, Output_Tag))



if __name__ == '__main__':
    if len(sys.argv) == 1:
        main()
    else:
        function_name = sys.argv[1]
        function_parameters = sys.argv[2:]
        if function_name in locals().keys():
            locals()[function_name](*function_parameters)
        else:
            sys.exit('ERROR: function_name=%s, parameters=%s' % (function_name, function_parameters))
    # if END: len(sys.argv)
# if END: __name__