#!/home/hkim/anaconda3/bin/python

import os, sys, pickle, time, subprocess, json, re, string, regex, random, collections, itertools
import subprocess as sp
import numpy as np
import scipy.stats as stats
import matplotlib as mpl
import pandas as pd
import multiprocessing as mp
from Bio import SeqIO
from tqdm import tqdm

## Original code: B_main_B206.py by JM Park (B206 server)
## Version 22.01.12.


sBASE_DIR = '/extdata1/Jinman/pe_screening'
sSRC_DIR = '/extdata1/GS/FASTQ_read_count'
sOut_BASE = '/extdata1/GS/FASTQ_read_count/output'
sDATA_DIR = '/extdata1/Jinman'

sTIME_STAMP = '%s' % (time.ctime().replace(' ', '-').replace(':', '_'))
sFLASH = '/home/hkim/bin/FLASH-1.2.11-Linux-x86_64/flash'


class cAbsoluteCNV: pass


class cBarcodeData: pass


class cPEData: pass

def main():

    ## For test ##
    # sSRC_DIR = '/home/hkim/src/pe_screening'
    # sBASE_DIR = r'C:\Users\gsyu9\Desktop\BC_sorting_test_220114'
    # sOut_BASE = r'C:\Users\gsyu9\Desktop\BC_sorting_test_220114\output'

    sAnalysis = 'GSY_211008'

    ## Multiplex Options ##
    nCores = 40
    nBins = 40

    ## Main Input Dir ##
    sInputDir = '%s/input' % sBASE_DIR

    ## Rename SamplesIDs -> SampleName ##
    sDataDir = '%s/%s' % (sInputDir, sAnalysis)
    sRenameList = '%s/RenameLists/%s.txt' % (sInputDir, sAnalysis)
    # rename_samples (sDataDir, sRenameList)

    ## Run FLASH ##
    sFileList = '%s/FileLists/%s.txt' % (sInputDir, sAnalysis)
    dict_sFiles = load_NGS_files(sFileList)
    # run_FLASH (sAnalysis, sDataDir, dict_sFiles, bTestRun)

    ## Target Barcode RE ##
    dict_sRE = {'ClinVar': '[ACGT]{6}[T]{6}[ACGT]{24}', 'Profiling': '[ACGT]{6}[T]{6}[ACGT]{24}',
                'Subpool': '[ACGT]{6}[T]{6}[ACGT]{24}', 'D4D21': '[ACGT]{6}[T]{6}[ACGT]{24}',
                'JustNGSReads': '[ACGT]{6}[T]{6}[ACGT]{23}', 'TTTTBC(18+2)': '[T]{6}[ACGT]{20}',
                'TTTTBC(14+2)': '[T]{6}[ACGT]{16}', 'PECV6K_ver220117': '[T]{6}[ACGT]{16}',
                'PE_off_partial_ver220117': '[T]{6}[ACGT]{20}',
                'Offtarget': '[ACGT]{6}[T]{6}[ACGT]{23}', 'Offtarget2': '[T]{6}[ACGT]{20}',
                'Offtarget2-Intended': '[T]{6}[ACGT]{20}', 'Offtarget3': '[T]{6}[ACGT]{20}',
                'Offtarget2-Intended-Test': '[T]{6}[ACGT]{20}',
                }

    ## Error Type: ErrorProne, ErrorFree
    # list_sErrorType = ['ErrorFree', 'ErrorProne']
    list_sErrorType = ['ErrorFree']
    # list_sErrorType = ['ErrorProne']

    ## Guide Run List: Profiling, ClinVar
    list_sGuideRun = ['PE_off_partial_ver220117']

    # list_sGuideRun  = ['D4D21'] # HKK_210405
    # list_sGuideRun  = ['JustNGSReads'] # HKK_201126 and HKK_210405

    for sRun in list_sGuideRun:
        sRE = dict_sRE[sRun]

        for sError in list_sErrorType:
            sRunName = '%s_%s' % (sRun, sError)
            sOutputDir = '%s/%s_%s' % (sOut_BASE, sAnalysis, sRun)
            os.makedirs(sOutputDir, exist_ok=True)

            ## Load Barcode Data ##
            sBarcodeFile = '%s/input/BarcodeTargets/%s_Barcode_Targets_%s.csv' % (sSRC_DIR, sAnalysis, sRun)
            dict_cPE, dict_sPar = load_PE_barcode_list(sBarcodeFile)

            ## Run Analysis ##
            list_sSamples = list(dict_sFiles.keys())

            for sSample in list_sSamples:
                print('Processing %s %s - %s' % (sAnalysis, sRunName, sSample))

                sInDir = '%s/%s/flash/%s' % (sInputDir, sAnalysis, sSample)
                sFastqTag = '%s.extendedFrags' % sSample
                # split_fq_file(sInDir, sFastqTag, bTestRun)
                # gzip_fastq_list (sInDir, sFastqTag, bTestRun)
                # nLineCnt        = get_line_cnt (sInDir, sFastqTag) # set On in real
                list_sSplitFile = get_split_list(sInDir, sFastqTag)
                mp_sort_by_barcode(nCores, sSample, sInDir, sOutputDir, dict_cPE, dict_sPar, list_sSplitFile, sRE, sError)

                nSplitNo = 10
                combine_output_pickle_new(sSample, sOutputDir, sError, dict_cPE, dict_sPar, list_sSplitFile, nSplitNo)


            # loop END: sSample

        # loop END: sError

    # loop END: sRun

# def END: main_test

def combine_output_pickle_new(sSample, sOutputDir, sError, dict_cPE, dict_sPar, list_sSplitFiles, nSplitNo):
    print('combine_output - %s' % sSample)
    sOutDir = '%s/%s' % (sOutputDir, sError)
    os.makedirs(sOutDir, exist_ok=True)

    sCombineOut = '%s/combinedFreq' % sOutDir
    sOthersOut = '%s/others' % sOutDir
    sTempDir = '%s/temp' % sOutDir

    os.makedirs(sCombineOut, exist_ok=True)
    os.makedirs(sOthersOut, exist_ok=True)

    nFileCnt = len(list_sSplitFiles)
    list_nBins = [[int(nFileCnt * (i + 0) / nSplitNo), int(nFileCnt * (i + 1) / nSplitNo)] for i in range(nSplitNo)]

    dict_combined_count = {}

    for sBarcode in dict_sPar.keys():
        cPE = dict_cPE[sBarcode]
        dict_combined_count[sBarcode] = {}

        for ParSeq in dict_sPar[sBarcode]:
            dict_combined_count[sBarcode][ParSeq] = 0


    for nStart, nEnd in list_nBins:
        list_sSubSplit = list_sSplitFiles[nStart:nEnd]

        for i, sSplitFile in enumerate(list_sSubSplit):
            # print('\nOutput Temp Pickle %s/%s -- %s' % ((i + nStart), nEnd, sSplitFile))
            sSplitTag = '_'.join(sSplitFile.split('.')[1].split('_')[-2:])
            sTempFile = '%s/%s/%s' % (sTempDir, sSample, sSplitTag)
            sInFile = '%s/%s.temp_result.csv' % (sTempFile, sSplitTag)

            df_temp_count = pd.read_csv(sInFile, names=['Barcode', 'OnSeq', 'WTSeq', 'SeqProduct', 'read_count'])

            for idx in df_temp_count.index:
                temp_count_bc = df_temp_count['Barcode'].iloc[idx]
                temp_count_seqProd = df_temp_count['SeqProduct'].iloc[idx]
                temp_count_count = df_temp_count['read_count'].iloc[idx]

                dict_combined_count[temp_count_bc][temp_count_seqProd] += temp_count_count

        # loop END: sSplitFile

    df_final = pd.concat({k: pd.Series(v) for k, v in dict_combined_count.items()}).reset_index()
    print('-' * 50, '\nMaking Final output file\n', '-' * 50)
    df_final.to_csv('%s/%s_count_result.csv' % (sCombineOut, sSample))


# def END: combine_output



def reverse_complement(sSeq):
    dict_sBases = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '.': '.', '*': '*',
                   'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    list_sSeq = list(sSeq)  # Turns the sequence in to a gigantic list
    list_sSeq = [dict_sBases[sBase] for sBase in list_sSeq]
    return ''.join(list_sSeq)[::-1]


# def END: reverse_complement

def load_NGS_files(sInFile):
    InFile = open(sInFile, 'r')
    list_sOutput = [sReadLine.strip('\n') for sReadLine in InFile if not sReadLine.startswith('#')]
    InFile.close()

    dict_sOutput = {}
    for sFile in list_sOutput:

        sFileName = sFile.split('.')[0]
        nSampleNo = int(sFileName.split('_')[-1])
        sSample = '_'.join(sFileName.split('_')[:-1])

        if sSample not in dict_sOutput:
            dict_sOutput[sSample] = {}
        if nSampleNo not in dict_sOutput[sSample]:
            dict_sOutput[sSample][nSampleNo] = ''
        dict_sOutput[sSample][nSampleNo] = sFile

    # loop END: sFile

    return dict_sOutput


# def END: load_NGS_files


def get_line_cnt(sWorkDir, sFastqTag):
    sInFile = '%s/%s.fastq' % (sWorkDir, sFastqTag)
    sOutFile = '%s/%s.linecnt.txt' % (sWorkDir, sFastqTag)

    if os.path.isfile(sOutFile):
        InFile = open(sOutFile, 'r')
        nLineCnt = int([sReadLine.strip('\n').split(' ') for sReadLine in InFile][0][0])
        InFile.close()
        print('Fastq Line Cnt', nLineCnt)
        return nLineCnt
    else:
        sCmd = 'wc -l %s > %s' % (sInFile, sOutFile)
        os.system(sCmd)
        InFile = open(sOutFile, 'r')
        nLineCnt = int([sReadLine.strip('\n').split(' ') for sReadLine in InFile][0][0])
        InFile.close()
        return nLineCnt


# def END: get_line_cnt


def get_split_list(sWorkDir, sFastqTag):
    sInFile = '%s/%s.split.txt' % (sWorkDir, sFastqTag)
    InFile = open(sInFile, 'r')
    list_sSplits = [sReadLine.strip('\n') for sReadLine in InFile if not sReadLine.startswith('#')]
    InFile.close()
    return list_sSplits


# def END: get_split_list


def mp_sort_by_barcode(nCores, sSample, sInDir, sOutputDir, dict_cPE, dict_sPar, list_sSplits, sRE, sError):
    sOutDir = '%s/%s' % (sOutputDir, sError)
    os.makedirs(sOutDir, exist_ok=True)

    list_sParameters = []
    for sSplitFile in list_sSplits:
        # PE_off_target_293T_1.extendedFrags_|fastq_0001|.fq
        sSplitTag = '_'.join(sSplitFile.split('.')[1].split('_')[-2:])  # fastq_0001
        sInFile = '%s/split/%s' % (sInDir, sSplitFile)  # PE_off_target_293T_1.extendedFrags_fastq_0001.fq
        sTempDir = '%s/temp/%s/%s' % (sOutDir, sSample, sSplitTag)
        os.makedirs(sTempDir, exist_ok=True)

        list_sParameters.append([sSplitTag, sInFile, sTempDir, dict_cPE, dict_sPar, sRE, sError])
    # loop END: sSplitFile

    p = mp.Pool(nCores)
    p.map_async(sort_by_barcode_all_partial_edit, list_sParameters).get()
    p.map_async(determine_output_for_partial_edit, list_sParameters).get()


# def END: mp_sort_by_barcode


def sort_by_barcode_all_partial_edit(list_sParameters):
    sSplitTag = list_sParameters[0]
    sInFile = list_sParameters[1]
    sTempOut = list_sParameters[2]
    dict_cPE = list_sParameters[3]
    dict_sPar = list_sParameters[4]
    sRE = list_sParameters[5]
    sError = list_sParameters[6]

    nRefBuffer = 24  # Barcode length to subtract from back of RefSeq  ## For offtarget
    nTargetBuffer = 2  # Barcode length to subtract from front of TarSeq ## For offtarget

    list_sBarcodes = list(dict_cPE.keys())
    nTotalCnt = len(list_sBarcodes)

    print('Bio.SeqIO Version - Sort by Barcode Running - %s' % (list_sParameters[0]))

    dict_sOutput = {}
    dict_sOutput2 = {}
    InFile = open(sInFile, 'r')
    for sSeqData in SeqIO.parse(InFile, 'fastq'):
        sReadID = str(sSeqData.id)
        sNGSSeq = str(sSeqData.seq)

        for sReIndex in regex.finditer(sRE, sNGSSeq, overlapped=True):
            nIndexStart = sReIndex.start()
            nIndexEnd = sReIndex.end()
            sBarcodeMatch = sNGSSeq[nIndexStart + nTargetBuffer:nIndexEnd]
            sRefSeqCheck = sNGSSeq[:nIndexStart]
            sTargetSeq = sNGSSeq[nIndexEnd:]

            if nRefBuffer == 0: nRefBuffer = (len(sNGSSeq) - nIndexStart - 1)

            ### Skip Non-barcodes ###
            try:
                cPE = dict_cPE[sBarcodeMatch]
            except KeyError:
                continue
            #########################

            ## Skip error in Refseq ##
            if sError == 'ErrorFree':
                # print('sRefSeqCheck', sRefSeqCheck)
                # print('cPE.sRefSeq[:-nRefBuffer]', cPE.sRefSeq[:-nRefBuffer])
                if not cPE.sRefSeq[:-nRefBuffer] in sRefSeqCheck: continue # 1-2 nt added sequence is correct!
            ##########################

            if sBarcodeMatch not in dict_sOutput:
                dict_sOutput[sBarcodeMatch] = []
            dict_sOutput[sBarcodeMatch].append([sReadID, sNGSSeq, nIndexEnd - nTargetBuffer])

            if sBarcodeMatch not in dict_sOutput2:
                dict_sOutput2[sBarcodeMatch] = []
            dict_sOutput2[sBarcodeMatch].append(sSeqData)
        # loop END: i, sReadLine
    # loop END: sSeqData
    InFile.close()
    print('%s Found= %s' % (sError, len(dict_sOutput)))

    ## Pickle Out #### dict_sOutput <- change not to dictionary ????
    sOutFile = '%s/%s.data' % (sTempOut, sSplitTag)
    OutFile = open(sOutFile, 'wb')
    pickle.dump(dict_sOutput, OutFile)
    OutFile.close()

    sOutFile = '%s/%s.vSeqIO.data' % (sTempOut, sSplitTag)
    OutFile = open(sOutFile, 'wb')
    pickle.dump(dict_sOutput2, OutFile)
    OutFile.close()

    sTempOut_bybar = '%s/bybarcodes' % sTempOut
    os.makedirs(sTempOut_bybar, exist_ok=True)

    '''
    ## Output By Barcode ##
    for sBarcode in dict_sOutput2:
        sOutFile = '%s/%s.fastq' % (sTempOut_bybar, sBarcode)
        OutFile = open(sOutFile, 'w')

        # sOut     = '%s\t%s\n' % (sBarcode, ','.join(dict_sOutput2[sBarcode]))
        # OutFile.write(sOut)

        for sSeqData in dict_sOutput2[sBarcode]:
            try:
                SeqIO.write(sSeqData, OutFile, 'fastq')
            except TypeError:
                continue
        # loop END: sReadID, sNGSRead

        OutFile.close()

    # loop END: sBarcode
    '''
    print('Bio.SeqIO Version - Sort by Barcode DONE - %s' % (list_sParameters[0]))

# def END: sort_by_barcode_v3


def determine_output_for_partial_edit(list_sParameters):

    print('Determine Output - %s' % list_sParameters[2])

    sSplitTag = list_sParameters[0]
    sInFastq = list_sParameters[1]
    sTempOut = list_sParameters[2]
    dict_cPE = list_sParameters[3]
    dict_sPar = list_sParameters[4]

    ## Pickle Load ##
    sInFile = '%s/%s.data' % (sTempOut, sSplitTag)
    InFile = open(sInFile, 'rb')
    dict_sBarcodes = pickle.load(InFile)

    InFile.close()
    print('%s dict_sBarcodes= %s' % (sSplitTag, len(dict_sBarcodes)))

    dict_sOutput = {}

    for sBarcode in dict_sBarcodes: # dict_sBarcodes = fastq file

        cPE = dict_cPE[sBarcode]

        if sBarcode not in dict_sOutput:
            dict_sOutput[sBarcode] = {}
            for partial_seq in dict_sPar[sBarcode]:
                dict_sOutput[sBarcode][partial_seq] = []

        for sReadID, sNGSSeq, sIndexS in dict_sBarcodes[sBarcode]: # python version > 3.6 for dictionary key order

            sTargetSeq = reverse_complement(sNGSSeq[sIndexS:-22])

            for SeqProduct in dict_sOutput[sBarcode].keys(): # dict_sOutput = barcode from fastq file
                if SeqProduct not in sTargetSeq:
                    continue
                elif SeqProduct in sTargetSeq:
                    dict_sOutput[cPE.sBarcode][SeqProduct].append(sNGSSeq)
                    break
            if SeqProduct == 'Other':
                dict_sOutput[sBarcode][SeqProduct].append(sNGSSeq)
        # loop END: sReadID, sNGSSeq
    # loop END: sBarcode

    list_barcode = []
    list_OnSeq = []
    list_WTSeq = []
    list_SeqProduct = []
    list_read_count = []

    for sBarcode in dict_sOutput:
        cPE = dict_cPE[sBarcode]
        for SeqProduct in dict_sOutput[sBarcode]:
            read_count = len(dict_sOutput[sBarcode][SeqProduct])
            list_barcode += [sBarcode]
            list_OnSeq += [cPE.sOn_Seq]
            list_WTSeq += [cPE.sWTSeq]
            list_SeqProduct += [SeqProduct]
            list_read_count += [read_count]

    All_data = zip(list_barcode, list_OnSeq, list_WTSeq, list_SeqProduct, list_read_count)
    sOutFile = pd.DataFrame(All_data)
    sOutFile.to_csv('%s/%s.temp_result.csv' % (sTempOut, sSplitTag), index=False, header= None)


# def END: determine_output_vOfftarget2_Intended

def load_PE_barcode_list(sInFile):  # sInFile = sBarcodeFile
    dict_sOutput = {}
    dict_bc_parSeq = {}
    dict_InFile_columns = {}

    df_barcode_list = pd.read_csv(sInFile, dtype=str)
    count = 0

    for column_name in list(df_barcode_list):
        dict_InFile_columns[count] = column_name
        count += 1

    for idx in tqdm(df_barcode_list.index, desc='Making barcode dictionary',
                    ncols=100): 
        cPE = cPEData()

        cPE.sBarcode = df_barcode_list[dict_InFile_columns[6]].iloc[idx]
        cPE.sRefSeq = df_barcode_list[dict_InFile_columns[7]].iloc[idx].upper().replace('N', '')
        cPE.sOn_Seq = df_barcode_list[dict_InFile_columns[8]].iloc[idx]
        cPE.sWTSeq = df_barcode_list[dict_InFile_columns[9]].iloc[idx]

        sKey = cPE.sBarcode

        list_partial_ED_seq = []
        for i in range(10, len(dict_InFile_columns)):
            sColumn_Seq = df_barcode_list[str(dict_InFile_columns[i])].iloc[idx]
            if isinstance(sColumn_Seq, str):
                list_partial_ED_seq.append(sColumn_Seq[2:-1].upper())
            else:
                break

        dict_bc_parSeq[sKey] = [cPE.sWTSeq[2:-1].upper()] + list_partial_ED_seq + ['Other']

        if sKey not in dict_sOutput:
            dict_sOutput[sKey] = ''
        dict_sOutput[sKey] = cPE
    # loop END:

    return dict_sOutput, dict_bc_parSeq


# def END: load_PE_barcode_list


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
