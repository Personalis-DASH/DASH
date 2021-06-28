import pandas as pd
import numpy as np
from functools import partial
import argparse
import os
import sys

AAs_and_mods = set(['C', 'G', 'A', 'T', 'S', 'N', 'D', 'E', 'Q', 'K', 'R', 'V', 'I', 'L', 'M', 'W', 'F', 'Y', 'H', 'P', '*'])

def check_tmt_pep(pep, AAs_and_mods):
    ''' '''
    is_tmt_pep = 'Y'
    pep_new = pep.replace('(+229.16)(+229.16)', '+')
    # print(pep_new)
    pep_new = pep_new.replace('(+229.16)', '*')
    # print(pep_new)
    status = 'ok'

    ## Check permitted AAs
    if pep_new.startswith('K+'):
        # print("replacing first K+ by K*")
        pep_new = pep_new[0] + '*' + pep_new[2:]
        # print(pep_new)

    if not set(pep_new).issubset(AAs_and_mods):
        is_tmt_pep = 'unknown_aas_or_mods'

    ## Check that n-term aa is modified
    if pep_new[1] != '*':
        is_tmt_pep = 'ntermaa_notmodified'

    ## check len
    # print(pep_new)
    pep_new = pep_new.replace('*', '')
    # print(pep_new)
    pep_len = len(pep_new)
    if (pep_len < 8) and (pep_len > 11):
        is_tmt_pep = 'too_long'

    ## check diversity
    breaks = [0]
    for i in np.arange(1, len(pep_new)):
        if pep_new[i - 1] != pep_new[i]:
            breaks.append(i)
    breaks.append(len(pep_new))

    seg_lens = []
    for i in np.arange(1, len(breaks)):
        seg_len = breaks[i] - breaks[i - 1]
        seg_lens.append(seg_len)

    max_seg_len = np.max(seg_lens)
    if max_seg_len >= 4:
        is_tmt_pep = 'long_repeats_gt4'

    ## print result
    # print("is_tmt_pep: %s" % is_tmt_pep)

    return is_tmt_pep


def check_signal_intensity(max_intensity, min_max_intensity):
    ''' '''
    has_min_signal_intensity = 'Y'
    if max_intensity < min_max_intensity:
        has_min_signal_intensity = 'N'

    return has_min_signal_intensity

def log_tx(value):
    ''' '''
    return np.log10(value+10)

def calc_fc(row):
    ''' '''
    #print(row)
    nr, dr = row[0], row[1]
    fc = (nr+0.01)/(dr+0.01)
    return fc

def get_args():
    ''' '''

    paramDict = {}

    parser = argparse.ArgumentParser()
    parser.add_argument("--in_file", help="infile", action='store', required=True)
    parser.add_argument("--out_file", help="outfile", action='store', required=True)
    parser.add_argument("--sep", help="separator in input file", action='store', default=',')
    parser.add_argument("--pep_colname", help="peptide colname", action='store', default='Peptide')
    parser.add_argument("--intensity_colnums", help="intensity columns (zero indexed)", action='store', default='')
    parser.add_argument("--min_max_intensity", help="min_max_intensity", action='store', default=9999)
    parser.add_argument("--adjusted_median_intensity", help="min_max_intensity", action='store', default=10)

    ###### Parse parameters
    args = parser.parse_args()
    for arg in vars(args):
        paramName, paramValue = arg, getattr(args, arg)
        paramDict[paramName] = paramValue

    ###### checks
    if not os.path.isfile(paramDict['in_file']):
        print("[ERROR] Can not find in_file: %s" % paramDict['in_file'])
        sys.exit()
    if os.path.isfile(paramDict['out_file']):
        print("[ERROR] out_file already exists. program does not overwrite. out_file: %s" % paramDict['out_file'])
        sys.exit()

    df_tmp = pd.read_csv(paramDict['in_file'], sep=paramDict['sep'], nrows=1)
    if paramDict['pep_colname'] not in df_tmp.columns:
        print("[ERROR] Can not find pep_colname: %s in in_file. sep= %s" % (paramDict['pep_colname'], paramDict['sep']))
        sys.exit()

    return paramDict



def main():
    ''' '''
    paramDict = get_args()
    print(paramDict)
    check_tmt_pep_partial = partial(check_tmt_pep, AAs_and_mods=AAs_and_mods)

    in_file = paramDict['in_file']
    df = pd.read_csv(in_file, sep=",")
    df['is_tmt_pep'] = df['Peptide'].apply(check_tmt_pep_partial)


    ## filter intensities
    if paramDict['intensity_colnums'] != '':
        intensity_colnums = paramDict['intensity_colnums'].strip().split(',')
        intensity_colnums = [int(item) for item in intensity_colnums]
    else:
        intensity_colnums = []

    ## filter/transform columns
    if len(intensity_colnums) >= 1:

        ## get intensity colnames
        intensity_colnames = []
        for intensity_colnum in intensity_colnums:
            intensity_colname = df.columns[intensity_colnum]
            intensity_colnames.append(intensity_colname)

        ## filter by sum of intensities
        df['max_intensity'] = df[intensity_colnames].max(axis=1)
        check_signal_intensity_partial = partial(check_signal_intensity, min_max_intensity=paramDict['min_max_intensity'])
        df['has_good_signal'] = df['max_intensity'].apply(check_signal_intensity_partial)

        ## transform intensities
        intensity_colnames_tx = []
        for intensity_colname in intensity_colnames:

            intensity_colname_tx = intensity_colname + ' normalized'
            intensity_colnames_tx.append(intensity_colname_tx)
            df[intensity_colname_tx] = df[intensity_colname].apply(log_tx)
            median_intensity = df[ (df['is_tmt_pep'] == 'Y')  & (df['has_good_signal'] == 'Y' )][intensity_colname_tx].median()
            df[intensity_colname_tx] = df[intensity_colname_tx] - median_intensity + paramDict['adjusted_median_intensity']

        ## calculate ratios
        df['fc'] = df[intensity_colnames_tx].apply(calc_fc, axis=1)

    ## Write output
    df.to_csv(paramDict['out_file'], sep=paramDict['sep'], index=False)

    return




if __name__ == '__main__':
    main()

