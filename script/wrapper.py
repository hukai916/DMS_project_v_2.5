"""
A Python wrapper of DMS pipeline.
"""
from pathlib import Path
from toolbox_kai import fastq_splitter_wrapper, quality_control_library, parse_bamfile, wt_codon, quality_control_library_double
from toolbox_enrich2 import enrich2_json_encoder, enrich2_hdf5_extractor, enrich2_hdf5_extractor2, tsv_plot_output, tsv_plot_output_aa, enrich2_count_wrapper, tsv_plot_output_double, tsv_plot_output_aa_double, enrich2_json_encoder_count
from ast import literal_eval as make_tuple
import shutil # complement pathlib with copy file utilities
import subprocess
import pandas as pd
import numpy as np
import pysam
import sys
import os
import glob # Unix style pathname pattern expansion, this seems to work better than Path().glob()
from Bio import SeqIO
from constants import CodonList, AA2CODON, AAList, CODON_GROUPS, CODON2AA
from Bio import SeqIO
import argparse

class ConfigParam(): # store the configuration paramters
    def __init__(self, configFile):
        self.ngs_data_local = []
        self.parseConfig(configFile)
    def parseConfig(self, filename):
        for line in open(filename):
            if line.startswith("NGS raw data(if NCBI)"):
                self.ngs_data_ncbi  = line.split(":")[1].split()
            elif line.startswith("Seq_file"):
                self.ngs_data_local.append(line.split(":")[1].strip())
            elif line.startswith("WT template sequence:"):
                self.wtfile = workdir.joinpath(line.split("WT template sequence:")[1].strip())
            elif line.startswith("Experimental conditions:"):
                self.condition = line.split("Experimental conditions:")[1].strip().split()
            elif line.startswith("Amplicon locations:"):
                self.amplicon = list(map(make_tuple, line.split("Amplicon locations:")[1].strip().split()))
            elif line.startswith("Mutation coordinates in WT template:"):
                self.mut_pos  = list(map(make_tuple, line.split("Mutation coordinates in WT template:")[1].strip().split()))
            elif line.startswith("WT masks:"):
                self.wt_mask  = list(line.split("WT masks:")[1].strip().split())
                self.wt_mask  = [[int(item.split('-')[0]), item.split('-')[1]] for item in self.wt_mask]
        self.mut_list = [item for subitem in self.mut_pos for item in subitem]

def get_ngs_data_ncbi(param, folder_data_sra):
    # retrieve NGS data from NCBI, save them to data_sra folder
    if not len(glob.glob(folder_data_sra.as_posix() + "/*")) == 0: # check to see if file already downloaded
        print("Step0: sra already downloaded!")
    else:
        Path(folder_data_sra).mkdir(parents=True, exist_ok=True) # create folder to store sra
        folder_sra = Path(folder_data_sra.as_posix().replace(' ', '\ '))
        print("Step0: downloading sra from NCBI ...")
        for sra in param.ngs_data_ncbi:
            command = ' '.join(['fastq-dump', '--split-files', sra, '-O', '\'' + folder_sra + '\''])
            try:
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'))
            except:
                print("ERROR in Step0! Fail to download data from NCBI!!")
                exit()
        print("Step0: download succeeded!")

    # Rename sra file to match temperature: NrdF_28_R1.fastq etc
    print("Step0: renaming ...")
    file_pair = list(zip(param.ngs_data_ncbi, param.condition))
    try:
        for pair in file_pair:
            p = Path(folder_data_sra.joinpath(Path(pair[0]+"_1.fastq"))).as_posix().replace(' ', '\ ')
            target = Path(folder_data_sra.joinpath(Path(pair[1]+'_R1.fastq'))).as_posix().replace(' ', '\ ')
            shutil.copy(p, target)
            p = Path(folder_data_sra.joinpath(Path(pair[0]+"_2.fastq"))).as_posix().replace(' ', '\ ')
            target = Path(folder_data_sra.joinpath(Path(pair[1]+'_R2.fastq'))).as_posix().replace(' ', '\ ')
            shutil.copy(p, target)
    except:
        print("Step0: rename failed")
        exit()
    print('Step0: rename succeeded!')

def get_ngs_data_local(param, folder_data_sra):
    # retrieve NGS data from local storage
    if not len(glob.glob(folder_data_sra.as_posix() + "/*")) == 0: # check to see if file already downloaded
        print("Step0: sra already downloaded!")
    else:
        print("Step0: preparing ngs files from local storage ...")
        Path(folder_data_sra).mkdir(parents=True, exist_ok=True) # create folder to store sra
        folder_sra = Path(folder_data_sra).as_posix().replace(' ', '\ ')

    # Rename sra file to match temperature: NrdF_28_R1.fastq etc
    name_list = [Path(file) for file in param.ngs_data_local]
    file_pair = list(zip(name_list, param.condition))
    try:
        for i in range(len(param.condition)):
            prior_file  = name_list[i*2]
            target_file = folder_data_sra.joinpath(param.condition[i] + '_R1.fastq')
            try:
                shutil.copy(prior_file, target_file)
            except: pass
            prior_file  = name_list[i*2 +1]
            target_file = folder_data_sra.joinpath(param.condition[i] + '_R2.fastq')
            try:
                shutil.copy(prior_file, target_file)
            except: pass
    except:
        print("Step0: data preparation failed")
        exit()
    print('Step0: ngs data preparaiont succeeded!')

def prep_ref(param, folder_ref):
    # prepare ref seq by copying it into ref folder and renaming seq name to match file name
    filename = param.wtfile.name
    Path(folder_ref).mkdir(parents=True, exist_ok=True)
    wt_filename = folder_ref.joinpath(filename)
    _wtfile = SeqIO.parse(param.wtfile.as_posix().replace(' ', '\ '), 'fasta')
    tem_file = []
    for record in _wtfile:
        tem_file.append('>' + filename)
        tem_file.append(str(record.seq))

    Path(wt_filename).touch(exist_ok=True)

    file   = open(wt_filename, 'w+')
    for item in tem_file:
        file.write(item + '\n')
    file.close()

    param.wtfile = wt_filename

def fastq_merger(param, folder_data_sra, folder_merge):
    # Merge paired-end sequences into one single read
    if not len(glob.glob(folder_merge.as_posix() + "/*.fastq")) == 0:
        print("Step1: merged file already in merge folder!")
    else:
        Path(folder_merge).mkdir(parents=True, exist_ok=True) # create folder to store merged data
        print("Step1: merging ...")
        try:
            for cond in param.condition:
                infile1 = folder_data_sra.joinpath(Path(cond+'_R1.fastq'))
                infile2 = folder_data_sra.joinpath(Path(cond+'_R2.fastq'))
                outfile = folder_merge.joinpath(Path(cond+'.fastq'))
                command = ' '.join(['bbmerge.sh', 'in1=\''+infile1.as_posix().replace(' ', '\ ') + '\'',
                'in2=\''+infile2.as_posix().replace(' ', '\ ') + '\'', 'out=\''+outfile.as_posix().replace(' ', '\ ') + '\'', '-Xmx250m'])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
                # stderr is to silence bbmerge.sh
                # also note that for path to work with bbmerge.sh, path(with spaces) must be quoted and all spaces must be replaced with '\ '.
        except:
            print("Step1: merging failed!")
            exit()
            # Note that if bbmerge.sh is not installed, the error will not be set to 1 because all error infor are ouptut to os.devnull.
        print("Step1: merge succeeded!")

def first_mapper(param, folder_merge, folder_first_mapper):
    # Map merged reads onto reference genome
    try:
        command = ' '.join(['bwa index', param.wtfile.as_posix().replace(' ', '\ ')])
        subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    except:
        print("Step2: indexing ref failed!")
        exit()

    if not len(glob.glob(folder_first_mapper.as_posix() + "/*.csi")) == 0:
        print("Step2: first-round mapped files already in first_mapper folder!")
    else:
        Path(folder_first_mapper).mkdir(parents=True, exist_ok=True) # create folder to store first-round mapped files
        print("Step2: first-round mapping ...")
        try: # bwa mem
            for cond in param.condition:
                # mapping
                infile  = Path(folder_merge).joinpath(cond+'.fastq')
                outfile = Path(folder_first_mapper).joinpath(cond+'_bwa.sam')
                command = ' '.join(['bwa mem', param.wtfile.as_posix().replace(' ', '\ '), infile.as_posix().replace(' ', '\ '),
                 '>', outfile.as_posix().replace(' ', '\ ')])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # sorting
                infile  = Path(folder_first_mapper).joinpath(cond+'_bwa.sam')
                outfile = Path(folder_first_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools sort -m 250M', infile.as_posix().replace(' ', '\ '),
                '>', outfile.as_posix().replace(' ', '\ ')])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # indexing
                infile  = Path(folder_first_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools index -m 250M', infile.as_posix().replace(' ', '\ ')])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

        except:
            print("Step2: first-round mapping failed!")
            exit()
        print("Step2: first-round mapping succeeded!")

def quality_control_INDEL(param, workdir, folder_qc_INDEL):
    # Quality control step, filter out unqualified reads with INDELs
    if not len(glob.glob(folder_qc_INDEL.as_posix() + "/*.fastq")) == 0:
        print("Step3: qc INDEL already performed!")
    else:
        try:
            Path(folder_qc_INDEL).mkdir(parents=True, exist_ok=True) # create folder to store qc INDEL files.
            print("Step3: qc step: remove queries with INDELs and undetermined nucleotides ...")
            for cond in param.condition:
                #print(cond, test, "!!!!")
                samfile = pysam.AlignmentFile(workdir.joinpath('first_mapper/' + cond + "_bwa.bam"), 'rb')
                wt_name = param.wtfile.name
                itr = list(samfile.fetch(wt_name)) # Note the ref seq name must match its file name.
                                                   # itr can only be referenced once, so use list to unpack it first.

                target = [[item.cigartuples, item.query_name, item.query_sequence] for item in itr]
                df = pd.DataFrame(target, columns=['Cigar', 'Query_name', 'Query_seq'])

                df['Cigar_len'] = df.apply(lambda row:len(row['Cigar']), axis=1)
                df['Cigar_first'] = df.apply(lambda row:row['Cigar'][0][0], axis=1)
                df['Query_name']  = [item.query_name for item in itr]
                df['Query_unknown_nucleotide'] = df.apply(lambda row: all(x in ['A', 'T', 'C', 'G'] for x in row['Query_seq']), axis=1)

                read_count_clean_df     = df[(df['Cigar_len'] == 1) & (df['Cigar_first'] ==0) & (df['Query_unknown_nucleotide'] == True)]
                # 'Query_unknown_nucleotide' == True means that there is not undetermined nucleotide at that row.

                read_count_clean_count  = read_count_clean_df.shape[0]

                read_count_with_unknown_count = df[df['Query_unknown_nucleotide'] == False].count()[0]
                read_count_with_INDEL_count = df[~((df['Cigar_len'] == 1) & (df['Cigar_first'] == 0))].count()[0]

                print(" ".join(["Step3:", str(read_count_with_INDEL_count), 'sequences with INDELs;', str(read_count_with_unknown_count), 'sequences with undetermined nucleotides;', str(read_count_clean_count), 'valid sequences ...' ]))

                # Below is to split valid fastq files out and save them into folder_qc_INDEL folder.
                read_count_clean_query_name_file = workdir.joinpath('TemFolder/tem.txt')
                np.savetxt(read_count_clean_query_name_file,
                            read_count_clean_df['Query_name'].values, fmt='%s')
                raw_fastq_file = workdir.joinpath('merge/' + cond + ".fastq")
                out_fastq_file = folder_qc_INDEL.joinpath('qc_INDEL_' + cond + '.fastq')
                fastq_splitter_wrapper(raw_fastq_file, read_count_clean_query_name_file, out_fastq_file)
        except:
            print("Step3: qc step failed!")
        print("Step3: qc step: queries with INDELs and undetermined nucleotides removed!")

def quality_control_library_wrapper(param, folder_qc_INDEL, folder_qc_library):
    # Add in library sequences
    if not len(glob.glob(folder_qc_library.as_posix() + "/*.fastq")) == 0:
        print("Step4: library seqs already in qc_library folder!")
    else:
        folder_qc_library.mkdir(parents=True, exist_ok=True)
        #print(param.amplicon, param.mut_pos)
        #exit()
        assert (len(param.amplicon) == len(param.mut_pos)), "Amplicon # doesn't match mut_pos #!"

        library_seq = quality_control_library(param.wtfile, param.amplicon, param.mut_pos)
        for cond in param.condition:
            INDEL_list = open(Path(folder_qc_INDEL.joinpath('qc_INDEL_' + cond + '.fastq'))).readlines()
            library_list = INDEL_list + library_seq
            library_list = [x.strip() for x in library_list]
            library_seq_df = pd.DataFrame({'fastq': library_list})
            library_file_name = Path(folder_qc_library.joinpath(Path('qc_library_' + cond + '.fastq')))
            Path(library_file_name).touch(exist_ok=True)
            np.savetxt(library_file_name, library_seq_df['fastq'].values, fmt='%s')
        print("Step4: qc step: library seqs added!")

def second_mapper(param, folder_qc_library, folder_second_mapper):
    """
    Performs essentially the same operation on quality-controlled
    sequences (removal of INDELs and addition of library queries).
    """
    if not len(glob.glob(folder_second_mapper.as_posix() + "/*.csi")) == 0:
        print("Step5: second-round mapped files already in second_mapper folder!")
    else:
        Path(folder_second_mapper).mkdir(parents=True, exist_ok=True) # create folder to store second-round mapped files
        print("Step5: second-round mapping ...")
        try: # bwa mem
            for cond in param.condition:
                # mapping
                infile  = Path(folder_qc_library).joinpath('qc_library_' + cond+'.fastq')
                outfile = Path(folder_second_mapper).joinpath(cond+'_bwa.sam')
                command = ' '.join(['bwa mem', param.wtfile.as_posix().replace(' ', '\ '), infile.as_posix().replace(' ', '\ '),
                 '>', outfile.as_posix().replace(' ', '\ ')])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # sorting
                infile  = Path(folder_second_mapper).joinpath(cond+'_bwa.sam')
                outfile = Path(folder_second_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools sort -m 250M', infile.as_posix().replace(' ', '\ '),
                '>', outfile.as_posix().replace(' ', '\ ')])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))

                # indexing
                infile  = Path(folder_second_mapper).joinpath(cond+'_bwa.bam')
                command = ' '.join(['samtools index -m 250M', infile.as_posix().replace(' ', '\ ')])
                subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        except:
            print("Step5: second-round mapping failed!")
            exit()
        print("Step5: second-round mapping succeeded!")

def bam2enrich_wrapper(param, folder_second_mapper, folder_enrich2_input):
    """
    This function converts bam files into Enrich2 required input format
    and stores the output into folder enrich2_input folder.
    """
    if not len(glob.glob(folder_enrich2_input.as_posix() + "/*.fastq")) == 0:
        print("Step6: formatted Enrich2_input files already in enrich2_input folder!")
    else:
        if 1:
            print("Step6: converting bam files into Enrich2 input format ...")
            folder_enrich2_input.mkdir(parents=True, exist_ok=True)
            for cond in param.condition:
                bamfile = folder_second_mapper.joinpath(cond + '_bwa.bam')
                df = parse_bamfile(bamfile, param.wtfile)

                """
                df columns:
                           'Cigar'
                           'Query_name'
                           'Query_seq'
                           'Query_quality'
                           'Ref_pos' # list object, mapped query coordinates onto reference.
                           'Query_md'
                           'Cigar_len'
                           'Cigar_first'
                           'Query_name'
                           'Mut_pos' # list object, mapped query mutation site onto reference.
                           ‘Which_amplicon’ # stores the amplicon id, starting from 1
                """

                count_INDEL = df[~((df['Cigar_len'] == 1) & (df['Cigar_first'] == 0))].count()

                def which_amplicon(row, amplicon_range_list): # Assign the amplicon id for each of the query
                    """
                    The amplicon which overlaps the most to the given sequence is
                    the correct amplicon that the sequence comes from.
                    """
                    set1 = set(row['Ref_pos'])
                    max_value = 0
                    which_amplicon = 0
                    for i in range(len(amplicon_range_list)):
                        tem  = len(set1.intersection(
                                set(np.arange(amplicon_range_list[i][0],amplicon_range_list[i][1]))
                                ))/len(set1)
                        if tem > max_value:
                            max_value = tem
                            which_amplicon = i + 1
                    return(which_amplicon)

                def which_mut_site(row, mut_pos):
                    which_mut_site = '*' # default is *

                    set_mask = set()

                    for item in param.wt_mask:
                        try:
                            query_position = row['Ref_pos'].index(item[0])
                            if row['Query_seq'][query_position - 1] == item[1]: # why -1 would work?
                                set_mask.add(item[0])
                        except:
                            pass

                    for i in range(len(mut_pos)): # control for amplicon
                        for pos in mut_pos[i]:
                            start, end = pos, pos + 2
                            set_real    = set(row['Mut_pos'])
                            set_design  = set(np.arange(start, end+1))
                            if set_real.union(set_design) == set_design and row['Which_amplicon'] == i + 1:
                                if len(set_real) == 0:
                                    which_mut_site = 'wt'
                                    return which_mut_site
                                else:
                                    which_mut_site = start
                                    return which_mut_site

                            if len(set_mask) > 0:
                                if set_real.union(set_mask).union(set_design) == set_design.union(set_mask) and row['Which_amplicon'] == i + 1:
                                    if set_real == set_mask:
                                        which_mut_site = 'wt'
                                        return which_mut_site
                                    else:
                                        which_mut_site = start
                                        return which_mut_site

                    return(which_mut_site)

                df['Which_amplicon'] = df.apply(which_amplicon, axis=1, amplicon_range_list=param.amplicon)
                df['Which_mut_site'] = df.apply(which_mut_site, axis=1, mut_pos=param.mut_pos)

                wtseq = SeqIO.parse(param.wtfile.as_posix(), 'fasta')
                wtseq = next(wtseq).seq
                for i in range(len(param.mut_pos)):
                    for site in param.mut_pos[i]:
                        df_subset = df[(df["Which_mut_site"]==site) |
                                       (df["Which_mut_site"]=='wt') &
                                       (df["Which_amplicon"]==i+1)].copy()
                       # df_subset keeps only the sequences that have a single mutated aa at the specified site or wt.
                       # This, however, will also keep wt fragment that don't cover target mutation site. This can cause trouble.
                       # Use _fragment_filter to get rid of those fragments.
                        def _fragment_filter(row): # to remove fragments that don't cover mutation site
                            for pos in param.mut_pos[i]:
                                if not pos in row['Ref_pos']:
                                    return("Yes")
                            return("No")

                        df_subset['Fragment'] = df_subset.apply(lambda row: _fragment_filter(row), axis=1)
                        df_subset = df_subset[df_subset['Fragment'] == "No"]

                        df_subset["Which_mut_site2"] = site

                        df_subset['Ref_index'] = df_subset.apply(lambda row: row['Ref_pos'].index(row['Which_mut_site2']), axis=1)

                        df_subset['Mut_seq']   = df_subset.apply(lambda row: row['Query_seq'][row['Ref_index']-1:row['Ref_index']+2], axis=1)

                        df_subset['Mut_seq_quality'] = df_subset.apply(lambda row: row['Query_quality'][row['Ref_index']-1:row['Ref_index']+2], axis=1)

                        df_subset['Mut_seq_fastq']   = df_subset.apply(lambda row: '\n'.join(['@'+row['Query_name'], row['Mut_seq'], '+', ''.join([chr(x+33) for x in row['Mut_seq_quality']])]), axis=1)

                        outfilename = '_'.join([cond, str(site), str(site+2), str(wtseq[site-1:site+2])])
                        enrich_infile_name = Path(folder_enrich2_input.joinpath(outfilename + '.fastq'))
                        Path(enrich_infile_name).touch(exist_ok=True)
                        np.savetxt(enrich_infile_name, df_subset['Mut_seq_fastq'].values, fmt='%s')
                        # Should add one info to log.txt stating the number of reads that fall into each category.
        else:
            print("Step6: Enrich2 input preparation failed!")
            exit()
        print("Step6: Enrich2 input preparation succeeded!")

def enrich2_json_encoder_wrapper(param, folder_enrich2_json, folder_enrich2_input, folder_enrich2_output):
    # Build Enrich2 json files
    cond_input = param.condition[-1]
    cond_experiment = param.condition[0]
    condition_list = param.condition[:-1]

    if not len(glob.glob(folder_enrich2_json.as_posix() + "/*.json")) == 0:
        print("Step7: configuration files already in enrich2_json folder!")
    else:
        print("Step7: building configuration json files ...")
        folder_enrich2_json.mkdir(parents=True, exist_ok=True)
        json_commandfile_name = folder_enrich2_json.joinpath("json.sh")
        json_commandfile_name.touch(exist_ok=True)
        json_commandfile = open(json_commandfile_name, 'w+')
        json_commandfile.write("source activate py2-dms\n")
        try:
            for cond in condition_list:
                for mut in param.mut_list:
                    #file_prefix = '_'.join([cond, str(mut), str(mut+2)])
                    wt_code   = wt_codon(param.wtfile, mut)
                    #print("out: ", cond, cond_input, mut, wt_code)
                    res_json = enrich2_json_encoder(cond, cond_input, folder_enrich2_input, mut, wt_code)
                    jsonfilename = Path(folder_enrich2_json.joinpath("_".join([cond, str(mut), str(mut+2), str(wt_code)+".json"])))
                    jsonfilename.touch(exist_ok=True)
                    jsonfile = open(jsonfilename, 'w+')
                    for x in res_json: jsonfile.write(x)
                    jsonfile.close()

                    json_command = ' '.join(['enrich_cmd', jsonfilename.as_posix().replace(' ', '\ '),
                                        "ratios complete --no-plots --output-dir ",
                                        folder_enrich2_output.as_posix().replace(' ', '\ ') + "/" + cond + "_" + "_".join([str(mut), str(mut+2), str(wt_code)])])
                    json_commandfile.write(json_command + '\n')
            json_commandfile.close()
        except:
            print("Step7: configuration jsons creation failed!")
            exit()
        print("Step7: configuration jsons created in enrich2_json folder!")

def enrich2_wrapper(json_bashfile, folder_enrich2_output):
    # Perform Enrich2
    if not len(glob.glob(folder_enrich2_output.as_posix() + "/*")) == 0:
        print("Step8: Enrich2 already performed, results saved in enrich2_output folder!")
    else:
        print("Step8: performing Enrich2 ...")
        try:
            folder_enrich2_output.mkdir(parents=True, exist_ok=True)
            command = "bash " + json_bashfile.as_posix().replace(' ', '\ ')
            subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        except:
            print("Step8: Enrich2 failed to perform!")
            exit()
        print("Step8: Enrich2 performed, results saved in enrich2_output folder!")

def enrich2_hdf5_extractor_wrapper(param, folder_enrich2_output, _type='codon'):
    """
    For each target mutation site, extract and combine information from
    raw_count file and main_score file. Detailed rational, see notes in
    enrich2_hdf5_extractor function. The returned DataFrame contains the
    following columns:
    count0, codon, count1, score, SE, logratio, variance, mut_aa, wt_aa
    The type parameter specifies the extraction type: either codon or aa.
    codon will use default enrich2_hdf5_extractor function whereas aa uses
    a slightly modified enrich2_hdf5_extractor2 to retrieve aa data.
    use _type instead of type, because type is built-in keyword.

    To deal with false positives:
    1. If both read counts (count0 and count1) are below5: use nan to substitute both values;
    2. If both read counts are below 0.5% of its total counts:
    3. If count1 is below 0.5% of the total counts, assign 'deplete' to it.

    """

    cond_input = param.condition[-1]
    condition_list = param.condition[:-1]

    # Check type, if codon, save results into tsv folder, if aa, save to aa_tsv folder.
    if not _type == 'codon':
        folder_enrich2_output_tsv = folder_enrich2_output.joinpath('aa_tsv')
    else: folder_enrich2_output_tsv = folder_enrich2_output.joinpath('tsv')

    if not len(glob.glob(folder_enrich2_output_tsv.as_posix() + "/*")) == 0:
        print("Step9: extracted tsv/aa_tsv files already in enrich2_output/" +
              folder_enrich2_output_tsv.name + " folder!")
    else:
        if _type == 'codon':
            print("Step9: extracting codon data from Enrich2 output hdf5 files ...")
        else: print("Step9: extracting aa data by re-running Enrich2 via count mode ...")
        try:
            folder_enrich2_output_tsv.mkdir(parents=True, exist_ok=True)
            for cond in condition_list:
                for mut in param.mut_list:
                    wt_code = str(wt_codon(param.wtfile, mut))
                    file_folder = folder_enrich2_output.joinpath(
                    "_".join([cond, str(mut), str(mut+2), wt_code]))

                    raw_count_file1 = file_folder.joinpath(cond_input + "_lib.h5")
                    raw_count_file2 = file_folder.joinpath(cond + "_lib.h5")
                    score_file      = file_folder.joinpath(cond + "_sel.h5") # In fact, score_file.h5 also contain raw counts data

                    if _type == 'codon':
                        df = enrich2_hdf5_extractor(raw_count_file1, raw_count_file2, wt_code, score_file)
                    else:
                        df = enrich2_hdf5_extractor2(wt_code, score_file)

                    outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond, str(mut), str(mut+2), wt_code + '.tsv']))

                    df_original = np.concatenate((df.columns.values.reshape((1,-1)), df.values), axis = 0)
                    np.savetxt(outfilename, df_original, fmt="%s", delimiter='\t')
                    """
                    Above saves the original df into tsv folder.
                    """

                    """
                    If type == aa, combine the SY and WT into one single row using the formula:
                    """
                    if _type == 'aa': # combine SY and WT into one row (if specified analysis type is 'aa') by using count file to perform Enrich2 again
                        df_filtered = df.dropna() # Enrich2 will keep rows without NaNs, and use that to peform calculation.
                        c0 = df_filtered['count0'].sum()
                        c1 = df_filtered['count1'].sum()
                        c0_wt = df_filtered[df_filtered['wt_aa'] == 'WT']['count0'].sum()
                        c1_wt = df_filtered[df_filtered['wt_aa'] == 'WT']['count1'].sum()
                        res = enrich2_count_wrapper(c0_wt, c1_wt, c0, c1)
                        df = df_filtered[~(df_filtered['wt_aa']=='WT')]
                        res['count0'] = c0_wt
                        res['count1'] = c1_wt
                        res['wt_aa']  = df_filtered.iloc[-1, -2]
                        res['mut_aa'] = res['wt_aa']
                        df = df.append(res, ignore_index=True)

                    df2 = df[(df['count0']/df['count0'].sum() > 0.005) | (df['count1']/df['count1'].sum() > 0.005)].copy()
                    # Keep the rows where at least one count is no less than 0.005 of the total counts

                    df2['score'] = df2.apply(lambda row: row['score'] if row['count1']/df2['count1'].sum() > 0.005 else "deplete", axis=1)
                    # Assign "deplete" to rows with count1 below 0.005

                    df2['score'] = df2.apply(lambda row: row['score'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1) # if both counts are below 5, we tend to believe that the score is not reliable and assigna nan to it.

                    df2['SE'] = df2.apply(lambda row: row['SE'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)
                    # Adjust SE column according to score column: if score is nan, assign nan to SE too. Same for other columns: logratio and variance. Note, here if using lambda row: row['SE'] if not row['SE'] == np.nan, it won't work for unknown reason.
                    df2['logratio'] = df2.apply(lambda row: row['logratio'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)
                    df2['variance'] = df2.apply(lambda row: row['variance'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)

                    outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond, str(mut), str(mut+2), wt_code + '.preprocessed.tsv']))

                    df_preprocessed = np.concatenate((df2.columns.values.reshape((1,-1)), df2.values), axis = 0)
                    np.savetxt(outfilename, df_preprocessed, fmt="%s", delimiter='\t')
                    """
                    Above save preprocessed df into tsv folder.
                    """
        except:
            print("Step9: information extraction failed!")
            exit()
        print("Step9: information extracted from Enrich2 output!")

def tsv_plot_input(param, folder_plot_input, folder_enrich2_output):
    # Prepare plot input data
    condition_list = param.condition[:-1]

    if not len(glob.glob(folder_plot_input.as_posix() + "/*")) == 0:
        print("Step10: plot input files already in plot_input folder!")
    else:
        print("Step10: preparing plot input from extracted tsv files ...")
        folder_plot_input.mkdir(parents=True, exist_ok=True)
        _CodonList = CodonList + ['pos']
        #AAList = [key for key in AA2CODON]
        _AAList    = AAList + ['pos']

        for _type in ('tsv', 'aa_tsv'):
            folder_enrich2_output_tsv = folder_enrich2_output.joinpath(_type)

            if _type == 'tsv':
                index = 'codon'
                ColumnList = _CodonList
            else:
                index = 'mut_aa'
                ColumnList = _AAList
            try: # change to try
                for cond in condition_list:
                    df = pd.DataFrame(columns=ColumnList)
                    df_se = df.copy()
                    for mut in param.mut_list:
                        row = pd.Series(index=df.columns)
                        tsvfile = "_".join([cond, str(mut), str(mut+2), str(wt_codon(param.wtfile, mut))+ ".preprocessed.tsv"])
                        df_tsv = pd.read_csv(folder_enrich2_output_tsv.joinpath(tsvfile), delimiter = '\t')

                        df_tsv.set_index([index], inplace=True)
                        #df_tsv['score_se'] = df_tsv.apply(lambda row: (row['score'], row['SE']), axis=1)
                        row_se = row.copy()

                        row.update(df_tsv['score']) # only keep Enrich2 score
                        row['pos'] = int((mut+2)/3)
                        df = df.append(row, ignore_index=True)

                        row_se.update(df_tsv['SE'])
                        row_se['pos'] = int((mut+2)/3)
                        df_se = df_se.append(row_se, ignore_index=True)

                    if _type == "aa_tsv":
                        filler = '_aa'
                        df = df.rename(columns={'Ter': 'Stop'})
                        df_se = df_se.rename(columns={'Ter': 'Stop'})
                    else: filler = ''

                    savefile = folder_plot_input.joinpath(cond+filler+".tsv")
                    savefile.touch(exist_ok=True)

                    df2array = np.concatenate((df.columns.values.reshape((1,-1)), df.values), axis = 0)
                    np.savetxt(savefile, df2array, fmt='%s', delimiter='\t')

                    savefile = folder_plot_input.joinpath(cond+filler+".se.tsv")
                    savefile.touch(exist_ok=True)

                    df2array = np.concatenate((df_se.columns.values.reshape((1,-1)), df_se.values), axis = 0)
                    np.savetxt(savefile, df2array, fmt='%s', delimiter='\t')
            except:
                print("Step10: tsv to plot_input failed!")
                exit()
        print("Step10: tsv to plot_input step succeeded!")

def tsv_plot_output_wrapper(param, folder_plot_input, folder_plot_output_codon, scale = 'max'):
    """
    Generate multiple plots using tsv_plot_output function.
    Output will be saved in plot_output folder.
    Multiple plots will be created:
    1. raw: contains every information: all codons (64 for codon and 21 for aa), SE
    2. simple1: all-NaN columns removed, SE
    3. simple2: all-NaN columns removed, SE removed
    """
    if not len(glob.glob(folder_plot_output_codon.as_posix() + "/*")) == 0:
        print("Step11: plot output files already in plot_output/codon folder!")
    else:
        print("Step11: preparing plot output files (codon version) ...")
        folder_plot_output_codon.mkdir(parents=True, exist_ok=True)
        try:
            # Plot1:
            for cond in param.condition[:-1]:
                df1      = pd.read_csv(folder_plot_input.joinpath(cond+'.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'.se.tsv'), delimiter='\t')
                outfile1 = folder_plot_output_codon.joinpath(cond + '_raw.pdf')
                tsv_plot_output(param.wtfile, cond, df1, df_se=df_se, outfilename=outfile1, scale = scale)
            # Plot2:
            for cond in param.condition[:-1]:
                df2      = pd.read_csv(folder_plot_input.joinpath(cond+'.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'.se.tsv'), delimiter='\t')
                df2.replace(to_replace='deplete', value=1000, inplace=True)
                df2 = df2.apply(lambda col: pd.to_numeric(col), axis=0)
                for col in df2.columns:
                    if sum(np.isnan(df2[col])) == len(df2.index):
                        df2.drop(columns=col, inplace=True)
                            # Note tat np.nan != np.nan
                outfile2 = folder_plot_output_codon.joinpath(cond + '_simple1.pdf')
                tsv_plot_output(param.wtfile, cond, df2, df_se=df_se, outfilename=outfile2, version=2, scale = scale)
            # Plot3:
                outfile3 = folder_plot_output_codon.joinpath(cond + '_simple2.pdf')
                tsv_plot_output(param.wtfile, cond, df2, outfilename=outfile3, version=2, scale = scale)
        except:
            print("Step11: plots (codon version) creation failed!")
            exit()
        print("Step11: plots (codon version) created!")

def tsv_plot_output_aa_wrapper(param, folder_plot_input, folder_plot_output_aa, scale = 'max'):
    """
    Same function as tsv_plot_output_wrapper but optimized for aa mode.
    There will
    """
    if not len(glob.glob(folder_plot_output_aa.as_posix() + "/*")) == 0:
        print("Step11: plot output files already in plot_output/aa folder!")
    else:
        print("Step11: preparing plot output files (aa version) ...")
        folder_plot_output_aa.mkdir(parents=True, exist_ok=True)
        try:
            # Plot1:
            for cond in param.condition[:-1]:
                df1      = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.se.tsv'), delimiter='\t')
                df1.drop(columns = ['???'], inplace=True) # After adding undetermined nucleotide removal step the '???' column is no longer useful.
                df_se.drop(columns = ['???'], inplace=True)

                outfile1 = folder_plot_output_aa.joinpath(cond + '_raw.pdf')
                tsv_plot_output_aa(param.wtfile, cond, df1, df_se=df_se, outfilename=outfile1, scale = scale)
            # Plot2:
            for cond in param.condition[:-1]:
                df2      = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.se.tsv'), delimiter='\t')
                df2.replace(to_replace='deplete', value=1000, inplace=True)
                df2 = df2.apply(lambda col: pd.to_numeric(col), axis=0)
                for col in df2.columns:
                    if sum(np.isnan(df2[col])) == len(df2.index):
                        df2.drop(columns=col, inplace=True)
                            # Note that np.nan != np.nan
                outfile2 = folder_plot_output_aa.joinpath(cond + '_simple1.pdf')
                tsv_plot_output_aa(param.wtfile, cond, df2, df_se=df_se, outfilename=outfile2, scale = scale)
            # Plot3:
                outfile3 = folder_plot_output_aa.joinpath(cond + '_simple2.pdf')
                tsv_plot_output_aa(param.wtfile, cond, df2, outfilename=outfile3, scale = scale)
        except:
            print("Step11: plots (aa version) creation failed!")
            exit()
        print("Step11: plots (aa version) created!")

"""
Below are adapted functions for double mutational pipeline.
"""
def quality_control_library_double_wrapper(param, folder_qc_INDEL, folder_qc_library):
    # Add in library sequence for double mutational pipeline
    if not len(glob.glob(folder_qc_library.as_posix() + "/*.fastq")) == 0:
        print("Step4: library seqs already in qc_library folder!")
    else:
        folder_qc_library.mkdir(parents=True, exist_ok=True)
        assert (len(param.amplicon) == len(param.mut_pos)), "Amplicon # doesn't match mut_pos #!"
        library_seq = quality_control_library_double(param.wtfile, param.amplicon, param.mut_pos)
        for cond in param.condition:
            INDEL_list = open(Path(folder_qc_INDEL.joinpath('qc_INDEL_' + cond + '.fastq'))).readlines()
            library_list = INDEL_list + library_seq
            library_list = [x.strip() for x in library_list]
            library_seq_df = pd.DataFrame({'fastq': library_list})
            library_file_name = Path(folder_qc_library.joinpath(Path('qc_library_' + cond + '.fastq')))
            Path(library_file_name).touch(exist_ok=True)
            np.savetxt(library_file_name, library_seq_df['fastq'].values, fmt='%s')
        print("Step4: qc step: library seqs added!")

def bam2enrich_double_wrapper(param, folder_second_mapper, folder_enrich2_input):
    """
    This function converts bam files into Enrich2 required input format
    and stores the output into folder enrich2_input folder. Optimized for double mutant scheme. Via count mode.
    """
    if not len(glob.glob(folder_enrich2_input.as_posix() + "/*.txt")) == 0:
        print("Step6: formatted Enrich2_input files already in enrich2_input folder!")
    else:
        try:
            print("Step6: converting bam files into Enrich2 input format ...")
            folder_enrich2_input.mkdir(parents=True, exist_ok=True)
            for cond in param.condition:
                bamfile = folder_second_mapper.joinpath(cond + '_bwa.bam')
                df = parse_bamfile(bamfile, param.wtfile)
                """
                df columns:
                           'Cigar'
                           'Query_name'
                           'Query_seq'
                           'Query_quality'
                           'Ref_pos' # list object, mapped query coordinates onto reference.
                           'Query_md'
                           'Cigar_len'
                           'Cigar_first'
                           'Query_name'
                           'Mut_pos' # list object, mapped query mutation site onto reference.
                           ‘Which_amplicon’ # stores the amplicon id, starting from 1
                           'Mut1_seq' # the sequence of mutation 1 on each query, used for splitting into 32 classes
                """

                count_INDEL = df[~((df['Cigar_len'] == 1) & (df['Cigar_first'] == 0))].count() # double check to see if any INDEL exists or not

                def which_amplicon(row, amplicon_range_list): # Assign the amplicon id for each of the query
                    """
                    The amplicon which overlaps the most to the given sequence is
                    the correct amplicon that the sequence comes from.
                    """
                    set1 = set(row['Ref_pos'])
                    max_value = 0
                    which_amplicon = 0
                    for i in range(len(amplicon_range_list)):
                        tem  = len(set1.intersection(
                                set(np.arange(amplicon_range_list[i][0],amplicon_range_list[i][1]))
                                ))/len(set1)
                        if tem > max_value:
                            max_value = tem
                            which_amplicon = i + 1
                    return(which_amplicon)

                def which_mut_site(row, mut_pos):
                    which_mut_site = '*' # default is *, if no more extra mutations, it will be set to the start position of designed site, else remains as *
                    set_mask = set()
                    for item in param.wt_mask:
                        try:
                            query_position = row['Ref_pos'].index(item[0])
                            if row['Query_seq'][query_position - 1] == item[1]: # why -1 would work?
                                set_mask.add(item[0])
                        except:
                            pass

                    for i in range(len(mut_pos)): # control for amplicon
                        set_design = set([mut_pos[0][0], mut_pos[0][0]+1, mut_pos[0][0]+2, mut_pos[0][1], mut_pos[0][1]+1, mut_pos[0][1]+2])
                        set_real   = set(row['Mut_pos'])

                        if set_real.union(set_design) == set_design and row['Which_amplicon'] == i + 1:
                            if len(set_real) == 0:
                                which_mut_site = 'wt'
                                return which_mut_site
                            else:
                                which_mut_site = set_real
                                return which_mut_site

                        if len(set_mask) > 0:
                            if set_real.union(set_mask).union(set_design) == set_design.union(set_mask) and row['Which_amplicon'] == i + 1:
                                if set_real == set_mask:
                                    which_mut_site = 'wt'
                                    return which_mut_site
                                else:
                                    which_mut_site = set_real
                                    return which_mut_site

                    return(which_mut_site)

                df['Which_amplicon'] = df.apply(which_amplicon, axis=1, amplicon_range_list=param.amplicon)
                df['Which_mut_site'] = df.apply(which_mut_site, axis=1, mut_pos=param.mut_pos)
                #print(df.shape) # To check how many sequences are in df
                df = df[df['Which_mut_site'] != '*']
                # Keep sequences that have mutations only at designed sites or wt
                # This, however, will also keep wt fragment that don't cover target mutation site. This can cause trouble.
                # Use _fragment_filter to get rid of those fragments.

                def _fragment_filter(row): # to remove fragments that don't cover mutation site
                    for pos in param.mut_pos[0]:
                        if not pos in row['Ref_pos']:
                            return("Yes")
                    return("No")
                # df_subset keeps only the sequences that have a single mutated aa at the specified site or wt.

                df['Fragment'] = df.apply(lambda row: _fragment_filter(row), axis=1)
                df = df[df['Fragment'] == "No"]
                #print(df.shape) # To check how many sequences left after filtering seq with unexpected mutations

                def site_seq(row, site):
                    """
                    Find and return the sequence for the target mutation site.
                    """
                    site = row['Ref_pos'].index(site-1)
                    return(row['Query_seq'][site: site+3])

                df['Mut1_seq'] = df.apply(site_seq, axis=1, site = param.mut_pos[0][0])
                df['Mut2_seq'] = df.apply(site_seq, axis=1, site = param.mut_pos[0][1])

                def site_aa(row, site):
                    """
                    Find the corresponding aa for the target mutation site.
                    """
                    codon = row[site]
                    aa = CODON2AA[codon][0]
                    return(aa)

                df['Mut1_aa']  = df.apply(site_aa, axis=1, site = 'Mut1_seq')
                df['Mut2_aa']  = df.apply(site_aa, axis=1, site = 'Mut2_seq')

                wtseq = SeqIO.parse(param.wtfile.as_posix(), 'fasta')
                wtseq = next(wtseq).seq

                Mut_pool = df['Mut1_seq'].unique()
                #print(Mut_pool) # Checkpoint for unexpected codon

                """
                Counting on codon basis.
                """
                outfilename = '_'.join([cond, 'codon', str(wtseq[param.mut_pos[0][0]-1:param.mut_pos[0][0]+2]), str(wtseq[param.mut_pos[0][1]-1:param.mut_pos[0][1]+2])])
                enrich_infile_name = Path(folder_enrich2_input.joinpath(outfilename + '.txt'))
                outfile = open(enrich_infile_name, 'w+')
                outfile.write("\t" + "count" + "\n")

                for site1 in CodonList:
                    for site2 in CodonList:
                        count_number = df[(df["Mut1_seq"]==site1) & (df['Mut2_seq']==site2)].shape[0]
                        outfile.write('_'.join([site1, site2]) + "\t" + str(count_number) + "\n")

                outfile.close()

                """
                Counting on aa basis.
                """
                outfilename = '_'.join([cond, 'aa', str(wtseq[param.mut_pos[0][0]-1:param.mut_pos[0][0]+2]), str(wtseq[param.mut_pos[0][1]-1:param.mut_pos[0][1]+2])])
                enrich_infile_name = Path(folder_enrich2_input.joinpath(outfilename + '.txt'))
                outfile = open(enrich_infile_name, 'w+')
                outfile.write("\t" + "count" + "\n")

                for site1 in CODON_GROUPS:
                    site1_aa = site1[0]
                    for site2 in CODON_GROUPS:
                        site2_aa = site2[0]
                        count_number = df[(df['Mut1_aa'] == site1_aa) & (df['Mut2_aa'] == site2_aa)].shape[0]
                        outfile.write('_'.join([site1_aa, site2_aa]) + '\t' + str(count_number) + "\n")

                outfile.close()
        except:
            print("Step6: Enrich2 input preparation failded!")
            exit()
        print("Step6: Enrich2 input preparation succeeded!")

def enrich2_json_encoder_double_wrapper(param, folder_enrich2_json, folder_enrich2_input, folder_enrich2_output):
    # Build Enrich2 json files for double mutational pipeline
    cond_input = param.condition[-1]

    if not len(glob.glob(folder_enrich2_json.as_posix() + "/*.json")) == 0:
        print("Step7: configuration files already in enrich2_json folder!")
    else:
        print("Step7: building configuration json files ...")
        folder_enrich2_json.mkdir(parents=True, exist_ok=True)
        json_commandfile_name = folder_enrich2_json.joinpath("json.sh")
        json_commandfile_name.touch(exist_ok=True)
        json_commandfile = open(json_commandfile_name, 'w+')
        json_commandfile.write("source activate py2-dms\n")
        try:
            mut_1 = str(wt_codon(param.wtfile, param.mut_pos[0][0]))
            mut_2 = str(wt_codon(param.wtfile, param.mut_pos[0][1]))

            for _type in ('codon', 'aa'):
                count_file = []
                for cond in param.condition:
                    filename = '_'.join([cond, _type, mut_1, mut_2]) + '.txt'
                    count_file.append(folder_enrich2_input.joinpath(filename))

                res_json = enrich2_json_encoder_count([count_file[1].as_posix(),
                                                       count_file[0].as_posix()],
                                                       folder_enrich2_json.as_posix())

                json_outfile = folder_enrich2_json.joinpath('_'.join([cond, _type, mut_1, mut_2]) + '.json')
                json_outfile.touch(exist_ok=True)
                jsonfile = open(json_outfile, 'w+')
                for x in res_json: jsonfile.write(x)
                jsonfile.close()

                json_command = ' '.join(['enrich_cmd', json_outfile.as_posix().replace(' ', '\ '),
                                    "ratios complete --no-plots --output-dir ",
                                    folder_enrich2_output.as_posix().replace(' ', '\ ') + "/" + "_".join([cond, _type, mut_1, mut_2])])

                json_commandfile.write(json_command + '\n')
            json_commandfile.close()

            command = "bash " + folder_enrich2_json.joinpath('json.sh').as_posix().replace(' ', '\ ')
            subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        except:
            print("Step7: configuration jsons creation failed!")
            exit()
        print("Step7: configuration jsons created in enrich2_json folder! Enrich2 performed!")

def enrich2_tsv_extractor_double(param, folder_enrich2_output):
    """
    Extract count data for codon and aa mode.
    """
    cond_input     = param.condition[-1]
    cond_experiment= param.condition[0]
    condition_list = param.condition[:-1]
    mut_1 = str(wt_codon(param.wtfile, param.mut_pos[0][0]))
    mut_2 = str(wt_codon(param.wtfile, param.mut_pos[0][1]))
    wt_code = str(wt_codon(param.wtfile, param.mut_pos[0][1]))

    """
    Count codon mode.
    """
    _path_count = '_'.join([cond_input, 'codon', mut_1, mut_2]) + "/tsv/Enrich2_count_wrapper_sel/main_variants_counts_unfiltered.tsv"
    _path_score = '_'.join([cond_input, 'codon', mut_1, mut_2]) + "/tsv/Enrich2_count_wrapper_sel/main_variants_scores.tsv"
    folder_enrich2_output_tsv = folder_enrich2_output.joinpath('tsv')
    folder_enrich2_output_tsv.mkdir(parents=True, exist_ok=True)

    if not len(glob.glob(folder_enrich2_output_tsv.as_posix() + "/*")) == 0:
        print("Step8: extracted tsv/aa_tsv files already in enrich2_output/" +
              folder_enrich2_output_tsv.name + " folder!")
    else:
        print("Step8: extracting info from enrich2_output files (codon)")
        wt_aa = CODON2AA[wt_code][0]
        count_file = folder_enrich2_output.joinpath(_path_count)
        score_file = folder_enrich2_output.joinpath(_path_score)

        df_count   = pd.read_csv(count_file, delimiter='\t')
        df_score   = pd.read_csv(score_file, delimiter='\t')
        df_combine = df_count.set_index('Unnamed: 0').join(df_score.set_index('Unnamed: 0'))
        df_combine = df_combine.reset_index()
        df_combine = df_combine.rename(columns={"Unnamed: 0": "Mut_combine"})

        def findMut1(row):
            return(row["Mut_combine"].split("_")[0])
        def findMut2(row):
            return(row["Mut_combine"].split("_")[1])

        df_combine['Mut1']   = df_combine.apply(lambda row: findMut1(row), axis=1)
        df_combine['codon']  = df_combine.apply(lambda row: findMut2(row), axis=1)
        df_combine['mut_aa'] = df_combine.apply(lambda row: CODON2AA[row['codon']][0], axis=1)
        df_combine['wt_aa']  = wt_aa
        df_combine = df_combine.rename(columns={'c_0': 'count0', 'c_1': 'count1'})

        for mut in range(1,65):
            outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond_experiment, str(mut), str(mut+2), wt_code + '.tsv']))
            df_subset = df_combine[(df_combine['Mut1'] == CodonList[mut-1])]
            df_subset = df_subset.drop(columns = ['Mut_combine', 'Mut1'])

            df_original = np.concatenate((df_subset.columns.values.reshape((1,-1)), df_subset.values), axis = 0)
            np.savetxt(outfilename, df_original, fmt="%s", delimiter='\t')

            df2 = df_subset[(df_subset['count0']/(df_subset['count0'].sum() + 1) > 0.005) | (df_subset['count1']/df_subset['count1'].sum() > 0.005)].copy()
            # Keep the rows where at least one count is no less than 0.005 of the total counts

            if df2.empty:
                continue
                # Skip empty df2

            #print(df2.head())
            df2['score'] = df2.apply(lambda row: row['score'] if row['count1']/(df2['count1'].sum() + 1) > 0.005 else "deplete", axis=1)
            # Assign "deplete" to rows with count1 below 0.005

            df2['score'] = df2.apply(lambda row: row['score'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1) # if both counts are below 5, we tend to believe that the score is not reliable and assigna nan to it.

            df2['SE'] = df2.apply(lambda row: row['SE'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1) # Adjust SE column according to score column: if score is nan, assign nan to SE too. Same for other columns: logratio and variance
            # Note, here if using lambda row: row['SE'] if not row['SE'] == np.nan, it won't work for unknown reason.
            df2['logratio'] = df2.apply(lambda row: row['logratio'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)
            df2['variance'] = df2.apply(lambda row: row['variance'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)

            outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond_experiment, str(mut), str(mut+2), wt_code + '.preprocessed.tsv']))

            df_preprocessed = np.concatenate((df2.columns.values.reshape((1,-1)), df2.values), axis = 0)
            np.savetxt(outfilename, df_preprocessed, fmt="%s", delimiter='\t')

    """
    Count aa mode.
    """
    _path_count = '_'.join([cond_input, 'aa', mut_1, mut_2]) + "/tsv/Enrich2_count_wrapper_sel/main_variants_counts_unfiltered.tsv"
    _path_score = '_'.join([cond_input, 'aa', mut_1, mut_2]) + "/tsv/Enrich2_count_wrapper_sel/main_variants_scores.tsv"
    folder_enrich2_output_tsv = folder_enrich2_output.joinpath('aa_tsv')
    folder_enrich2_output_tsv.mkdir(parents=True, exist_ok=True)

    if not len(glob.glob(folder_enrich2_output_tsv.as_posix() + "/*")) == 0:
        print("Step8: extracted tsv/aa_tsv files already in enrich2_output/" +
              folder_enrich2_output_tsv.name + " folder!")
    else:
        print("Step8: extracting info from enrich2_output files (aa)")
        try:
            wt_aa = CODON2AA[wt_code][0]
            count_file = folder_enrich2_output.joinpath(_path_count)
            score_file = folder_enrich2_output.joinpath(_path_score)
            df_count   = pd.read_csv(count_file, delimiter='\t')
            df_score   = pd.read_csv(score_file, delimiter='\t')

            df_combine = df_count.set_index('Unnamed: 0').join(df_score.set_index('Unnamed: 0'))
            df_combine = df_combine.reset_index()
            df_combine = df_combine.rename(columns={"Unnamed: 0": "Mut_combine"})

            def findMut1(row):
                return(row['Mut_combine'].split("_")[0])

            def findMut2(row):
                return(row['Mut_combine'].split("_")[1])

            df_combine['mut1_aa']   = df_combine.apply(lambda row: findMut1(row), axis=1)
            df_combine['mut2_aa']  = df_combine.apply(lambda row: findMut2(row), axis=1)

            df_combine['wt_aa']  = wt_aa

            df_combine = df_combine.rename(columns={'c_0': 'count0', 'c_1': 'count1'})

            for group in CODON_GROUPS: # combine codons that translate to the same aa
                aa, start, end = group[0], group[1], group[2]

                outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond_experiment, str(start), str(end), wt_code + '.tsv']))
                df_subset = df_combine[(df_combine['mut1_aa'] == aa)]
                df_subset = df_subset.drop(columns = ['Mut_combine', 'mut1_aa'])
                df_subset = df_subset.rename(columns={'mut2_aa': 'mut_aa'})

                # Save the original aa based Enrich2 scores:
                df_output = np.concatenate((df_subset.columns.values.reshape((1,-1)), df_subset.values), axis = 0)
                np.savetxt(outfilename, df_output, fmt="%s", delimiter='\t')

                df2 = df_subset[(df_subset['count0']/(df_subset['count0'].sum() + 1) > 0.005) | (df_subset['count1']/(df_subset['count1'].sum() + 1) > 0.005)].copy()
                # Keep the rows where at least one count is no less than 0.005 of the total counts

                if df2.empty:
                    continue
                    # Skip empty df2

                df2['score'] = df2.apply(lambda row: row['score'] if row['count1']/df2['count1'].sum() > 0.005 else "deplete", axis=1)
                # Assign "deplete" to rows with count1 below 0.005

                df2['score'] = df2.apply(lambda row: row['score'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1) # if both counts are below 5, we tend to believe that the score is not reliable and assigna nan to it.

                df2['SE'] = df2.apply(lambda row: row['SE'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1) # Adjust SE column according to score column: if score is nan, assign nan to SE too. Same for other columns: logratio and variance
                # Note, here if using lambda row: row['SE'] if not row['SE'] == np.nan, it won't work for unknown reason.
                df2['logratio'] = df2.apply(lambda row: row['logratio'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)
                df2['variance'] = df2.apply(lambda row: row['variance'] if ((row['count1'] > 5) or (row['count0'] > 5)) else np.nan, axis=1)

                outfilename = folder_enrich2_output_tsv.joinpath("_".join([cond_experiment, str(start), str(end), wt_code + '.preprocessed.tsv']))

                df_preprocessed = np.concatenate((df2.columns.values.reshape((1,-1)), df2.values), axis = 0)
                np.savetxt(outfilename, df_preprocessed, fmt="%s", delimiter='\t')
        except:
            print("Step8: information extraction failed!")
            exit()
        print("Step8: information extracted from Enrich2 output!")

def tsv_plot_double_input(param, folder_plot_input, folder_enrich2_output):
    # Prepare plot input, optimized for double mutation pipeline
    condition_list = param.condition[:-1]
    if not len(glob.glob(folder_plot_input.as_posix() + "/*")) == 0:
        print("Step9: plot input files already in plot_input folder!")
    else:
        print("Step9: preparing plot input from extracted tsv files ...")
        folder_plot_input.mkdir(parents=True, exist_ok=True)
        _CodonList = CodonList + ['pos']
        _AAList    = AAList + ['pos']

        for _type in ('tsv', 'aa_tsv'):
            folder_enrich2_output_tsv = folder_enrich2_output.joinpath(_type)

            if _type == 'tsv':
                index = 'codon'
                ColumnList = _CodonList
            else:
                index = 'mut_aa'
                ColumnList = _AAList
            try: # change to try
                for cond in condition_list:
                    df = pd.DataFrame(columns=ColumnList)
                    df_se = df.copy()
                    if _type == 'tsv': # Deal with codon mode:
                        for mut in range(1,65):
                            row = pd.Series(index=df.columns)
                            tsvfile = "_".join([cond, str(mut), str(mut+2), str(wt_codon(param.wtfile, param.mut_pos[0][1]))+ ".preprocessed.tsv"])

                            if not Path(folder_enrich2_output_tsv.joinpath(tsvfile)).is_file(): # Skip empty files
                                continue
                            df_tsv = pd.read_csv(folder_enrich2_output_tsv.joinpath(tsvfile), delimiter = '\t')
                            df_tsv.set_index([index], inplace=True)
                            #df_tsv['score_se'] = df_tsv.apply(lambda row: (row['score'], row['SE']), axis=1)
                            row_se = row.copy()
                            row.update(df_tsv['score']) # only keep Enrich2 score
                            row['pos'] = int(mut)
                            df = df.append(row, ignore_index=True)
                            row_se.update(df_tsv['SE'])
                            row_se['pos'] = int(mut)
                            df_se = df_se.append(row_se, ignore_index=True)

                    elif _type == 'aa_tsv': # Deal with aa mode:
                        df = df.rename(columns={'Ter': 'Stop'})
                        df_se = df.rename(columns={'Ter': 'Stop'})
                        for group in CODON_GROUPS:
                            aa, start, end = group[0], group[1], group[2]
                            row = pd.Series(index=df.columns)
                            tsvfile = "_".join([cond, str(start), str(end), str(wt_codon(param.wtfile, param.mut_pos[0][1]))+ ".preprocessed.tsv"])

                            if not Path(folder_enrich2_output_tsv.joinpath(tsvfile)).is_file(): # Skip empty files
                                continue
                            df_tsv = pd.read_csv(folder_enrich2_output_tsv.joinpath(tsvfile), delimiter = '\t')
                            df_tsv.set_index([index], inplace=True)
                            #df_tsv['score_se'] = df_tsv.apply(lambda row: (row['score'], row['SE']), axis=1)
                            row_se = row.copy()
                            row.update(df_tsv['score']) # only keep Enrich2 score
                            row['pos'] = int(start)
                            df = df.append(row, ignore_index=True)
                            row_se.update(df_tsv['SE'])
                            row_se['pos'] = int(start)
                            df_se = df_se.append(row_se, ignore_index=True)

                    if _type == "aa_tsv":
                        filler = '_aa'
                        df = df.rename(columns={'Ter': 'Stop'})
                        df_se = df_se.rename(columns={'Ter': 'Stop'})
                    else: filler = ''

                    savefile = folder_plot_input.joinpath(cond+filler+".tsv")

                    savefile.touch(exist_ok=True)

                    df2array = np.concatenate((df.columns.values.reshape((1,-1)), df.values), axis = 0)
                    np.savetxt(savefile, df2array, fmt='%s', delimiter='\t')

                    savefile = folder_plot_input.joinpath(cond+filler+".se.tsv")
                    savefile.touch(exist_ok=True)

                    df2array = np.concatenate((df_se.columns.values.reshape((1,-1)), df_se.values), axis = 0)
                    np.savetxt(savefile, df2array, fmt='%s', delimiter='\t')
            except:
                print("Step9: tsv to plot_input failed!")
                exit()
        print("Step9: tsv to plot_input step succeeded!")

def tsv_plot_output_double_wrapper(param, folder_plot_input, folder_plot_output_codon, scale = 'max'):
    # Same as tsv_plot_output_wrapper but optimized for double mutation scheme.
    if not len(glob.glob(folder_plot_output_codon.as_posix() + "/*")) == 0:
        print("Step10: plot output files already in plot_output/codon folder!")
    else:
        print("Step10: preparing plot output files (codon version) ...")
        folder_plot_output_codon.mkdir(parents=True, exist_ok=True)
        try:
            # Plot1:
            for cond in param.condition[:-1]:
                df1      = pd.read_csv(folder_plot_input.joinpath(cond+'.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'.se.tsv'), delimiter='\t')
                outfile1 = folder_plot_output_codon.joinpath(cond + '_raw.pdf')

                # Add in missing rows to make it full 64 rows: data
                row_insert_template = pd.Series([np.nan for i in range(1,66)], index=df1.columns.values)
                for i in range(1,65):
                    if not i in df1['pos'].values:
                        row_insert_template['pos'] = i
                        df1 = df1.append(row_insert_template, ignore_index=True)
                df1.sort_values('pos', inplace=True, ascending=False)
                # Add in missing row rows to make it full 64: se
                row_insert_template = pd.Series([np.nan for i in range(1,66)], index=df_se.columns.values)
                for i in range(1,65):
                    if not i in df_se['pos'].values:
                        row_insert_template['pos'] = i
                        df_se = df_se.append(row_insert_template, ignore_index=True)
                df_se.sort_values('pos', inplace=True, ascending=False)

                def findCodon(row):
                    res = CodonList[int(row['pos']) - 1]
                    return(res)

                df1['Codon'] = df1.apply(lambda row: findCodon(row), axis=1)
                tsv_plot_output_double(param.wtfile, param.mut_pos, cond, df1, df_se=df_se, outfilename=outfile1, scale = scale)

            # Plot2:
            for cond in param.condition[:-1]:
                df2      = pd.read_csv(folder_plot_input.joinpath(cond+'.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'.se.tsv'), delimiter='\t')
                # Add in missing rows to make it full 64 rows: data
                row_insert_template = pd.Series([np.nan for i in range(1,66)], index=df2.columns.values)
                for i in range(1,65):
                    if not i in df2['pos'].values:
                        row_insert_template['pos'] = i
                        df2 = df2.append(row_insert_template, ignore_index=True)
                df2.sort_values('pos', inplace=True, ascending=False)

                # Add in missing rows to make it full 64: se
                row_insert_template = pd.Series([np.nan for i in range(1,66)], index=df_se.columns.values)
                for i in range(1,65):
                    if not i in df_se['pos'].values:
                        row_insert_template['pos'] = i
                        df_se = df_se.append(row_insert_template, ignore_index=True)
                df_se.sort_values('pos', inplace=True, ascending=False)

                # Drop columns with all NAN:
                df2.replace(to_replace='deplete', value=1000, inplace=True)
                df2 = df2.apply(lambda col: pd.to_numeric(col), axis=0)
                df2.dropna(axis = 1, how='all', inplace=True)
                # Drop rows with all NAN:
                cutoff = df2.shape[0]
                colNameSet = df2.columns.values[:-1] # ignore "pos" column when dropping nan
                df2.dropna(subset=colNameSet, axis = 0, how='all', inplace=True)
                outfile2 = folder_plot_output_codon.joinpath(cond + '_simple1.pdf')

                def findCodon(row):
                    res = CodonList[int(row['pos']) - 1]
                    return(res)
                df2['Codon'] = df2.apply(lambda row: findCodon(row), axis=1)

                # Plot:
                tsv_plot_output_double(param.wtfile, param.mut_pos, cond, df2, df_se=df_se, outfilename=outfile2, version=2, scale = scale)

            # Plot3:
                outfile3 = folder_plot_output_codon.joinpath(cond + '_simple2.pdf')
                tsv_plot_output_double(param.wtfile, param.mut_pos, cond, df2, outfilename=outfile3, version=2, scale = scale)
        except:
            print("Step10: plots (codon version) creation failed!")
            exit()
        print("Step10: plots (codon version) created!")

def tsv_plot_output_aa_double_wrapper(param, folder_plot_input, folder_plot_output_codon, scale = 'max'):
    # Same as tsv_plot_output_aa_wrapper, optimized for double mutational scheme.
    if not len(glob.glob(folder_plot_output_codon.as_posix() + "/*")) == 0:
        print("Step10: plot output files already in plot_output/codon folder!")
    else:
        print("Step10: preparing plot output files (aa version) ...")
        folder_plot_output_codon.mkdir(parents=True, exist_ok=True)
        try:
            # Plot1:
            for cond in param.condition[:-1]:
                df1      = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.se.tsv'), delimiter='\t')
                outfile1 = folder_plot_output_codon.joinpath(cond + '_raw.pdf')

                # Add in missing columns to make it full 21 columns: data
                for aa in AAList:
                    if aa == 'Ter': aa = 'Stop' # make naming consistent
                    if not aa in df1.columns.values: df1[aa] = np.nan
                # Add in missing rows to make it full 64 rows: data
                row_insert_template = pd.Series([np.nan for i in range(0,df1.shape[1])], index=df1.columns.values)
                for group in CODON_GROUPS:
                    start = group[1]
                    if not start in df1['pos'].values:
                        row_insert_template['pos'] = start
                        df1 = df1.append(row_insert_template, ignore_index=True)
                df1.sort_values('pos', inplace=True, ascending=False)

                # Add in missing columns to make it full 21 columns: se
                for aa in AAList:
                    if not aa in df_se.columns.values: df_se[aa] = np.nan
                # Add in missing rows to make it full 64 rows: se
                row_insert_template = pd.Series([np.nan for i in range(0,df_se.shape[1])], index=df_se.columns.values)
                for group in CODON_GROUPS:
                    start = group[1]
                    if not start in df_se['pos'].values:
                        row_insert_template['pos'] = start
                        df_se = df_se.append(row_insert_template, ignore_index=True)
                df_se.sort_values('pos', inplace=True, ascending=False)

                df1.drop(columns = ['???'], inplace=True) # After adding undetermined nucleotide removal step the '???' column is no longer useful.
                df_se.drop(columns = ['???'], inplace=True)
                tsv_plot_output_aa_double(param.wtfile, param.mut_pos, cond, df1, df_se=df_se, outfilename=outfile1, scale = scale)

            # Plot2:
            for cond in param.condition[:-1]:
                df2      = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.tsv'), delimiter='\t')
                df_se    = pd.read_csv(folder_plot_input.joinpath(cond+'_aa.se.tsv'), delimiter='\t')

                # Add in missing columns to make it full 21 columns: data
                for aa in AAList:
                    if aa == 'Ter': aa = 'Stop' # make naming consistent
                    if not aa in df2.columns.values: df2[aa] = np.nan
                # Add in missing rows to make it full 64 rows: data
                row_insert_template = pd.Series([np.nan for i in range(0,df2.shape[1])], index=df2.columns.values)
                for group in CODON_GROUPS:
                    start = group[1]
                    if not start in df1['pos'].values:
                        row_insert_template['pos'] = start
                        df2 = df2.append(row_insert_template, ignore_index=True)
                df2.sort_values('pos', inplace=True, ascending=False)

                # Add in missing columns to make it full 21 columns: se
                for aa in AAList:
                    if not aa in df_se.columns.values: df_se[aa] = np.nan
                # Add in missing rows to make it full 64 rows: se
                row_insert_template = pd.Series([np.nan for i in range(0,df_se.shape[1])], index=df_se.columns.values)
                for group in CODON_GROUPS:
                    start = group[1]
                    if not start in df_se['pos'].values:
                        row_insert_template['pos'] = start
                        df_se = df_se.append(row_insert_template, ignore_index=True)
                df_se.sort_values('pos', inplace=True, ascending=False)

                # Drop columns with all NAN:
                df2.replace(to_replace='deplete', value=1000, inplace=True)
                df2 = df2.apply(lambda col: pd.to_numeric(col), axis=0)
                df2.dropna(axis = 1, how='all', inplace=True)
                # Drop rows with all NAN:
                cutoff = df2.shape[0]
                colNameSet = df2.columns.values[:-1] # ignore "pos" column when dropping nan
                df2.dropna(subset=colNameSet, axis = 0, how='all', inplace=True)
                outfile2 = folder_plot_output_codon.joinpath(cond + '_simple1.pdf')
                # Plot:
                tsv_plot_output_aa_double(param.wtfile, param.mut_pos, cond, df2, df_se=df_se, outfilename=outfile2, scale=scale)

            # Plot3:
                outfile3 = folder_plot_output_codon.joinpath(cond + '_simple2.pdf')
                tsv_plot_output_aa_double(param.wtfile, param.mut_pos, cond, df2, outfilename=outfile3, scale=scale)
        except:
            print("Step10: plots (aa version) creation failed!")
            exit()
        print("Step10: plots (aa version) created!")

"""
Below are wrappers for sinle and double mutational pipeline.
"""
def func_single_wrapper(param, workdir, scale = 'max'):
    # Wrapper for single mutation pipeline
    folder_data_sra = workdir.joinpath('data_sra')

    if len(param.ngs_data_local) == 0:
        get_ngs_data_ncbi(param, folder_data_sra)
    elif param.ngs_data_local[0] in ('to', 'NA', 'N/A'):
        get_ngs_data_ncbi(param, folder_data_sra)
    else:
        get_ngs_data_local(param, folder_data_sra)
    # If both NCBI sra and local storage provided, prefer to use local copy.

    folder_ref      = workdir.joinpath('ref')
    prep_ref(param, folder_ref)

    folder_merge    = workdir.joinpath('merge')
    fastq_merger(param, folder_data_sra, folder_merge)

    folder_first_mapper = workdir.joinpath('first_mapper')
    first_mapper(param, folder_merge, folder_first_mapper)

    folder_qc_INDEL     = workdir.joinpath('qc_INDEL')
    quality_control_INDEL(param, workdir, folder_qc_INDEL)

    folder_qc_library   = workdir.joinpath('qc_library')
    quality_control_library_wrapper(param, folder_qc_INDEL, folder_qc_library)

    folder_second_mapper = workdir.joinpath('second_mapper')
    second_mapper(param, folder_qc_library, folder_second_mapper)

    folder_enrich2_input = workdir.joinpath('enrich2_input')
    bam2enrich_wrapper(param, folder_second_mapper, folder_enrich2_input)

    folder_enrich2_json = workdir.joinpath('enrich2_json')
    folder_enrich2_output = workdir.joinpath('enrich2_output')
    enrich2_json_encoder_wrapper(param, folder_enrich2_json, folder_enrich2_input, folder_enrich2_output)
    enrich2_wrapper(folder_enrich2_json.joinpath('json.sh'), folder_enrich2_output)

    enrich2_hdf5_extractor_wrapper(param, folder_enrich2_output, _type='codon')
    enrich2_hdf5_extractor_wrapper(param, folder_enrich2_output, _type='aa')

    folder_plot_input = workdir.joinpath('plot_input')
    tsv_plot_input(param, folder_plot_input, folder_enrich2_output)

    folder_plot_output_codon = workdir.joinpath('plot_output/codon')
    tsv_plot_output_wrapper(param, folder_plot_input, folder_plot_output_codon, scale = scale)

    folder_plot_output_aa = workdir.joinpath('plot_output/aa')
    tsv_plot_output_aa_wrapper(param, folder_plot_input, folder_plot_output_aa, scale = scale)

    print("Congrats! All analyisis finished! Check plot_output folder!")

def func_double_wrapper(param, workdir, scale = 'max'):
    # Wrapper for double mutational pipeline
    folder_data_sra = workdir.joinpath('data_sra')

    if len(param.ngs_data_local) == 0:
        get_ngs_data_ncbi(param, folder_data_sra)
    elif param.ngs_data_local[0] in ('to', 'NA', 'N/A'):
        get_ngs_data_ncbi(param, folder_data_sra)
    else:
        get_ngs_data_local(param, folder_data_sra)
    # If both NCBI sra and local storage provided, prefer to use local copy.

    folder_ref      = workdir.joinpath('ref')
    prep_ref(param, folder_ref)

    folder_merge    = workdir.joinpath('merge')
    fastq_merger(param, folder_data_sra, folder_merge)

    folder_first_mapper = workdir.joinpath('first_mapper')
    first_mapper(param, folder_merge, folder_first_mapper)

    folder_qc_INDEL     = workdir.joinpath('qc_INDEL')
    quality_control_INDEL(param, workdir, folder_qc_INDEL)

    folder_qc_library   = workdir.joinpath('qc_library')
    quality_control_library_double_wrapper(param, folder_qc_INDEL, folder_qc_library)

    folder_second_mapper = workdir.joinpath('second_mapper')
    second_mapper(param, folder_qc_library, folder_second_mapper)

    folder_enrich2_input = workdir.joinpath('enrich2_input')
    bam2enrich_double_wrapper(param, folder_second_mapper, folder_enrich2_input)

    folder_enrich2_json = workdir.joinpath('enrich2_json')
    folder_enrich2_output = workdir.joinpath('enrich2_output')
    enrich2_json_encoder_double_wrapper(param, folder_enrich2_json, folder_enrich2_input, folder_enrich2_output)

    enrich2_tsv_extractor_double(param, folder_enrich2_output)

    folder_plot_input = workdir.joinpath('plot_input')
    tsv_plot_double_input(param, folder_plot_input, folder_enrich2_output)

    folder_plot_output_codon = workdir.joinpath('plot_output/codon')
    tsv_plot_output_double_wrapper(param, folder_plot_input, folder_plot_output_codon, scale = scale)

    folder_plot_output_aa = workdir.joinpath('plot_output/aa')
    tsv_plot_output_aa_double_wrapper(param, folder_plot_input, folder_plot_output_aa, scale = scale)

    print("Congrats! All analyisis finished! Check plot_output folder!")


if __name__ == '__main__':
    """
    In Version2.0, a double mutational pipeline was added.
    """

    def usg():
        return(
        """
        bash wrapper.sh
        -c/--config=path to configuration file
        -m/--mode=either single or double, default is single
        -s/--scale=Enrich2 score margin, default is max

        Example:
        bash wrapper.sh -c=/Path/XXX/config.txt -m=double -s=10

        Note, configuration file is required.
        """
        )

    parser = argparse.ArgumentParser(description="Pipeline for the analysis of DMS Data", usage=usg(), argument_default=argparse.SUPPRESS, add_help=False)
    if len(sys.argv) == 1:
        print(usg())
        exit()
    parser.add_argument('-c', '--config', required=True, help=argparse.SUPPRESS)
    parser.add_argument('-m', '--mode', type=str, help=argparse.SUPPRESS)
    parser.add_argument('-s', '--scale', type=str, help=argparse.SUPPRESS)

    args = parser.parse_args()

    workdir    = Path(Path.cwd()).parents[0]
    param = ConfigParam(args.config)

    Path(workdir.joinpath('TemFolder')).mkdir(parents=True, exist_ok=True) # create a temperate folder to contain tem files.

    if args.mode == 'single' and args.scale == 'max':
        print("Entering Single Mutation Mode ...")
        func_single_wrapper(param, workdir)
    elif args.mode == 'double' and args.scale == 'max':
        print("Entering Double Mutation Mode ...")
        func_double_wrapper(param, workdir)
    elif args.mode == 'single' and not args.scale == 'max':
        print("Entering Single Mutation Mode, score range set to " + str(args.scale) + " ...")
        func_single_wrapper(param, workdir, scale = args.scale)
    elif args.mode == 'double' and not args.scale == 'max':
        print("Entering Double Mutation Mode, score range set to " + str(args.scale) + " ...")
        func_double_wrapper(param, workdir, scale = args.scale)
    else:
        print(args.mode, args.scale)
