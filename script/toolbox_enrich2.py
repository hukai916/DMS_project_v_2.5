import matplotlib
matplotlib.use('Agg')
# Otherwise, python complains about not being built as a framework
# Must be put at the very beginning in case that other packages import matplotlib

from pathlib import Path
import json
import pandas as pd
import numpy as np
import pylab
import subprocess
import os
import shutil # use to remove non empty folders
from toolbox_kai import wt_codon
from constants import AA2CODON, CODON_GROUPS, CODON2AA, Property2AA, AA_GROUPS, CODON321, CodonList
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
from toolbox_matplotlib import recentered_cmap

## Note that pd.HDFStore requires: conda install -c conda-forge pytables

def enrich2_json_encoder(cond, cond_input, fastq_folder, mut, wt_code):
    """
    Customize enrich2 json configuration file by specifying:
    condition, mutation position, wt_codon and output folder name.
    """
    wt_code = str(wt_code)

    json_root = {"libraries": [],
                 "name"     : cond,
                 "output directory": ""
                 }

    json_library_lib    = {"fastq":
                                {
                                "filters":{},
                                "reads": Path(fastq_folder.joinpath('_'.join(
                                [cond_input, str(mut), str(mut+2), wt_code + '.fastq']))).as_posix(),
                                "reverse": False
                                },
                            "name": cond_input,
                            "report filtered reads": False,
                            "timepoint": 0,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": True,
                                    "reference offset": 0,
                                    "sequence": wt_code
                                    }
                                }
                            }

    json_library_cond    = {"fastq":
                                {
                                "filters":{},
                                "reads": Path(fastq_folder.joinpath('_'.join(
                                [cond, str(mut), str(mut+2), wt_code + '.fastq']))).as_posix(),
                                "reverse": False
                                },
                            "name": cond,
                            "report filtered reads": False,
                            "timepoint": 1,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": True,
                                    "reference offset": 0,
                                    "sequence": wt_code
                                    }
                                }
                            }

    json_root["libraries"] = [json_library_lib, json_library_cond]

    return(json.dumps(json_root, indent=4))

def enrich2_json_encoder_count(countfile_pathlist, output_folder):
    json_root = {"libraries": [],
                 "name"     : "Enrich2_count_wrapper",
                 "output directory": output_folder
                 }

    json_library_lib    = {"counts file": countfile_pathlist[0],
                           "fastq":
                                {
                                "filters":{},
                                "length": 'null',
                                "reads": 'null',
                                "reverse": False
                                },
                            "name": 'time0',
                            "report filtered reads": False,
                            "timepoint": 0,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": False,
                                    "reference offset": 0,
                                    "sequence": 'ATG'
                                    }
                                }
                            }

    json_library_cond    = {"counts file": countfile_pathlist[1],
                           "fastq":
                                {
                                "filters":{},
                                "length": 'null',
                                "reads": 'null',
                                "reverse": False
                                },
                            "name": 'time1',
                            "report filtered reads": False,
                            "timepoint": 1,
                            "variants":
                                {
                                "use aligner": False,
                                "wild type":
                                    {
                                    "coding": False,
                                    "reference offset": 0,
                                    "sequence": 'ATG'
                                    }
                                }
                            }

    json_root["libraries"] = [json_library_lib, json_library_cond]

    return(json.dumps(json_root, indent=4))

def enrich2_count_wrapper(c0_t, c1_t, c0, c1):
    """
    input 4 count integers and retrieve the Enrich2 output.
    c0_t: target count at time 0
    c1_t: target count at time 1
    c0: all count at time 0 including c0_t
    c1: all count at time 1 including c1_t
    """
    #print("inside count_wrapper: \n", c0_t, c1_t, c0, c1)

    workdir    = Path(Path.cwd()).parents[0]
    folder_json_tem   = workdir.joinpath("json_tem")
    if os.path.isdir(folder_json_tem): shutil.rmtree(folder_json_tem, ignore_errors=False, onerror=None) # Clean potential folders

    folder_json_tem.mkdir(parents=True, exist_ok=True)

    file_count0_tem   = folder_json_tem.joinpath('c0.tsv')
    file_count0_tem.touch(exist_ok=True)
    count0 = open(file_count0_tem, 'w+')
    count0.write("\tcount\n")
    count0.write('c_t\t' +  str(c0_t) + "\n")
    count0.write('c_other\t' + str(c0-c0_t))
    count0.close()

    file_count1_tem   = folder_json_tem.joinpath('c1.tsv')
    file_count1_tem.touch(exist_ok=True)
    count1 = open(file_count1_tem, 'w+')
    count1.write("\tcount\n")
    count1.write('c_t\t' + str(c1_t) + '\n')
    count1.write('c_other\t' +  str(c1-c1_t))
    count1.close()

    res_json = enrich2_json_encoder_count([file_count0_tem.as_posix(),
                                           file_count1_tem.as_posix()],
                                           folder_json_tem.as_posix())

    file_json_tem     = folder_json_tem.joinpath('tem.json')

    file_json_tem.touch(exist_ok=True)
    jsonfile = open(file_json_tem, 'w+')
    for x in res_json: jsonfile.write(x)
    jsonfile.close()

    json_command = ' '.join(['enrich_cmd', file_json_tem.as_posix().replace(' ', '\ '),
                        "ratios complete --no-plots --output-dir ",
                        folder_json_tem.as_posix().replace(' ', '\ ')])

    json_commandfile = open(folder_json_tem.joinpath('json.sh'), 'w+')
    json_commandfile.write("CONDA_BASE=$(conda info --base)\n")
    json_commandfile.write("source $CONDA_BASE/etc/profile.d/conda.sh\n")
    json_commandfile.write("conda activate py2-dms\n")
    json_commandfile.write(json_command)
    json_commandfile.close()


    command = "bash " + folder_json_tem.joinpath('json.sh').as_posix().replace(' ', '\ ')
    subprocess.run(command, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    my_store = pd.HDFStore(folder_json_tem.joinpath('Enrich2_count_wrapper_sel.h5'))
    df = my_store.select("/main/variants/scores")
    my_store.close()

    shutil.rmtree(folder_json_tem) # clean the tem folder
    return(df.iloc[0,:])

def extract_codon(sentence, wt_codon):
    wt_codon = list(wt_codon)
    tem = sentence.split()
    tem = ''.join(tem[::2])
    for i in range(len(tem)):
        if tem[i] == str(1):
            wt_codon[0] = tem[i+3]
        elif tem[i] == str(2):
            wt_codon[1] = tem[i+3]
        elif tem[i] == str(3):
            wt_codon[2] = tem[i+3]
    return("".join(wt_codon))

def enrich2_hdf5_extractor(raw_count1, raw_count2, wt_codon, score_file):
    """
    This function generates a big DataFrame that combines data from two raw_counts file and
    one score_file. We need to incorporate "raw count files" because "main" files of Enrich2 output
    only keep items that show in both raw files.
    """

    my_store = pd.HDFStore(raw_count1)
    df_raw_count1 = my_store.select("/raw/variants/counts")
    df_raw_count1.reset_index(inplace=True)
    my_store.close()

    my_store = pd.HDFStore(raw_count2)
    df_raw_count2 = my_store.select("/raw/variants/counts")
    df_raw_count2.reset_index(inplace=True)
    my_store.close()

    my_store = pd.HDFStore(score_file)
    df_score_file = my_store.select("/main/variants/scores")
    df_score_file.reset_index(inplace=True)
    my_store.close()

    df_raw_count1['codon'] = df_raw_count1.apply(lambda row: extract_codon(row['index'], wt_codon), axis=1)
    df_raw_count1.drop(['index'], axis=1, inplace=True)
    df_raw_count1.rename(columns={'count': 'count0'}, inplace=True)

    df_raw_count2['codon'] = df_raw_count2.apply(lambda row: extract_codon(row['index'], wt_codon), axis=1)
    df_raw_count2.drop(['index'], axis=1, inplace=True)
    df_raw_count2.rename(columns={'count': 'count1'}, inplace=True)

    df_score_file['codon'] = df_score_file.apply(lambda row: extract_codon(row['index'], wt_codon), axis=1)
    df_score_file.drop(['index'], axis=1, inplace=True)

    df_raw = pd.merge(df_raw_count1, df_raw_count2, how='outer', left_on='codon', right_on='codon')
    df_all = pd.merge(df_raw, df_score_file, how='outer', left_on='codon', right_on='codon')

    df_all['mut_aa'] = df_all.apply(lambda row: CODON2AA[row['codon']][0] if row['codon'] in CODON2AA else np.nan, axis=1)
    df_all['wt_aa']  = df_all.apply(lambda row: CODON2AA[wt_codon][0], axis=1)

    return(df_all)

def extract_aa(sentence):
    """
    wt_codon = list(wt_codon)
    tem = sentence.split()
    tem = ''.join(tem[::2])
    for i in range(len(tem)):
        if tem[i] == str(1):
            wt_codon[0] = tem[i+3]
        elif tem[i] == str(2):
            wt_codon[1] = tem[i+3]
        elif tem[i] == str(3):
            wt_codon[2] = tem[i+3]
    """
    wt_aa = sentence
    mut_aa = sentence
    if sentence == '_sy':
        wt_aa = 'WT'
        mut_aa = 'SY'
    elif sentence == '_wt':
        wt_aa = 'WT'
        mut_aa = 'WT'
    else:
        wt_aa  = sentence[2:5]
        mut_aa = sentence[-3:]



    return([wt_aa, mut_aa])

def enrich2_hdf5_extractor2(wt_codon, score_file):
    """
    Perform a similar operation and generate amino acid based data extraction.
    We can actually obtain all required info from _sel.h5 file instead of separately
    retrieving them from _lib.h5 files. _sel.h5 might be actually based on the two _lib.h5
    files.
    Score_file is the sole _sel.h5 file.
    """
    my_store = pd.HDFStore(score_file)
    df_raw_count = my_store.select('/main/synonymous/counts_unfiltered')
    my_store.close()

    my_store = pd.HDFStore(score_file)
    df_raw_score =  my_store.select('/main/synonymous/scores')
    my_store.close()

    df_raw = pd.merge(df_raw_count, df_raw_score, how='outer', left_index=True, right_index=True)
    df_raw.rename(columns={'c_0':'count0', 'c_1':'count1'}, inplace=True)

    df_raw.reset_index(inplace=True)

    df_raw['aa_list'] = df_raw.apply(lambda row: extract_aa(row['index']), axis=1)
    df_raw['wt_aa']   = df_raw.apply(lambda row: row['aa_list'][0], axis=1)
    df_raw['mut_aa']   = df_raw.apply(lambda row: row['aa_list'][1], axis=1)
    df_raw.drop(['aa_list', 'index'], axis=1, inplace=True)
    return(df_raw)

def tsv_plot_output(wtfile, cond, df, df_se=pd.DataFrame(), outfilename='test.pdf', version=1, scale = 'max'):
    """
    Read in score and se dataframes, output plots accordingly.
    wtfile: WT sequence file, used for y-labeling
    cond: The experiment condition, used for title info
    df: enrich score dataframe, from plot_input folder, can also use customized df
    df_se: enrich score dataframe, from plot_input folder, will match df (if customized) automatically
    outfilename: outputfile path + filename, where to store the output.
    Version controls the display of amino acid labels: in figure1(raw figure), use 3-letter; in figure2 and 2 (simple1 and simple2), use 1-letter.
    """

    row_number, col_number = df.shape

    # Create Figure object that consists of two axes:
    grid = GridSpec(2,1, height_ratios=[row_number,1], hspace=1/(row_number+1)*2)
    fig  = plt.figure()
    fig.set_size_inches((col_number+8)*0.2, (10+row_number)*0.2)
    # Later, will use subplots_adjust to make sure each square is 1*0.2^2 inches^2


    ### Up to now, there is no axes object created, therefore, no plot will be shown.

    # Create two axes object subplot.
    mesh_ax = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    # Adjust subplots to make room for add-on information:
        ## top: (title + amino acid grouping): 5 unit (each unit is 1*0.2 inch)
        ## padding between main figure and bar: 1 unit
        ## bottom: 3 units
        ## left 2.5 units
        ## right: 3.5 units
    plt.gcf().subplots_adjust(top=(row_number+5)/(row_number+10),
                              bottom=3/(row_number+10),
                              left=3/(col_number+8),
                              right=(col_number+3)/(col_number+8))

    # Replace 'deplete' with an arbituraliy large value, say, 1000
    df_raw = df.set_index(['pos'])
    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: float(row[col])
                                   if not row[col] == 'deplete'
                                   else 1000, axis=1)

    # Find score range:
    ls = df_raw.values.reshape(1,-1).tolist()[0]
    ls = [x for x in ls if not np.isnan(x)]

    _ls = sorted(set(ls))
    vmin = _ls[0]
    if not _ls[-1] == 1000:
        vmax = _ls[-1]
    else:
        vmax = _ls[-2]

### Try to fix the color overflow issue:
    if scale == 'max':
        colorlim = max(-vmin, vmax)
    else:
        colorlim = float(scale)

    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: colorlim if row[col] > colorlim and not row[col] == 1000 else row[col], axis=1)

    # Prepare a numpy array and mask NaNs for later plotting:
    arr_masked = np.ma.array(df_raw, mask=np.isnan(df_raw))

    # Get color map and set colors:
    cmap = plt.get_cmap("RdBu_r")

    # Recenter color map:
    cmap = recentered_cmap(cmap, -colorlim, colorlim)

    # Rescale the color by cutting a fragment from it:
    colors = [cmap(i) for i in range(65,200)]  # R -> G -> B
    cmap_new = LinearSegmentedColormap.from_list('place_holder', colors)
    # The new cmap_new will take advantage of old cmap, "shrink" the spectrum to make it lighter for better printing.

    # Set grey to NaN values, and black to totally depleted ones:
    cmap_new.set_bad("#808080",1)
    rgba = cmap_new(0) # The darkest in the new cmap_new, equal to cmap(65)
    cmap_new.set_over(rgba)

    # Plot the heatmap: return the value as mesn_pcolor as a mappable object in order to add color_bar:
    mesh_pcolor = mesh_ax.pcolormesh(arr_masked, cmap=cmap_new, vmin=-colorlim, vmax=colorlim)

    ## Below are modifications to plotting:

    # Add in the color bar
    cbar = fig.colorbar(mesh_pcolor, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Enrich2 Score")

    # Set mesh_ax y labels:
    mesh_wt_ax = mesh_ax.twinx()
    mesh_ax.set_ylabel("Position in WT")
    mesh_wt_ax.set_ylabel("WT sequence")

    # Set mesh_ax title:
    mesh_ax.set_title("Codon Enrichment for Experiment: " + cond, pad=50)
    # pad uses points, not sure the relationship between point and inch.
    # 50 works nicely it seems.

    # Add in column information:
    for i, x in enumerate(list(df_raw.columns)):
        mesh_ax.text(i + 0.5, len(df_raw.index) + 1, x,
        horizontalalignment="center",
        verticalalignment="center",
        rotation = 90)

    # Add in amino acid grouping information:
    new_CODON_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(CODON_GROUPS)):
        aaName  = CODON_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in AA2CODON[aaName]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaName, cstart, cend)
            new_CODON_GROUPS.append(newTuple)

    for codon, start, end in new_CODON_GROUPS:
        if version == 2: codon = CODON321[codon]
        mesh_ax.text((end - start + 1) / 2 + start,
            row_number + 2.5, codon,
            horizontalalignment="center",
            verticalalignment="center")

        bar = Line2D([start + 0.125, end + 1 - 0.125],
                [row_number + 2, row_number + 2], color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)

        # Add in the deliminator for amino acid groups:
        delimBar = Line2D([end + 1, end + 1],
                [0, len(df_raw.index) + 1], color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    WT_codon  = [str(wt_codon(wtfile, int(x*3-2))) for x in df['pos'].values]
    wt_aa     = [str(CODON2AA[codon][0]) for codon in WT_codon]

    # Add in mesh_ax y-label: aa coordinates in WT sequence:
    ypos = np.arange(len(df['pos'])) + 0.5
    mesh_ax.set_yticks(ypos)
    mesh_ax.set_yticklabels(list(map(int, df['pos'])), ha='right')

    # Add in wt_ax y-label: WT aa and codon:
    mesh_wt_ax.set_ylim(0, len(df['pos']))
    mesh_wt_ax.set_yticks(ypos)
    label = ['-'.join(element) for element in zip(WT_codon, wt_aa)]
    mesh_wt_ax.set_yticklabels(label, ha='left')

    # Add in WT label onto corresponding cell:
    x_coordinate = [list(df.columns).index(x) if x in df.columns else np.nan for x in WT_codon]
    y_coordinate = list(range(df.shape[1]))
    wt_xy        = zip(x_coordinate, y_coordinate)
    for x, y in wt_xy:
        if not x == np.nan:
            mesh_ax.add_patch(Circle((x + 0.5, y + 0.5), .1666,
            fill=True, facecolor="black",
            edgecolor="none", alpha=0.5))

    # Make the figure cleaner by removing ticks:
    mesh_ax.tick_params(bottom=False, left=False)
    mesh_wt_ax.tick_params(right=False)
    mesh_ax.get_xaxis().set_visible(False)
    cbar_ax.tick_params(bottom=False)

    # Add in SE if specified:
    if not df_se.shape == pd.DataFrame().shape: # Add in SE
        # Subset df_se to match df that might be truncated:
        for col in df_se.columns:
            if not col in df.columns: df_se.drop(columns=col, inplace=True)
        df_se = df_se[df_se['pos'].isin(df['pos'])]

        # Below is to find max se value:
        se_df = df_se.drop(columns=['pos'])
        tem = sorted(se_df.values.reshape(1,-1).tolist()[0])
        tem = [x for x in tem if not np.isnan(x)]
        se_max = tem[-1]

        for row in range(len(df_se['pos'])):
            for col in range(len(df_se.columns)-1):
                se_value = df_se.iloc[row, col]
                corner_dist = (se_max - se_value)/(2 * se_max)
                corner_dist = (1 - se_value) / 2
                diag = Line2D([col + corner_dist, col + 1 - corner_dist],
                        [row + corner_dist, row + 1 - corner_dist], color="black")

                if se_value > 0.02 and df_raw.iloc[row, col] != 1000: # se_value below 0.02 will not be displayed so as totally depleted ones
                    mesh_ax.add_line(diag)

    pylab.savefig(outfilename)

def tsv_plot_output_aa(wtfile, cond, df, df_se=pd.DataFrame(), outfilename='test.pdf', scale = 'max'):
    """
    Read in score and se dataframes, output plots accordingly.
    wtfile: WT sequence file, used for y-labeling
    cond: The experiment condition, used for title info
    df: enrich score dataframe, from plot_input folder, can also used customized df
    df_se: enrich score dataframe, from plot_input folder, will match df (if customized) automatically
    outfilename: outputfile path + filename, where to store the output
    This function is customized for aa.
    """

    row_number, col_number = df.shape

    # Create Figure object that consists of two axes:
    grid = GridSpec(2,1, height_ratios=[row_number,1], hspace=1/(row_number+1)*2)
    fig  = plt.figure()
    fig.set_size_inches((col_number+8)*0.2, (10+row_number)*0.2)
    # Later, will use subplots_adjust to make sure each square is 1*0.2^2 inches^2


    ### Up to now, there is no axes object created, therefore, no plot will be shown.

    # Create two axes object subplot.
    mesh_ax = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    # Adjust subplots to make room for add-on information:
        ## top: (title + amino acid grouping): 5 unit (each unit is 1*0.2 inch)
        ## padding between main figure and bar: 1 unit
        ## bottom: 3 units
        ## left 2.5 units
        ## right: 3.5 units
    plt.gcf().subplots_adjust(top=(row_number+5)/(row_number+10),
                              bottom=3/(row_number+10),
                              left=5/(col_number+8),
                              right=(col_number+5)/(col_number+8))

    # Replace 'deplete' with an arbituraliy large value, say, 1000
    df_raw = df.set_index(['pos'])
    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: float(row[col])
                                   if not row[col] == 'deplete'
                                   else 1000, axis=1)
    # Find score range:
    ls = df_raw.values.reshape(1,-1).tolist()[0]
    ls = [x for x in ls if not np.isnan(x)]
    _ls = sorted(set(ls))
    vmin = _ls[0]
    if not _ls[-1] == 1000:
        vmax = _ls[-1]
    else:
        vmax = _ls[-2]

### Try to fix the color overflow issue:
    if scale == 'max':
        colorlim = max(-vmin, vmax)
    else:
        colorlim = float(scale)

    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: colorlim if row[col] > colorlim and not row[col] == 1000 else row[col], axis=1)

    # Prepare a numpy array and mask NaNs for later plotting:
    arr_masked = np.ma.array(df_raw, mask=np.isnan(df_raw))

    # Get color map and set colors:
    cmap = plt.get_cmap("RdBu_r")

    # Recenter color map:
    cmap = recentered_cmap(cmap, -colorlim, colorlim)

    colors = [cmap(i) for i in range(65,200)]  # R -> G -> B
    cmap_new = LinearSegmentedColormap.from_list('place_holder', colors)
    # The new cmap_new will take advantage of old cmap, "shrink" the spectrum to make it lighter for better printing.

    # Set grey to NaN values, and black to totally depleted ones:
    cmap_new.set_bad("#808080",1)
    rgba = cmap_new(0) # The darkest in the new cmap_new, equal to cmap(65)
    cmap_new.set_over(rgba)

    # Plot the heatmap: return the value as mesn_pcolor as a mappable object in order to add color_bar:
    mesh_pcolor = mesh_ax.pcolormesh(arr_masked, cmap=cmap_new, vmin=-colorlim, vmax=colorlim)

    ## Below are modifications to plotting:

    # Add in the color bar
    cbar = fig.colorbar(mesh_pcolor, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Enrich2 Score")

    # Set mesh_ax y labels:
    mesh_ax.set_ylabel("Position in WT")

    # Set mesh_ax title:
    mesh_ax.set_title("Amino Acid Enrichment for Experiment: " + cond, pad=50)
    # pad uses points, not sure the relationship between point and inch.
    # 50 works nicely it seems.

    # Add in column information:
    for i, x in enumerate(list(df_raw.columns)):
        mesh_ax.text(i + 0.5, len(df_raw.index) + 1, x,
        horizontalalignment="center",
        verticalalignment="center",
        rotation = 90)

    # Add in amino acid grouping information:
    new_AA_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(AA_GROUPS)):
        aaProperty  = AA_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in Property2AA[aaProperty]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaProperty, cstart, cend)
            new_AA_GROUPS.append(newTuple)

    for codon, start, end in new_AA_GROUPS:
        mesh_ax.text((end - start + 1) / 2 + start,
            row_number + 2.5, codon,
            horizontalalignment="center",
            verticalalignment="center")

        bar = Line2D([start + 0.125, end + 1 - 0.125],
                [row_number + 2, row_number + 2], color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)


    WT_codon  = [str(wt_codon(wtfile, int(x*3-2))) for x in df['pos'].values]
    wt_aa     = [str(CODON2AA[codon][0]) for codon in WT_codon]

    # Add in mesh_ax y-label: aa coordinates in WT sequence:
    ypos = np.arange(len(df['pos'])) + 0.5
    mesh_ax.set_yticks(ypos)
    labelPos = list(map(str,list(map(int, df['pos']))))
    labelAA  = wt_aa
    labelCombine = ['-'.join(item) for item in zip(labelPos, labelAA)]
    mesh_ax.set_yticklabels(labelCombine, ha='right')

    # Add in deliminator horizontally:
    for i in range(0, df_raw.shape[0]):
        delimBar = Line2D([0, df_raw.shape[1]],[i, i],
                transform=mesh_ax.transData, color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    # Add in WT label onto corresponding cell:
    WT_aa = [CODON2AA[x][0] for x in WT_codon]
    x_coordinate = [list(df.columns).index(x) if x in df.columns else np.nan for x in WT_aa]
    y_coordinate = list(range(df.shape[1]))

    wt_xy        = zip(x_coordinate, y_coordinate)
    for x, y in wt_xy:
        if not x == np.nan:
            mesh_ax.add_patch(Circle((x + 0.5, y + 0.5), .1666,
            fill=True, facecolor="black",
            edgecolor="none", alpha=0.5))

    # Make the figure cleaner by removing ticks:
    mesh_ax.tick_params(bottom=False, left=False)
    mesh_ax.get_xaxis().set_visible(False)
    cbar_ax.tick_params(bottom=False)

    # Add in SE if specified:
    if not df_se.shape == pd.DataFrame().shape: # Add in SE
        # Subset df_se to match df that might be truncated:
        for col in df_se.columns:
            if not col in df.columns: df_se.drop(columns=col, inplace=True)
        df_se = df_se[df_se['pos'].isin(df['pos'])]

        # Below is to find max se value:
        se_df = df_se.drop(columns=['pos'])
        tem = sorted(se_df.values.reshape(1,-1).tolist()[0])
        tem = [x for x in tem if not np.isnan(x)]
        se_max = tem[-1]

        for row in range(len(df_se['pos'])):
            for col in range(len(df_se.columns)-1):
                se_value = df_se.iloc[row, col]
                corner_dist = (se_max - se_value)/(2 * se_max)
                corner_dist = (1 - se_value) / 2
                diag = Line2D([col + corner_dist, col + 1 - corner_dist],
                        [row + corner_dist, row + 1 - corner_dist], color="black")

                if se_value > 0.02 and df_raw.iloc[row, col] != 1000: # se_value below 0.02 will not be displayed so as totally depleted ones
                    mesh_ax.add_line(diag)

    pylab.savefig(outfilename)

def tsv_plot_output_double(wtfile, wt_mut, cond, df, df_se=pd.DataFrame(), outfilename='test.pdf', version=1, scale = 'max'):
    """
    Read in score and se dataframes, output plots accordingly. Optimized for double mutation scheme.
    wtfile: WT sequence file, used for y-labeling
    wt_mut: mutation positions
    cond: The experiment condition, used for title info
    df: enrich score dataframe, from plot_input folder, can also use customized df
    df_se: enrich score dataframe, from plot_input folder, will match df (if customized) automatically
    outfilename: outputfile path + filename, where to store the output.
    Version controls the display of amino acid labels: in figure1(raw figure), use 3-letter; in figure2 and 2 (simple1 and simple2), use 1-letter.
    Optimized for double mutation assays.
    """

    row_number, col_number = df.shape
    col_number = col_number - 1

    # Create Figure object that consists of two axes:
    grid = GridSpec(2,1, height_ratios=[row_number,1], hspace=1/(row_number+1)*2)
    fig  = plt.figure()
    fig.set_size_inches((col_number+8)*0.2, (10+row_number)*0.2)
    # Later, will use subplots_adjust to make sure each square is 1*0.2^2 inches^2

    ### Up to now, there is no axes object created, therefore, no plot will be shown.

    # Create two axes object subplot.
    mesh_ax = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    # Adjust subplots to make room for add-on information:
        ## top: (title + amino acid grouping): 5 unit (each unit is 1*0.2 inch)
        ## padding between main figure and bar: 1 unit
        ## bottom: 3 units
    plt.gcf().subplots_adjust(top=(row_number+4.2)/(row_number+10),
                              bottom=2.2/(row_number+10),
                              left=2/(col_number+8),
                              right=(col_number+2)/(col_number+8))

    # Replace 'deplete' with an arbituraliy large value, say, 1000
    df_raw0 = df.set_index(['pos'])
    df_raw  = df_raw0.drop(columns=['Codon'])

    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: float(row[col])
                                   if not row[col] == 'deplete'
                                   else 1000, axis=1)
    # Find score range:
    ls = df_raw.values.reshape(1,-1).tolist()[0]
    ls = [x for x in ls if not np.isnan(x)]

    _ls = sorted(set(ls))
    vmin = _ls[0]
    if not _ls[-1] == 1000:
        vmax = _ls[-1]
    else:
        vmax = _ls[-2]

### Try to fix the color overflow issue:
    if scale == 'max':
        colorlim = max(-vmin, vmax)
    else:
        colorlim = float(scale)

    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: colorlim if row[col] > colorlim and not row[col] == 1000 else row[col], axis=1)

    # Prepare a numpy array and mask NaNs for later plotting:
    arr_masked = np.ma.array(df_raw, mask=np.isnan(df_raw))

    # Get color map and set colors:
    cmap = plt.get_cmap("RdBu_r")

    # Recenter color map:
    cmap = recentered_cmap(cmap, -colorlim, colorlim)

    # Rescale the color by cutting a fragment from it:
    colors = [cmap(i) for i in range(65,200)]  # R -> G -> B
    cmap_new = LinearSegmentedColormap.from_list('place_holder', colors)
    # The new cmap_new will take advantage of old cmap, "shrink" the spectrum to make it lighter for better printing.

    # Set grey to NaN values, and black to totally depleted ones:
    cmap_new.set_bad("#808080",1)
    rgba = cmap_new(0) # The darkest in the new cmap_new, equal to cmap(65)
    cmap_new.set_over(rgba)

    # Plot the heatmap: return the value as mesn_pcolor as a mappable object in order to add color_bar:
    mesh_pcolor = mesh_ax.pcolormesh(arr_masked, cmap=cmap_new, vmin=-colorlim, vmax=colorlim)

    ## Below are modifications to plotting:

    # Add in the color bar
    cbar = fig.colorbar(mesh_pcolor, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Enrich2 Score")

    # Set mesh_ax title:
    WT_codon  = [str(wt_codon(wtfile, wt_mut[0][0])), str(wt_codon(wtfile, wt_mut[0][1]))]
    WT_site1 = np.nan
    WT_site2 = np.nan
    # These lines for determining WT label coordinates: site1 determines row position while site2 determines column position

    Col_site = " (Column: " + str(int((wt_mut[0][1] + 2)/3)) + "-" + CODON2AA[WT_codon[1]][0] + " (" + WT_codon[1] + "))"
    mesh_ax.set_title("Codon Enrichment (Double) for: " + cond + "\n" + Col_site, pad=50)
    # pad uses points, not sure the relationship between point and inch.
    # 50 works nicely it seems.

    # Add in codon information: columns
    for i, x in enumerate(list(df_raw.columns)):
        mesh_ax.text(i + 0.5, len(df_raw.index) + 1, x,
        horizontalalignment="center",
        verticalalignment="center",
        rotation = 90)
        if x == WT_codon[1]: WT_site2 = i

    # Add in amino acid grouping information: columns
    new_CODON_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(CODON_GROUPS)):
        aaName  = CODON_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in AA2CODON[aaName]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaName, cstart, cend)
            new_CODON_GROUPS.append(newTuple)

    for codon, start, end in new_CODON_GROUPS:
        if version == 2: codon = CODON321[codon]
        mesh_ax.text((end - start + 1) / 2 + start,
            row_number + 2.5, codon,
            horizontalalignment="center",
            verticalalignment="center")

        bar = Line2D([start + 0.125, end + 1 - 0.125],
                [row_number + 2, row_number + 2], color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)

    # Add in the deliminator for amino acid groups: columns
        delimBar = Line2D([end + 1, end + 1],
                [0, len(df_raw.index) + 1], color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    # Add in codon information: row
    for i, x in enumerate(list(df['pos'][::-1])):
        mesh_ax.text(col_number, i+1-0.5, list(df['Codon'])[i], # Note that if -0.5, visually not perfect but this is the best achievable
                horizontalalignment="center",
                verticalalignment="center",)
        if CodonList[::-1][int(x) - 1] == WT_codon[0]: WT_site1 = i

    # Add in amino acid grouping information: rows
    new_CODON_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(CODON_GROUPS)):
        aaName  = CODON_GROUPS[i][0]
        count = 0
        for item in df['pos'].values:
            if CodonList[int(item)-1] in AA2CODON[aaName]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaName, cstart, cend)
            new_CODON_GROUPS.append(newTuple)

    for codon, start, end in new_CODON_GROUPS:
        if version == 2: codon = CODON321[codon]
        y = row_number - (end - start)/2 - start - 0.5
        mesh_ax.text(col_number+1.85, y, codon,
            horizontalalignment="center",
            verticalalignment="center")

        x = [col_number + 1, col_number + 1]
        y = [row_number - start - 1 + 0.875, row_number - end - 0.875]
        bar = Line2D(x, y, color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)

    # Add in the deliminator for amino acid groups: rows
        x = [0, col_number]
        y = [row_number - end - 1, row_number - end - 1]
        delimBar = Line2D(x, y, color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

        """
    # Legacy, change the position to figure title and right side.
    # Add in mutation site info onto the left y-axis:
        x = -2
        y = row_number / 3 * 2

        mesh_ax.text(x, y, "Row:\n" + str(int((wt_mut[0][0] + 2)/3)) + "-" + CODON2AA[WT_codon[0]][0] + "\n(" + WT_codon[0] + ")",
        horizontalalignment="center",
        verticalalignment="center")

        mesh_ax.text(x, y/2, "Column:\n" + str(int((wt_mut[0][1] + 2)/3)) + "-" + CODON2AA[WT_codon[1]][0] + "\n(" + WT_codon[1] + ")",
        horizontalalignment="center",
        verticalalignment="center")
        """

    # Add in first mutation site (per row) info onto figure right margin
    x = col_number + 3.5
    y = row_number / 2
    mesh_ax.text(x, y, "(Row: " + str(int((wt_mut[0][0] + 2)/3)) + "-" + CODON2AA[WT_codon[0]][0] + " (" + WT_codon[0] + "))",
    horizontalalignment="center",
    verticalalignment="center",
    rotation = 270)

    # Add in WT label onto corresponding cells:
    wt_aa     = [str(CODON2AA[codon][0]) for codon in WT_codon]

    if not WT_site1 == np.nan and not WT_site2 == np.nan:
            mesh_ax.add_patch(Circle((WT_site2 + 0.5, WT_site1 + 0.5), .1666,
            fill=True, facecolor="black",
            edgecolor="none", alpha=0.5))

    # Make the figure cleaner by removing ticks:
    mesh_ax.tick_params(bottom=False, left=False)
    mesh_ax.get_xaxis().set_visible(False)
    mesh_ax.get_yaxis().set_visible(False)
    cbar_ax.tick_params(bottom=False)

    # Add in SE if specified:
    if not df_se.shape == pd.DataFrame().shape: # Add in SE
        # Subset df_se to match df that might be truncated:
        for col in df_se.columns:
            if not col in df.columns: df_se.drop(columns=col, inplace=True)
        df_se = df_se[df_se['pos'].isin(df['pos'])]

        # Below is to find max se value:
        se_df = df_se.drop(columns=['pos'])
        tem = sorted(se_df.values.reshape(1,-1).tolist()[0])
        tem = [x for x in tem if not np.isnan(x)]
        se_max = tem[-1]

        for row in range(len(df_se['pos'])):
            for col in range(len(df_se.columns)-1):
                se_value = df_se.iloc[row, col]
                corner_dist = (se_max - se_value)/(2 * se_max)
                corner_dist = (1 - se_value) / 2
                diag = Line2D([col + corner_dist, col + 1 - corner_dist],
                        [row + corner_dist, row + 1 - corner_dist], color="black")

                if se_value > 0.02 and df_raw.iloc[row, col] != 1000: # se_value below 0.02 will not be displayed so as totally depleted ones
                    mesh_ax.add_line(diag)

    pylab.savefig(outfilename)

def tsv_plot_output_aa_double(wtfile, wt_mut, cond, df, df_se=pd.DataFrame(), outfilename='test.pdf', version=1, scale = 'max'):
    # Same as tsv_plot_output_double but for aa mode.
    row_number, col_number = df.shape

    # Create Figure object that consists of two axes:
    grid = GridSpec(2,1, height_ratios=[row_number,1], hspace=1/(row_number+1)*2)
    fig  = plt.figure()
    fig.set_size_inches((col_number+10)*0.2, (10+row_number)*0.2)
    # Later, will use subplots_adjust to make sure each square is 1*0.2^2 inches^2

    ### Up to now, there is no axes object created, therefore, no plot will be shown.

    # Create two axes object subplot.
    mesh_ax = plt.subplot(grid[0])
    cbar_ax = plt.subplot(grid[1])
    # Adjust subplots to make room for add-on information:
        ## top: (title + amino acid grouping): 5 unit (each unit is 1*0.2 inch)
        ## padding between main figure and bar: 1 unit
        ## bottom: 3 units
    plt.gcf().subplots_adjust(top=(row_number+4.2)/(row_number+10),
                              bottom=2.2/(row_number+10),
                              left=2/(col_number+10),
                              right=(col_number+2)/(col_number+10))

    # Replace 'deplete' with an arbituraliy large value, say, 1000
    df_raw = df.set_index(['pos'])
    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: float(row[col])
                                   if not row[col] == 'deplete'
                                   else 1000, axis=1)
    # Find score range:
    ls = df_raw.values.reshape(1,-1).tolist()[0]
    ls = [x for x in ls if not np.isnan(x)]

    _ls = sorted(set(ls))
    vmin = _ls[0]
    if not _ls[-1] == 1000:
        vmax = _ls[-1]
    else:
        vmax = _ls[-2]

### Try to fix the color overflow issue:
    if scale == 'max':
        colorlim = max(-vmin, vmax)
    else:
        colorlim = float(scale)

    for col in df_raw.columns:
        df_raw[col] = df_raw.apply(lambda row: colorlim if row[col] > colorlim and not row[col] == 1000 else row[col], axis=1)

    # Prepare a numpy array and mask NaNs for later plotting:
    arr_masked = np.ma.array(df_raw, mask=np.isnan(df_raw))

    # Get color map and set colors:
    cmap = plt.get_cmap("RdBu_r")

    # Recenter color map:
    cmap = recentered_cmap(cmap, -colorlim, colorlim)

    # Rescale the color by cutting a fragment from it:
    colors = [cmap(i) for i in range(65,200)]  # R -> G -> B
    cmap_new = LinearSegmentedColormap.from_list('place_holder', colors)
    # The new cmap_new will take advantage of old cmap, "shrink" the spectrum to make it lighter for better printing.

    # Set grey to NaN values, and black to totally depleted ones:
    cmap_new.set_bad("#808080",1)
    rgba = cmap_new(0) # The darkest in the new cmap_new, equal to cmap(65)
    cmap_new.set_over(rgba)

    # Plot the heatmap: return the value as mesn_pcolor as a mappable object in order to add color_bar:
    mesh_pcolor = mesh_ax.pcolormesh(arr_masked, cmap=cmap_new, vmin=-colorlim, vmax=colorlim)

    ## Below are modifications to plotting:

    # Add in the color bar
    cbar = fig.colorbar(mesh_pcolor, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("Enrich2 Score")

    # Set mesh_ax title:
    WT_codon  = [str(wt_codon(wtfile, wt_mut[0][0])), str(wt_codon(wtfile, wt_mut[0][1]))]
    WT_site1 = np.nan
    WT_site2 = np.nan
    # These lines for determining WT label coordinates: site1 determines row position while site2 determines column position

    Col_site = " (Column: " + str(int((wt_mut[0][1] + 2)/3)) + "-" + CODON2AA[WT_codon[1]][0] + " (" + WT_codon[1] + "))"
    mesh_ax.set_title("Amino Acid Enrichment (Double) for : " + cond + "\n" + Col_site, pad=50)
    # pad uses points, not sure the relationship between point and inch.
    # 50 works nicely it seems.

    # Reoder the rows according to Property order list:
    PropertyOrderList = Property2AA['Polar'] + Property2AA['Charged'] + Property2AA['Non-polar'] + ('Stop',)
    df_row_order = pd.DataFrame()
    for aa in PropertyOrderList:
        pos = [group[1] for group in CODON_GROUPS if group[0] == aa][0]
        row = df[df['pos'] == pos]
        df_row_order = df_row_order.append(row)
    df_raw = df_row_order.copy()
    df_raw.set_index('pos')

    # Add in aa information: columns
    for i, x in enumerate(list(df_raw.columns[:-1])):
        mesh_ax.text(i + 0.5, len(df_raw.index) + 1, x,
        horizontalalignment="center",
        verticalalignment="center",
        rotation = 90)
        if x == CODON2AA[WT_codon[1]][0]: WT_site2 = i

    # Add in amino acid property grouping information: columns
    new_AA_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(AA_GROUPS)):
        aaProperty = AA_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in Property2AA[aaProperty]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaProperty, cstart, cend)
            new_AA_GROUPS.append(newTuple)

    for codon, start, end in new_AA_GROUPS:
        if version == 2: codon = CODON321[codon]
        mesh_ax.text((end - start + 1) / 2 + start,
            row_number + 2.5, codon,
            horizontalalignment="center",
            verticalalignment="center")

        bar = Line2D([start + 0.125, end + 1 - 0.125],
                [row_number + 2, row_number + 2], color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)

    # Add in the deliminator for amino acid property groups: columns
        delimBar = Line2D([end + 1, end + 1],
                [0, len(df_raw.index) + 1], color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    # Add in aa information: row
    for i, x in enumerate(list(df_raw['pos'][::-1])):
        aa_text = [group[0] for group in CODON_GROUPS if group[1] == x][0] # determine the corresponding aa
        mesh_ax.text(col_number, i+1-0.5, aa_text, # Note that if -0.5, visually not perfect but this is the best achievable
                horizontalalignment="center",
                verticalalignment="center",)
        if CODON2AA[CodonList[int(x) - 1]][0] == CODON2AA[WT_codon[0]][0]:
            posInPropertyList = [i for i, x in enumerate(PropertyOrderList[::-1]) if x == CODON2AA[WT_codon[0]][0]][0]
            WT_site1 = posInPropertyList

    # Add in amino acid property grouping information: rows
    new_AA_GROUPS = []
    cstart = 0
    cend   = -1
    for i in range(0, len(AA_GROUPS)):
        aaProperty  = AA_GROUPS[i][0]
        count = 0
        for item in df.columns.values:
            if item in Property2AA[aaProperty]:
                count = count + 1
        if count == 1: cstart = cend = cend + 1
        else:
            gap = count - 1
            cstart = cend + 1
            cend = cstart + gap
        if cend >= cstart:
            newTuple = (aaProperty, cstart, cend)
            new_AA_GROUPS.append(newTuple)

    for codon, start, end in new_AA_GROUPS:
        if version == 2: codon = CODON321[codon]
        y = row_number - (end - start)/2 - start - 0.5
        mesh_ax.text(col_number+1.25, y, codon,
            horizontalalignment="left",
            verticalalignment="center")

        x = [col_number + 1, col_number + 1]
        y = [row_number - start - 1 + 0.875, row_number - end - 0.875]
        bar = Line2D(x, y, color="black")
        bar.set_clip_on(False)
        mesh_ax.add_line(bar)

    # Add in the deliminator for amino acid property groups: rows
        x = [0, col_number]
        y = [row_number - end - 1, row_number - end - 1]
        delimBar = Line2D(x, y, color="white")
        delimBar.set_clip_on(False)
        mesh_ax.add_line(delimBar)

    """
    # Legacy, add to figure title (colomn site) and to the right margin (row site)
    # Add in mutation site info onto the left y-axis:
        x = -2
        y = row_number / 3 * 2

        mesh_ax.text(x, y, "Row:\n" + str(int((wt_mut[0][0] + 2)/3)) + "-" + CODON2AA[WT_codon[0]][0] + "\n(" + WT_codon[0] + ")",
        horizontalalignment="center",
        verticalalignment="center")

        mesh_ax.text(x, y/2, "Column:\n" + str(int((wt_mut[0][1] + 2)/3)) + "-" + CODON2AA[WT_codon[1]][0] + "\n(" + WT_codon[1] + ")",
        horizontalalignment="center",
        verticalalignment="center")
    """

    # Add in first mutation site (per row) info onto figure right margin
    x = col_number + 5
    y = row_number / 2
    mesh_ax.text(x, y, "(Row: " + str(int((wt_mut[0][0] + 2)/3)) + "-" + CODON2AA[WT_codon[0]][0] + " (" + WT_codon[0] + "))",
    horizontalalignment="center",
    verticalalignment="center",
    rotation = 270)

    # Add in WT label onto corresponding cells:
    wt_aa     = [str(CODON2AA[codon][0]) for codon in WT_codon]
    if not WT_site1 == np.nan and not WT_site2 == np.nan:
        mesh_ax.add_patch(Circle((WT_site2 + 0.5, WT_site1 + 0.5), .1666,
        fill=True, facecolor="black",
        edgecolor="none", alpha=0.5))

    # Make the figure cleaner by removing ticks:
    mesh_ax.tick_params(bottom=False, left=False)
    mesh_ax.get_xaxis().set_visible(False)
    mesh_ax.get_yaxis().set_visible(False)
    cbar_ax.tick_params(bottom=False)

    # Add in SE if specified:
    if not df_se.shape == pd.DataFrame().shape: # Add in SE
        # Subset df_se to match df that might be truncated:
        for col in df_se.columns:
            if not col in df.columns: df_se.drop(columns=col, inplace=True)
        df_se = df_se[df_se['pos'].isin(df['pos'])]

        # Below is to find max se value:
        se_df = df_se.drop(columns=['pos'])
        tem = sorted(se_df.values.reshape(1,-1).tolist()[0])
        tem = [x for x in tem if not np.isnan(x)]
        se_max = tem[-1]

        for row in range(len(df_se['pos'])):
            for col in range(len(df_se.columns)-1):
                se_value = df_se.iloc[row, col]
                corner_dist = (se_max - se_value)/(2 * se_max)
                corner_dist = (1 - se_value) / 2
                diag = Line2D([col + corner_dist, col + 1 - corner_dist],
                        [row + corner_dist, row + 1 - corner_dist], color="black")

                if se_value > 0.02 and df_raw.iloc[row, col] != 1000: # se_value below 0.02 will not be displayed so as totally depleted ones
                    mesh_ax.add_line(diag)

    pylab.savefig(outfilename)
