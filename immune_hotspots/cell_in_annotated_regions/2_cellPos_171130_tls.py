#Assign annotated regions to T, B cells
import numpy as np
import pickle
import os
import math
import pandas as pd
import glob

cws_folder = os.path.normpath('/Volumes/COPAINGE/hzhang/21_01_06_tracerx_serial_3_pairs/serial_3pairs/data/cws_reg')
wsi_files = sorted(glob.glob(os.path.join(cws_folder, '*_Tcell.ndpi')))
results_dir = os.path.normpath('/Users/hzhang/Documents/project/sum/final/TLS_annotation/Region/T_region')
output_dir = os.path.normpath('/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/T')

if not os.path.isdir(output_dir):
    os.makedirs(output_dir, exist_ok=True)


for wsi_n in range(0, len(wsi_files)):
    current_wsi = wsi_files[wsi_n]
    param = pickle.load(open(os.path.join(current_wsi, 'param.p'), 'rb'))
    print(param)
    current_results_dir = os.path.join(results_dir, os.path.basename(current_wsi))

    slide_dimension = np.array(param['slide_dimension']) / param['rescale']
    slide_h = slide_dimension[1]
    slide_w = slide_dimension[0]
    cws_read_size = param['cws_read_size']
    cws_h = cws_read_size[0]
    cws_w = cws_read_size[1]
    divisor = np.float64(16)

    # Initialize Pandas Data Frame
    cellPos = pd.DataFrame(columns=['class', 'x', 'y', 'region'])
    iter_tot_tiles = 0
    for h in range(int(math.ceil((slide_h - cws_h) / cws_h + 1))):
        for w in range(int(math.ceil((slide_w - cws_w) / cws_w + 1))):
            # print('Processing Da_' + str(iter_tot_tiles))
            start_h = h * cws_h
            start_w = w * cws_w

            if os.path.isfile(os.path.join(current_results_dir, 'Da' + str(iter_tot_tiles) + '.csv')):
                csv = pd.read_csv(os.path.join(current_results_dir, 'Da' + str(iter_tot_tiles) + '.csv'))
                # print(csv)
                csv.columns = ['class', 'x', 'y', 'region', 'region_idx', 'join_idx']
                csv.x = csv.x + start_w
                csv.y = csv.y + start_h
                # detection = np.divide(np.float64(detection), divisor)
                cellPos = cellPos.append(csv)

            iter_tot_tiles += 1

    cellPos.x = np.round(np.divide(np.float64(cellPos.x), divisor)).astype('int')
    cellPos.y = np.round(np.divide(np.float64(cellPos.y), divisor)).astype('int')
    cellPos.loc[cellPos.x == 0, 'x'] = 1
    cellPos.loc[cellPos.y == 0, 'y'] = 1
    cellPos.to_csv(os.path.join(output_dir, os.path.splitext(os.path.basename(current_wsi))[0] + '_cellPos.csv'),
                   index=False)
