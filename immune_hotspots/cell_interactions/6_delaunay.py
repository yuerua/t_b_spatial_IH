"""Construct Delaunay Triangulation"""
import os
# import networkx as nx
# from PIL import Image
# from shutil import copyfile
# import cv2
import numpy as np
from scipy.spatial import Delaunay, delaunay_plot_2d
import pandas as pd
# import matplotlib.pyplot as plt
# from glob import glob


#convert to edges
def less_first(a, b):
    return [a,b] if a < b else [b,a]

def delaunay2edges(df):

    points = np.array(list(zip(df['x'], df['y'])))

    tri = Delaunay(points)

    list_of_edges = []

    for triangle in tri.simplices:
        for e1, e2 in [[0,1],[1,2],[2,0]]: # for all edges of triangle
            list_of_edges.append(less_first(triangle[e1],triangle[e2])) # always lesser index first

    array_of_edges = np.unique(list_of_edges, axis=0) # remove duplicates

    #length of edges
    list_of_lengths = []

    for p1,p2 in array_of_edges:
        x1, y1 = tri.points[p1]
        x2, y2 = tri.points[p2]
        list_of_lengths.append((x1-x2)**2 + (y1-y2)**2)

    array_of_lengths = np.sqrt(np.array(list_of_lengths))

    # p_all = array_of_lengths.shape[0]

    if pixel_th:
        # print(array_of_lengths.shape, array_of_lengths[array_of_lengths <= pixel_th].shape)
        array_of_edges = np.take(array_of_edges, np.where(array_of_lengths <= pixel_th)[0], 0)
        array_of_lengths = np.take(array_of_lengths, np.where(array_of_lengths <= pixel_th)[0], 0)

    # if p_all - array_of_lengths.shape[0] > 0:
    #     print(p_all - array_of_lengths.shape[0])

    #cell types of edges
    cell_edge_dic = {key: 0 for key in cell_dic_columns}
    for p1, p2 in array_of_edges:
        c1 = df['class'][p1]
        c2 = df['class'][p2]
        label = '-'.join(sorted([c1, c2]))
        cell_edge_dic[label] += 1

    cell_num = np.unique(df['class'], return_counts=True)
    for i in range(len(cell_num[0])):
        cell_edge_dic[cell_num[0][i]] = cell_num[1][i]



    return array_of_edges, array_of_lengths, cell_edge_dic


if __name__ == "__main__":

    clump_dir = "/Users/hzhang/Documents/project/sum/final/clump_l_cl"
    clump_ext = "_hs_sum_clump_df.csv"

    cellPosDir = "/Users/hzhang/Documents/project/sum/final/TLS_annotation/CellPos/CellPos_region_and_hs"
    cellPos_ext = ".csv"

    save_dir = "/Users/hzhang/Documents/project/sum/final/graph/clump_l_cl_TLS_LAG_th"

    res_40x = 228
    res_ss1 = res_40x * 2 * 16 / 1000
    distance_th = 250
    pixel_th = distance_th / res_ss1


    select_cell_types = sorted(['cd20', 'cd20cxcr5', 'cd4', 'cd79bCoexp', 'cd8', 'foxp3', 'p40'])
    select_hs = ['cl', 'l']
    cell_interaction = ["-".join(sorted([i, j])) for i in select_cell_types for j in select_cell_types]
    cell_interaction = list(set(cell_interaction))

    cell_dic_columns = cell_interaction + select_cell_types

    cellPos_f = [f for f in os.listdir(cellPosDir) if f.endswith(cellPos_ext)]
    print(cellPos_f)

    if os.path.exists(save_dir) is False:
        os.makedirs(save_dir)

    for f in cellPos_f:
        fname = f.replace(cellPos_ext, '')
        print('Processing', fname)
        cellPos_df_ = pd.read_csv(os.path.join(cellPosDir, f))
        clump_df = pd.read_csv(os.path.join(clump_dir, fname + clump_ext))
        #area of clumps
        clump_df_count = clump_df.groupby(['hs', 'clump_id'])['hs'].count().reset_index(name="clump_area")
        #
        # cellPos_df_ = pd.read_csv("/Users/hzhang/Documents/project/sum/TLS_annotation/region_morisita/CellPos_region_and_hs/82110_Cellpos_add_hs.csv")
        # clump_df = pd.read_csv("/Users/hzhang/Documents/report/lusc_b/scripts/ITLR/results/clump_l_cl_v1/82110_cl_df_hs_sum_clump_df.csv")

        #use the reassigned row index to avoid conflicts
        cellPos_df_.drop('Unnamed: 0', axis=1, inplace=True)
        clump_df.drop('Unnamed: 0', axis=1, inplace=True)

        cellPos_df = cellPos_df_[(cellPos_df_['class'].isin(select_cell_types))]
        cellPos_df = cellPos_df[(cellPos_df['hs'].isin(select_hs))]

        #without TLS
        # cellPos_df = cellPos_df[cellPos_df['region'] != 'tls']
        # cellPos_df = cellPos_df[~cellPos_df['region'].isin(['tls', 'lag'])]
        cellPos_df = cellPos_df[cellPos_df['region'].isin(['tls', 'lag'])]



        df_all = pd.merge(cellPos_df, clump_df, on = ['x.s', 'y.s', 'hs'])
        df_all.reset_index(drop=True, inplace=True)

        slide_df = pd.DataFrame(columns= cell_dic_columns + ['hs', 'clump_id', 'slide'])

        slide_dic = {key: [] for key in ['hs', 'clump_id', 'slide', 'edges', 'lengths', 'cell_edge_dic']}

        for hs in select_hs:
            print(hs)
            df_hs = df_all[df_all['hs'] == hs]
            clump_id = np.unique(df_hs['clump_id'])
            for clump in clump_id:
                df = df_hs[df_hs['clump_id'] == clump]
                df.reset_index(drop=True, inplace=True)

                if len(df) <= 3:
                    continue

                edges, lengths, cell_edge_dic = delaunay2edges(df)

                #slide_df: cell interaction counts
                cell_edge_dic['hs'] = hs
                cell_edge_dic['clump_id'] = clump
                cell_edge_dic['slide'] = fname

                slide_df = slide_df.append(cell_edge_dic, ignore_index=True)

                #slide_dic: other variables
                slide_dic['hs'].append(hs)
                slide_dic['clump_id'].append(clump_id)
                slide_dic['slide'].append(fname)
                slide_dic['edges'].append(edges)
                slide_dic['lengths'].append(lengths)
                slide_dic['cell_edge_dic'].append(cell_edge_dic)

        slide_df = pd.merge(slide_df, clump_df_count, on = ['hs', 'clump_id'])
        slide_df.to_csv(os.path.join(save_dir, fname + "_cell_interaction.csv"))
        np.save(os.path.join(save_dir, fname + ".npy"), slide_dic, allow_pickle=True)