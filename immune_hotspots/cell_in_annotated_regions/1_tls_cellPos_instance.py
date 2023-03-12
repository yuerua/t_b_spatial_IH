import os
import numpy as np
from glob import glob
import scipy.io as sio
import pandas as pd
import cv2
import re
import scipy.ndimage as ndi
import pickle

# import matplotlib.pyplot as plt
# import matplotlib
# import platform
# if platform.system() == 'Darwin':
#     matplotlib.use('MacOSX')

# TODO: Pick cells in the same region, give them a unique label
# TODO: merge continues regions into one

np.random.seed(655291)

def put_markers(csv, im):
    """
    Put markers on img
    :param csv: cell csv
    :param im: Da*
    :param pre_labels: labels, list
    :return: marked img
    """
    im_c = im.copy()
    r = 4 #change radius if needed
    for i, row in csv.iterrows():
        label_c = row['V1']
        x = row['V2']
        y = row['V3']

        # if label_c == 'None':
        color_c = (0,255,0) #default as green
        # else:
        #     h = color_dict[label_c].lstrip('#')
        #     color_c = tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))

        cv2.circle(im_c, (x, y), r, color=color_c, thickness=-1, lineType=cv2.LINE_AA)
    return im_c

def cell_region_label_instance(mask, img_name, cell_csv):
    """
    Attach region label to cell coordinates
    :param mask: binary mask, region type split by channels
    :param img_name: Da*
    :return: cell csv matched with region, region idx
    """
    inst_mask, num_inst = ndi.label(mask)

    cell_csv['region'] = None
    cell_csv['region_idx'] = None
    try:
        for c in [c for c in range(inst_mask.shape[-1]) if np.max(inst_mask[..., c]) > 0]:
            for _, row in cell_csv.iterrows():
                if mask[row['V3'], row['V2'],c] > 0:
                    cell_csv.loc[_,'region'] = region_type[c]
                    cell_csv.loc[_, 'region_idx'] = img_name + '_i' + str(int(inst_mask[row['V3'], row['V2'],c]))

    except Exception as e:
        print(e)
        pass

    return cell_csv

def get_cws_w_h(path_to_param, cws_i):
    """
    Transform cws index to width and height
    :param path_to_param: path to param.p
    :param cws_i: index of Da*
    :return: h_i: nth row, w_i: nth column
    """
    #cws_i = int(re.search(r'\d+', cws_name).group())
    cws_i = int(cws_i)
    param = pickle.load(open(path_to_param, 'rb'))
    slide_dimension = np.array(param['slide_dimension']) / param['rescale']

    slide_w, slide_h = slide_dimension
    cws_w, cws_h = param['cws_read_size']

    divisor_w = np.ceil(slide_w / cws_w)
    divisor_h = np.ceil(slide_h / cws_h)

    # exact pixel
    # h_i = np.floor(cws_i / divisor_w) * cws_h
    # w_i = (cws_i - h_i / cws_h * divisor_w) * cws_w
    # index only
    h_i = int(np.floor(cws_i / divisor_w))
    w_i = int(cws_i - h_i * divisor_w)

    # img = cv2.imread(path)
    # img_all[h_i : h_i + int(img.shape[0]), w_i : w_i + int(img.shape[1]),:] = img
    # eg. (17, 16): the 18th row, the 17th column
    return h_i, w_i

def get_instance(mask, img_name):
    """
    Find the min&max x,y, check if they equal to 0/img width
    # If an annotation is integrate, it shouldn't touch the boarder
    :param mask: binary mask, region type split by channels
    :param img_name: Da*
    :return: dataframe of region properties
    # columns - instance index; channel; patch index; touch top; right; bottom; left
    """
    img_idx = int(re.search(r'\d+', img_name).group())
    # get instance mask for each channel
    inst_mask, num_inst = ndi.label(mask)
    # # change instance label of each channel to continuous value
    # for i in range(inst_mask.shape[2]):
    #     mask_c = inst_mask[...,i]
    #     inst = np.unique(mask_c)
    #     for idx in range(len(inst)):
    #         mask_c[mask_c==inst[idx]] = idx

    # generate the variable list
    channel, patch_idx, top, right, bottom, left = ([] for i in range(6))

    for i in range(num_inst):
        X, Y, C = np.where(inst_mask == i + 1)
        patch_idx.append(int(img_idx))
        channel.append(int(np.unique(C)[0]))
        if np.min(Y) == 0:
            # left.append((min(X[Y == 0]), max(X[Y == 0])))
            left.append(extract_border_ends(X, Y, 0))
        else:
            left.append(None)

        if np.max(Y) == mask.shape[0] - 1:
            # right.append((min(X[Y == mask.shape[0] - 1]), max(X[Y == mask.shape[0] - 1])))
            right.append(extract_border_ends(X, Y, mask.shape[0] - 1))
        else:
            right.append(None)

        if np.min(X) == 0:
            # top.append((min(Y[X == 0]), max(Y[X == 0])))
            top.append(extract_border_ends(Y, X, 0))
        else:
            top.append(None)

        if np.max(X) == mask.shape[1] - 1:
            # bottom.append((min(Y[X == mask.shape[1] - 1]), max(Y[X == mask.shape[1] - 1])))
            bottom.append(extract_border_ends(Y, X, mask.shape[1] - 1))
        else:
            bottom.append(None)

    inst_idx = [int(i + 1) for i in range(num_inst)]

    df_patch = pd.DataFrame({"patch_idx": patch_idx,
                             "inst_idx": inst_idx,
                             "channel": channel,
                             "top": top,
                             "right": right,
                             "bottom": bottom,
                             "left": left})

    return df_patch

def get_join_idx(df_all, path_to_param):
    """
    Label regions spanning on adjacent patches
    # eg. join_idx = Da408_i1_Da407_i2
    :param df_all: dataframe of regions from the whole slide
    :param path_to_param: path to param.p
    :return: dataframe of regions with join_idx
    # join_idx: list of idx of connected regions
    """
    df_all["row_i"] = None
    df_all["col_i"] = None
    # allow some mismatch
    err = 10

    for i, row in df_all.iterrows():
        patch_idx = row['patch_idx']
        row_i, col_i = get_cws_w_h(path_to_param, patch_idx)
        df_all.loc[i, "row_i"] = row_i
        df_all.loc[i, "col_i"] = col_i
    #     # almost no single region
    #     df_all["region_type"] = pd.np.where(df_all[['top', 'right', 'bottom', 'left']].eq(None).all(), 'Single', None)
    #     print(np.where(pd.notnull(row[['top', 'right', 'bottom', 'left']])))

    df_all["join_idx"] = None
    df_all["join_idx"] = df_all["join_idx"].astype('object')
    df_right = df_all[df_all.right.notnull()]
    df_bottom = df_all[df_all.bottom.notnull()]

    for i, row in df_all.iterrows():
        #touch on left-right
        if row["left"] is not None:
            suspec_idxs = []
            if not df_right.empty:
                for j, s_row in df_right.iterrows():
                    for ends_i in row["left"]:
                        for ends_j in s_row["right"]:
                            if (ends_i[0] - err <= ends_j[0] <= ends_i[0] + err) or \
                                             (ends_i[1] - err <= ends_j[1] <= ends_i[1] + err):
                                suspec_idxs.append(j)

            if suspec_idxs:
                for s_i in suspec_idxs:
                    if (df_all.loc[s_i, "col_i"] == row["col_i"] - 1 and
                        df_all.loc[s_i, "row_i"] == row["row_i"] and
                        df_all.loc[s_i, "channel"] == row["channel"]):

                        join_idx = ["Da" + str(int(row["patch_idx"])) + "_i" + str(int(row["inst_idx"])),
                                    "Da" + str(int(df_all.loc[s_i, "patch_idx"])) + "_i" + str(int(df_all.loc[s_i, "inst_idx"]))]
                        if df_all.loc[i, "join_idx"] is None:
                            df_all.at[i, "join_idx"] = join_idx
                        else:
                            df_all.at[i, "join_idx"] += join_idx

                        if df_all.loc[s_i, "join_idx"] is None:
                            df_all.at[s_i, "join_idx"] = join_idx
                        else:
                            df_all.at[s_i, "join_idx"] += join_idx

        # touch on bottom-top
        if row["top"] is not None:
            suspec_idxs = []
            if not df_bottom.empty:
                for j, s_row in df_bottom.iterrows():
                    for ends_i in row["top"]:
                        for ends_j in s_row["bottom"]:
                            if (ends_i[0] - err <= ends_j[0] <= ends_i[0] + err) or \
                                             (ends_i[1] - err <= ends_j[1] <= ends_i[1] + err):
                                suspec_idxs.append(j)

            if suspec_idxs:
                for s_i in suspec_idxs:
                    if (df_all.loc[s_i, "row_i"] == row["row_i"] - 1 and
                            df_all.loc[s_i, "col_i"] == row["col_i"] and
                            df_all.loc[s_i, "channel"] == row["channel"]):

                        join_idx = ["Da" + str(int(row["patch_idx"])) + "_i" + str(int(row["inst_idx"])),
                                    "Da" + str(int(df_all.loc[s_i, "patch_idx"])) + "_i" + \
                                   str(int(df_all.loc[s_i, "inst_idx"]))]

                        if df_all.loc[i, "join_idx"] is None:
                            df_all.at[i, "join_idx"] = join_idx
                        else:
                            df_all.at[i, "join_idx"] += join_idx

                        if df_all.loc[s_i, "join_idx"] is None:
                            df_all.at[s_i, "join_idx"] = join_idx
                        else:
                            df_all.at[s_i, "join_idx"] += join_idx

    return df_all

def get_common_idx(r_i):
    """
    Iterate through a list of list, filter the common list
    # https://stackoverflow.com/questions/42036188/merging-tuples-if-they-have-one-common-element
    :param r_i: list of list or tuple
    :return: list of sets of common elements
    """
    data = list(map(set, r_i))

    oldlen = len(data) + 1
    while len(data) < oldlen:
        oldlen = len(data)
        for i in range(len(data)):
            for j in range(i + 1, len(data)):
                if len(data[i] & data[j]):
                    data[i] = data[i] | data[j]
                    data[j] = set()
        data = [data[i] for i in range(len(data)) if data[i] != set()]
    return data

def extract_border_ends(MC, FC, border):
    """
    Get ends of border, return in tuple of tuples
    #Credit: https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
    :param MC: array of moving axis (eg. to get the left border, MC = Y)
    :param FC: array of fixed axis (eg. to get the left border, FC = X)
    :param border: int, the border value of X or Y (eg. 0 or 1999)
    :return: tuple of ends for the moving axis eg. ((0,3), (19, 1999))
    """
    nums = MC[FC == border]

    nums = sorted(set(nums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 1 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    border_ends = tuple(zip(edges, edges))

    return border_ends


def run():
    """
    Label cells with idx of freehand regions
    :return: csv columns - V1 V2 V3 region region_idx join_idx
    """
    slides = [d for d in os.listdir(annotation_dir) if d.endswith('.ndpi')]

    for slide in slides:
        if os.path.exists(os.path.join(save_path, slide)) is True:
            print("Already processed", slide)
            continue
        print(slide)
        slide_no = str(re.findall(r'\d+', slide)[0])

        file_paths = glob(os.path.join(annotation_dir, slide, '*.mat'))

        param_slide = os.path.join(cws_path, prefix + slide_no + ext, 'param.p')
        region_list = []#store region dfs

        if os.path.exists(os.path.join(save_path, prefix + slide_no + ext)) is False:
            os.makedirs(os.path.join(save_path, prefix + slide_no + ext))

        if os.path.exists(os.path.join(save_path, prefix + slide_no + ext, 'tmp')) is False:
            os.makedirs(os.path.join(save_path, prefix + slide_no + ext, 'tmp'))

        print("Processing regions")
        for file_path in file_paths:

            img_name = os.path.splitext(os.path.basename(file_path))[0]
            if os.path.exists(os.path.join(save_path, slide, img_name + '.csv')) is True:
                print("Already processed", img_name)
                continue

            print(img_name)

            f = sio.loadmat(file_path)
            mask = f['map']
            # first label cells with region id, and save to tmp dir
            img_csv = os.path.join(csv_path, prefix + slide_no + ext, img_name + '.csv')

            try:
                cell_csv = pd.read_csv(img_csv)
            except:
                print("Doesn't exit for", img_csv)
                continue

            if divisor:
                cell_csv.V2 = np.round(np.divide(np.float64(cell_csv.V2), divisor)).astype('int')
                cell_csv.V3 = np.round(np.divide(np.float64(cell_csv.V3), divisor)).astype('int')
                cell_csv.loc[cell_csv.V2 == 2000, 'V2'] = 1999
                cell_csv.loc[cell_csv.V3 == 2000, 'V3'] = 1999

            cell_csv = cell_region_label_instance(mask, img_name, cell_csv)
            cell_csv.to_csv(os.path.join(save_path, prefix + slide_no + ext, 'tmp', img_name + '.csv'), index=False)

            if check_img:

                if os.path.exists(os.path.join(save_path, 'img', prefix + slide_no + ext)) is False:
                    os.makedirs(os.path.join(save_path, 'img', prefix + slide_no + ext))

                img_ori_f = os.path.join(annotated_tiles_dir, slide, img_name + '.jpg')
                img_ori = cv2.imread(img_ori_f)
                img_marked = put_markers(cell_csv[cell_csv.region.notnull()], img_ori)
                cv2.imwrite(os.path.join(save_path, 'img', prefix + slide_no + ext, img_name + '.jpg'), img_marked)

            # then deal with regions, find index of each region
            df_patch = get_instance(mask, img_name)
            region_list.append(df_patch)

        # convert list of dfs to single df
        df_region = region_list[0]
        for df_ in region_list[1:]:
            df_region = pd.concat([df_region, df_], ignore_index=True)

        # label the continuous regions on dif patches
        df_region = get_join_idx(df_region, param_slide)
        df_region["patch_idx"] = df_region["patch_idx"].astype(int)
        df_region["inst_idx"] = df_region["inst_idx"].astype(int)
        df_region["region_idx"] = "Da" + df_region["patch_idx"].astype(str) + "_i" + df_region["inst_idx"].astype(str)

        r_i = df_region["join_idx"].tolist()
        single_no = sum(1 for i in r_i if i is None)
        r_i = [x for x in r_i if x is not None]
        data = get_common_idx(r_i)

        for i, row in df_region.iterrows():
            if row["join_idx"] is not None:
                for common in data:
                    # sort the string
                    common = sorted(list(common))
                    if any(item in common for item in row["join_idx"]):
                        df_region.at[i, "join_idx"] = "_".join(common)
            else:
                df_region.loc[i, "join_idx"] = (df_region.loc[i, "region_idx"])

        df_region.to_csv(os.path.join(save_path, prefix + slide_no + ext + '_region.csv'), index=False)
        print("Total no. of region:", len(data)+single_no)

        print('Joining cell csv')
        for file_path in file_paths:
            img_name = os.path.splitext(os.path.basename(file_path))[0]
            img_csv = os.path.join(save_path, prefix + slide_no + ext, 'tmp', img_name + '.csv')
            try:
                cell_csv = pd.read_csv(img_csv)
            except:
                continue
            # join the two df and write
            cell_csv['region_idx'] = cell_csv['region_idx'].astype(object)
            df_join = pd.merge(cell_csv, df_region[["region_idx", "join_idx"]], how='left', on='region_idx')

            df_join.to_csv(os.path.join(save_path, prefix + slide_no + ext, img_name + '.csv'), index=False)
            os.remove(os.path.join(save_path, prefix + slide_no + ext, 'tmp', img_name + '.csv'))

        try:
            os.rmdir(os.path.join(save_path, prefix + slide_no + ext, 'tmp'))
        except Exception as e:
            print(e)
            pass

        # Copy other csv to dir
        print("Copying other csv")
        all_files = glob(os.path.join(csv_path, prefix + slide_no + ext, 'Da*'))
        for src in all_files:
            csv_name = os.path.basename(src)
            dst = os.path.join(save_path, prefix + slide_no + ext, csv_name)
            if os.path.exists(dst) is False:
                csv_ori = pd.read_csv(src)
                if divisor:
                    csv_ori.V2 = np.round(np.divide(np.float64(csv_ori.V2), divisor)).astype('int')
                    csv_ori.V3 = np.round(np.divide(np.float64(csv_ori.V3), divisor)).astype('int')
                    cell_csv.loc[cell_csv.V2 == 2000, 'V2'] = 1999
                    cell_csv.loc[cell_csv.V3 == 2000, 'V3'] = 1999

                csv_ori['region'] = None
                csv_ori['region_idx'] = None
                csv_ori['join_idx'] = None
                csv_ori.to_csv(dst, index=False)


if __name__=='__main__':
    annotation_dir = '/Users/hzhang/Documents/project/sum/final/Misc/freehand/freehandlabels'#mat files
    annotated_tiles_dir = '/Users/hzhang/Documents/project/sum/final/Misc/freehand/AnnotatedTiles'
    save_path = '/Users/hzhang/Documents/project/sum/final/TLS_annotation/Region/T_region_new'
    csv_path = '/Users/hzhang/Documents/project/sum/CellPos_all/T' # csv for classification
    # csv_path = '/Volumes/proj4/TracerX/SqCCs_fi_20x/serial_3pairs/results/classification/csv'
    cws_path = '/Volumes/COPAINGE/hzhang/21_01_06_tracerx_serial_3_pairs/serial_3pairs/data/cws_reg' # only param.p is required
    divisor = ''#np.float64(2022/2000)
    check_img = False
    # slide name: prefix + idx + ext
    # divisor = 2022/2000
    prefix = 'Reg_'
    ext = '_Tcell.ndpi'
    # region type: channel
    region_type = {
        0: 'lag',
        1: 'tls'
    }

    run()

