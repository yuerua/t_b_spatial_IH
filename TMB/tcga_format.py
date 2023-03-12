#Adopted from https://github.com/zengxiaofei/blog_scripts/tree/master/01.archives_9
#Format the 'Tumor_Sample_Barcode' column in maf files to Project-TSS-Participant
#Usage: python tcga_format.py LUSC_549.maf > LUSC_549_sort.maf

# Author: Xiaofei Zeng
# Email: xiaofei_zeng@whu.edu.cn
# Created Time: 2019-05-15 11:23
# Python 2.7

from __future__ import print_function
import sys
import gzip

def handle_maf(maf_file):
    count_dict = {}
    if maf_file.endswith('.gz'):
        fopen = gzip.open
    else:
        fopen = open
    with fopen(maf_file) as f:
        for line in f:
            if line.startswith('#'):
                print(line, end='')
                continue
            ls = line.split('\t')
            # header
            if line.startswith('Hugo'):
                index = ls.index('Tumor_Sample_Barcode')
                print(line, end='')
                continue
            # content
            prefix = ls[index].rsplit('-', 4)[0]
            new_line = '\t'.join(ls[:index] + [prefix] + ls[index+1:])
            print(new_line, end='')

def main():
    handle_maf(sys.argv[1])

if __name__ == '__main__':
    main()