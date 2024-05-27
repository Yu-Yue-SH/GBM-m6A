"""
decompress the tar file and the gz file
"""


import tarfile
import gzip


def untar(tar_file, output_dic):
    with tarfile.open(tar_file, 'r') as f:
        f.extractall(output_dic)


def ungz(gz_file, output_file):
    with gzip.open(gz_file, 'rb') as f:
        with open(output_file, 'wb') as out:
            out.write(f.read())


if __name__ == '__main__':
    from pathlib import Path

    # untar raw data
    raw_path = '../data/GSE165037_RAW.tar'
    dic_path = '../data/GSE165037_RAW'
    untar(raw_path, dic_path)

    # ungz all files
    ungz_dic = '../data/GSE165037_UNGZ'
    Path(ungz_dic).mkdir(exist_ok=True)
    for i in Path(dic_path).iterdir():
        out_path = Path(ungz_dic, i.stem)
        ungz(i, out_path)
