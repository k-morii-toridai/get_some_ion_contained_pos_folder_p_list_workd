import os
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
from tqdm import tqdm

from package_bond_search_algorithm.algolithm_bond_search_for_BrO3 import concat_filter
from package_bond_search_algorithm.package_file_conversion.nnlist2df import nnlist2df


def wrap_nnlist2df_and_concat_filter(nnlist_p):
    """
    This func() executes two func()s.: nnlistdf() and concat_filter().

    Usage:
    ------
    wrap_nnlist2df_and_concat_filter(nnlist_p=nnlist_p)

    Parameter:
    -----------
    nnlist_p: str or PosixPath

    Return:
    -------
    bool: True or False
    """
    df_nnlist = nnlist2df(nnlist_p)
    bool_ = concat_filter(df_nnlist=df_nnlist)

    return bool_


# ターゲットとなるイオンの元素種を含むPOSCARのディレクトリパス一覧を取得
target_npy_p = '../get_some_speceis_existed_poscar_path_list/Br_O_existed_poscar_folder_path_list.npy'
Br_O_existed_poscar_folder_path_list = np.load(target_npy_p, allow_pickle=True)
print(f"len(Br_O_existed_poscar_folder_path_list): {len(Br_O_existed_poscar_folder_path_list)}")

# ターゲットとなるイオンの元素種を含むPOSCAR.nnlistのディレクトリパス一覧を取得
add_str = '/nnlist_5/POSCAR.nnlist'
Br_O_existed_poscar_nnlist_path_list = [str(p) + add_str for p in Br_O_existed_poscar_folder_path_list]

# 並列化
pp = Pool(cpu_count() - 1)
total = len(Br_O_existed_poscar_nnlist_path_list)
try:
    bool_ion_contain = list(tqdm(pp.imap(wrap_nnlist2df_and_concat_filter, Br_O_existed_poscar_nnlist_path_list), total=total))
finally:
    pp.close()
    pp.join()

# apply filter(:bool_ion_contain) to Br_O_existed_poscar_nnlist_path_list
ion_contained_poscar_nnlist_path_list = np.array(Br_O_existed_poscar_nnlist_path_list)[bool_ion_contain]
# get ion contained folder path list
ion_contained_poscar_folder_path_list = [Path(os.path.split(Path(os.path.split(Path(p))[0]))[0]) for p in ion_contained_poscar_nnlist_path_list]
# save ion contained folder path list as .npy
saved_npy_fname = 'BrO3_contained_poscar_folder_path_list.npy'
np.save(saved_npy_fname, ion_contained_poscar_folder_path_list)
# print parcentage
print(f"len(ion_contained_poscar_folder_path_list)/len(Br_O_existed_poscar_folder_path_list) :\
{len(ion_contained_poscar_folder_path_list)}/{total}")
