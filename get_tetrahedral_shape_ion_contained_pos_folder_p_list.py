import os
import argparse
from pathlib import Path
from multiprocessing import Pool, cpu_count

import numpy as np
from tqdm import tqdm

from package_bond_search_algorithm.algolithm_bond_search_for_tetrahedral_shape import concat_filter
from package_bond_search_algorithm.package_file_conversion.nnlist2df import nnlist2df


# コマンドライン引数を受け取る
parser = argparse.ArgumentParser(description='This script takes five arguments: arg1, arg2, arg3, arg4, arg5 and arg6.',
                                 usage='%(prog)s <arg1> <arg2> <arg3> <arg4> <arg5> <arg6> \
                                 \nexample: python3 %(prog)s NH4 N H 0.82 1.24 ../get_some_speceis_existed_poscar_path_list/N_H_existed_poscar_folder_path_list.npy')
parser.add_argument('arg1', help='target_ion_name: NH4')
parser.add_argument('arg2', help='central_atom_symbol: N')
parser.add_argument('arg3', help='neighboring_atom_symbol: H')
parser.add_argument('arg4', help='bond_length_lower_end: 0.82')
parser.add_argument('arg5', help='bond_length_upper_end: 1.24')
parser.add_argument('arg6', help='npy_file_path: ../get_some_speceis_existed_poscar_path_list/N_H_existed_poscar_folder_path_list.npy')
args = parser.parse_args()
target_ion_name = args.arg1
central_atom_symbol = args.arg2
neighboring_atom_symbol = args.arg3
bond_length_lower_end = args.arg4
bond_length_upper_end = args.arg5
target_npy_p = args.arg6
print(f'target_ion_name: {target_ion_name}')
print(f'central_atom_symbol: {central_atom_symbol}')
print(f'neighboring_atom_symbol: {neighboring_atom_symbol}')
print(f'bond_length_lower_end: {bond_length_lower_end}')
print(f'bond_length_upper_end: {bond_length_upper_end}')
print(f'target_npy_p: {target_npy_p}')
print(f'os.path.exists(target_npy_p): {os.path.exists(target_npy_p)}')

some_species_existed_poscar_folder_path_list = np.load(target_npy_p, allow_pickle=True)
print(f"len(some_species_existed_poscar_folder_path_list): {len(some_species_existed_poscar_folder_path_list)}")
# テスト用
# some_species_existed_poscar_folder_path_list = some_species_existed_poscar_folder_path_list[0:100]


def nnlist2df_and_concat_filter(nnlist_path, central_atom_symbol, neighboring_atom_symbol, bond_length_lower_end, bond_length_upper_end):
    df_nnlist = nnlist2df(nnlist_path=nnlist_path)
    return_bool = concat_filter(df_nnlist=df_nnlist,
                                central_atom_symbol=central_atom_symbol,
                                neighboring_atom_symbol=neighboring_atom_symbol,
                                bond_length_lower_end=float(bond_length_lower_end),
                                bond_length_upper_end=float(bond_length_upper_end))
    return return_bool


def wrap_nnlist2df_and_concat_filter(args):
    return nnlist2df_and_concat_filter(*args)


def mk_job_args(some_species_existed_nnlist_folder_p_list, central_atom_symbol, neighboring_atom_symbol, bond_length_lower_end, bond_length_upper_end):
    # ターゲットとなるイオンの元素種を含むPOSCAR.nnlistのディレクトリパス一覧を取得
    add_str = '/nnlist_5/POSCAR.nnlist'
    some_species_existed_nnlist_p_list = [str(p) + add_str for p in some_species_existed_nnlist_folder_p_list]
    number_of_nnlist = len(some_species_existed_nnlist_folder_p_list)
    central_atom_symbol_list = [central_atom_symbol for i in range(number_of_nnlist)]
    neighboring_atom_symbol_list = [neighboring_atom_symbol for i in range(number_of_nnlist)]
    bond_length_lower_end_list = [bond_length_lower_end for i in range(number_of_nnlist)]
    bond_length_upper_end_list = [bond_length_upper_end for i in range(number_of_nnlist)]
    job_args = zip(some_species_existed_nnlist_p_list,
                   central_atom_symbol_list,
                   neighboring_atom_symbol_list,
                   bond_length_lower_end_list,
                   bond_length_upper_end_list)
    return job_args


job_args = mk_job_args(some_species_existed_poscar_folder_path_list,
                       central_atom_symbol=central_atom_symbol,
                       neighboring_atom_symbol=neighboring_atom_symbol,
                       bond_length_lower_end=bond_length_lower_end,
                       bond_length_upper_end=bond_length_upper_end)

# 並列化
pp = Pool(cpu_count() - 1)
total = len(some_species_existed_poscar_folder_path_list)
try:
    bool_ion_contain = list(tqdm(pp.imap(wrap_nnlist2df_and_concat_filter, job_args), total=total))
finally:
    pp.close()
    pp.join()

# 保存
# apply filter(:bool_ion_contain) to N_H_existed_poscar_nnlist_path_list
ion_contained_poscar_nnlist_path_list = np.array(some_species_existed_poscar_folder_path_list)[bool_ion_contain]
# save ion contained folder path list as .npy
saved_npy_fname = f'{target_ion_name}_contained_poscar_folder_path_list_ver2.npy'
np.save(saved_npy_fname, ion_contained_poscar_nnlist_path_list)
# print parcentage
print(f"len({target_ion_name}_contained_poscar_folder_path_list)/len({os.path.split(Path(target_npy_p))[1]}) :\
{len(ion_contained_poscar_nnlist_path_list)}/{total}")
