import argparse


# コマンドライン引数を受け取る
parser = argparse.ArgumentParser(description='This script takes five arguments: arg1, arg2, arg3, arg4 and arg5.',
                                 usage='%(prog)s <arg1> <arg2> <arg3> <arg4> <arg5> \
                                 \nexample: python3 %(prog)s N H 0.82 1.24 ../get_some_speceis_existed_poscar_path_list/N_H_existed_poscar_folder_path_list.npy')
parser.add_argument('arg1', help='central_atom_symbol')
parser.add_argument('arg2', help='neighboring_atom_symbol')
parser.add_argument('arg3', help='bond_length_lower_end')
parser.add_argument('arg4', help='bond_length_upper_end')
parser.add_argument('arg5', help='npy_file_path')
args = parser.parse_args()
central_atom_symbol = args.arg1
neighboring_atom_symbol = args.arg2
bond_length_lower_end = args.arg3
bond_length_upper_end = args.arg4
target_npy_p = args.arg5

print(f'central_atom_symbol: {central_atom_symbol}')
print(f'neighboring_atom_symbol: {neighboring_atom_symbol}')
print(f'bond_length_lower_end: {bond_length_lower_end}')
print(f'bond_length_upper_end: {bond_length_upper_end}')
print(f'target_npy_p: {target_npy_p}')
