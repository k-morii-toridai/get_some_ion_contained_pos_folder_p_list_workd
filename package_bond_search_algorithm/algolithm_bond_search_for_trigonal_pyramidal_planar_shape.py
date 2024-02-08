import numpy as np


def filter_2(df_nnlist, central_atom_symbol, neighboring_atom_symbol, bond_length_lower_end, bond_length_upper_end):
    """
    2．POSCAR.nnlistにおいて，原子'C'から0-2．のNH結合距離以内に，原子'O'を4つ含む，中心原子'C'が存在するかどうか判定．
    → 存在する場合、True値，{中心原子'C'の'central_atom_id': そのneighborsの'central_atom_id'}の辞書の2つを返す．
    → 存在しない場合，False値，空の辞書を返す．

    Usage:
    ------

    bool_2, dict_2 = filter_2(df_nnlist=df_nnlist, bond_length_lower_end=bond_length_lower_end, bond_length_upper_end=bond_length_upper_end)

    Parameters:
    -----------
    df_nnlist: pd.DataFrame

    Returns:
    --------
    bool_2: bool
    dict_2: dict
    """
    # df_nnlist_group_dict = df_nnlist[df_nnlist['central_atom_symbol'] == 'C'].groupby('central_atom_id').groups
    df_nnlist_group_dict = df_nnlist[df_nnlist['central_atom_symbol'] == central_atom_symbol].groupby('central_atom_id').groups
    df_nnlist_central_atom_ids = np.array(list(df_nnlist_group_dict.keys()))
    bool_list = []
    for key in df_nnlist_central_atom_ids:
        # bool_sort_dist_list = bond_length_lower_end < df_nnlist.iloc[df_nnlist_group_dict[key]].sort_values('rel_distance')['rel_distance'] < bond_length_upper_end
        bond_length_lower_end_c_eq = bond_length_lower_end < df_nnlist.iloc[df_nnlist_group_dict[key]].sort_values('rel_distance')['rel_distance']
        bond_length_upper_end_c_eq = df_nnlist.iloc[df_nnlist_group_dict[key]].sort_values('rel_distance')['rel_distance'] < bond_length_upper_end
        bool_sort_dist_list = bond_length_lower_end_c_eq & bond_length_upper_end_c_eq
        # bool_list.append(df_nnlist.iloc[df_nnlist_group_dict[key]].sort_values('rel_distance')[bool_sort_dist_list]['neighboring_atom_symbol'].tolist().count('O') == 3)
        bool_list.append(df_nnlist.iloc[df_nnlist_group_dict[key]].sort_values('rel_distance')[bool_sort_dist_list]['neighboring_atom_symbol'].tolist().count(neighboring_atom_symbol) == 3)
    df_nnlist_central_atom_ids_fillterd = df_nnlist_central_atom_ids[bool_list]
    bool_filter_2 = len(df_nnlist_central_atom_ids_fillterd) >= 1
    filtered_df_nnlist_group_dict = {key: df_nnlist_group_dict[key] for key in df_nnlist_group_dict.keys() if key in df_nnlist_central_atom_ids_fillterd}

    return bool_filter_2, filtered_df_nnlist_group_dict


def filter_3(df_nnlist, dict_2, neighboring_atom_symbol):
    """
    3．2．の中心原子'C'に対して，1番近い原子が'O'であり，かつ2番目・3番目・4番目に近い原子も'O'である，中心原子'C'が存在するかどうか判定．
    → 存在する場合，True値，{中心原子Cの'central_atom_id': そのneighborsの'central_atom_id'}の辞書の2つを返す．
    → 存在しない場合，False値，空の辞書を返す．

    Usage:
    ------
    bool_3, dict_3 = filter_3(df_nnlist=df_nnlist, dict_2=dict_2, neighboring_atom_symbol=neighboring_atom_symbol)

    Parameters:
    -----------
    df_nnlist: pd.DataFrame
    dict_2: dict
    neighboring_atom_symbol: str

    Returns:
    --------
    bool_3: bool
    dict_3: dict
    """
    bool_list_3 = []
    for k, v in dict_2.items():
        # bool_list_3.append(set(df_nnlist.iloc[dict_2[k]].sort_values(by='rel_distance')['neighboring_atom_symbol'].tolist()[1:5]) == {'O'})
        bool_list_3.append(set(df_nnlist.iloc[dict_2[k]].sort_values(by='rel_distance')['neighboring_atom_symbol'].tolist()[1:4]) == {neighboring_atom_symbol})
    df_nnlist_central_atom_ids_fillterd_3 = np.array(list(dict_2.keys()))[bool_list_3]
    filtered_3_df_nnlist_group_dict = {key: dict_2[key] for key in dict_2.keys() if key in df_nnlist_central_atom_ids_fillterd_3}
    bool_filter_3 = len(df_nnlist_central_atom_ids_fillterd_3) >= 1

    return bool_filter_3, filtered_3_df_nnlist_group_dict


def filter_4(df_nnlist, dict_3):
    """
    4．2．の中心原子'C'に対して5番目に近い原子が存在しないかを判定．
    → 存在しない場合，True値，{中心原子'C'の'central_atom_id': そのneighborsの'central_atom_id'}の辞書の2つを返す．
    → 存在する場合，False値，空の辞書の2つを返す．

    Usage:
    ------
    bool_4, dict_4 = filter_4(df_nnlist=df_nnlist, dict_3=dict_3)

    Parameters:
    -----------
    df_nnlist: pd.DataFrame
    dict_3: dict

    Returns:
    --------
    bool_4: bool
    dict_4: dict
    """
    bool_list_4 = []
    for k, v in dict_3.items():
        bool_list_4.append(len(df_nnlist.iloc[dict_3[k]].sort_values(by='rel_distance')['neighboring_atom_symbol'].tolist()) == 4)
    df_nnlist_central_atom_ids_fillterd_4 = np.array(list(dict_3.keys()))[bool_list_4]
    filtered_4_df_nnlist_group_dict = {key: dict_3[key] for key in dict_3.keys() if key in df_nnlist_central_atom_ids_fillterd_4}
    bool_filter_4 = len(df_nnlist_central_atom_ids_fillterd_4) >= 1

    return bool_filter_4, filtered_4_df_nnlist_group_dict


def filter_5(df_nnlist, dict_3):
    """
    5．2．の中心原子'C'に対して5番目に近い原子が，'C'に4番目に近い原子'O'と'C'のNH距離より大きい'C'が存在するどうかを判定．
    → 存在する場合，True値，{中心原子'C'の'central_atom_id': そのneighborsの'central_atom_id'}の辞書の2つを返す．
    → 存在しない場合，False値，空の辞書を返す．

    Usage:
    ------
    bool_5, dict_5 = filter_5(df_nnlist=df_nnlist, dict_3=dict_3)

    Parameters:
    -----------
    df_nnlist: pd.DataFrame
    dict_3: dict

    Returns:
    --------
    bool_5: bool
    dict_5: dict
    """
    bool_list_5 = []
    for k, v in dict_3.items():
        third_bond_dist = df_nnlist.iloc[dict_3[k]].sort_values(by='rel_distance')['rel_distance'].tolist()[3]
        fourth_dist = df_nnlist.iloc[dict_3[k]].sort_values(by='rel_distance')['rel_distance'].tolist()[4]
        bool_list_5.append(third_bond_dist < fourth_dist)
    df_nnlist_central_atom_ids_fillterd_5 = np.array(list(dict_3.keys()))[bool_list_5]
    filtered_5_df_nnlist_group_dict = {key: dict_3[key] for key in dict_3.keys() if key in df_nnlist_central_atom_ids_fillterd_5}
    bool_filter_5 = len(df_nnlist_central_atom_ids_fillterd_5) >= 1

    return bool_filter_5, filtered_5_df_nnlist_group_dict


def filter_6(df_nnlist, dict_3, central_atom_symbol):
    """
    6．3．の4つの原子'O'全てに対して，3．の中心の原子'C'との距離以内に，中心原子'C'以外の別の原子が存在しないかどうかを判定．
    → 存在しない場合，True値，中心原子'C'の'central_atom_id'のndarrayの2つを返す．
    → 存在する場合，False値，空のndarrayを返す．

    Usage:
    ------
    bool_6, N_ids = filter_6(df_nnlist=df_nnlist, dict_3=dict_3, central_atom_symbol=central_atom_symbol)

    Parameters:
    -----------
    df_nnlist: pd.DataFrame
    dict_3: dict
    central_atom_symbol: str

    Returns:
    --------
    bool_6: bool
    N_ids: ndarray
    """
    bool_list_6 = []
    for k, v in dict_3.items():
        # N周りのH4つのindex
        indices = df_nnlist.iloc[dict_3[k]].sort_values(by='rel_distance').index[1:4]
        O_ids = df_nnlist.iloc[indices].apply(lambda row: row['neighboring_atom_id'], axis=1).tolist()
        bool_list_temp = []
        for O_id in O_ids:
            bool_temp = df_nnlist[df_nnlist['central_atom_id'] == O_id].sort_values('rel_distance')['neighboring_atom_symbol'].tolist()[1] == central_atom_symbol
            bool_list_temp.append(bool_temp)
        if set(bool_list_temp) == {True}:
            bool_list_6.append(True)
        else:
            bool_list_6.append(False)
    O_ids = np.array(list(dict_3.keys()))[bool_list_6]
    bool_filter_6 = bool_list_6.count(True) >= 1

    return bool_filter_6, O_ids


def concat_filter(df_nnlist, central_atom_symbol, neighboring_atom_symbol, bond_length_lower_end, bond_length_upper_end):
    """
    filter_2()~filter_6()の関数を用いて，POSCAR.nnlistを用いて，POSCARファイルにアンモニウムイオンを含むかどうかの判定algolithmを作成．
    → True値が返された場合，アンモニウムイオンを含む．
    → False値が返され場合，アンモニウムイオンを含まない．

    Usage:
    ------
    concat_filter(df_nnlist=df_nnlist,
        central_atom_symbol=central_atom_symbol,
        neighboring_atom_symbol=neighboring_atom_symbol,
        bond_length_lower_end=bond_length_lower_end,
        bond_length_upper_end=bond_length_upper_end)

    Parameters:
    -----------
    df_nnlist: pd.DataFrame
    central_atom_symbol: str: ex.) 'C'
    neighboring_atom_symbol: str: ex.) 'O'
    bond_length_lower_end: float or int: ex.) 0.82 ← Unit is Å
    bond_length_upper_end: float or int: ex.) 1.24 ← Unit is Å

    Returns:
    --------
    bool: True or False
    """
    bool_2, dict_2 = filter_2(df_nnlist, central_atom_symbol, neighboring_atom_symbol, bond_length_lower_end, bond_length_upper_end)
    if bool_2:
        bool_3, dict_3 = filter_3(df_nnlist=df_nnlist, dict_2=dict_2, neighboring_atom_symbol=neighboring_atom_symbol)
        if bool_3:
            bool_4, dict_4 = filter_4(df_nnlist=df_nnlist, dict_3=dict_3)
            if bool_4:
                return True
            else:
                bool_5, dict_5 = filter_5(df_nnlist=df_nnlist, dict_3=dict_3)
                if bool_5:
                    bool_6, N_ids = filter_6(df_nnlist=df_nnlist, dict_3=dict_3, central_atom_symbol=central_atom_symbol)
                    if bool_6:
                        return True
                    else:
                        return False
                else:
                    return False
        else:
            return False
    else:
        return False


if __name__ == '__main__':
    from package_file_conversion.nnlist2df import nnlist2df
    nnlist_path = 'sample_test_files/1000033/nnlist_5/POSCAR.nnlist'
    df_nnlist = nnlist2df(nnlist_path=nnlist_path)
    print(concat_filter(df_nnlist, central_atom_symbol='C', neighboring_atom_symbol='O', bond_length_lower_end=0.99, bond_length_upper_end=1.66))
