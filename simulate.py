"""
Note:
    AA: Amino Acid
    Distance unit: same as the one in PDB file (1 * 10^-10 m, i.e. 100pm)
Assumption:
    1. behaviors between two types of AA are same
        type:https://upload.wikimedia.org/wikipedia/commons/a/a9/Amino_Acids.svg
(Possible)Truth:
    1. distance between two continuous atoms is constant
        N-CA: 1.459
        CA-C: 1.525
        C-N: 1.330
    2. angle for certain atom is constant:
        N-CA-C: 111.56
        CA-C-N: 116.53
        C-N-CA: 121.23
"""
import numpy as np
import pandas as pd
import math
import json
from read_PDB import read_pdb, vector_angle_vector, plane_angle_vector, plane_angle_plane, cursor, draw_backbone

# for rough simulating
POSITIVE = ['R', 'K', 'H']
NEGATIVE = ['D', 'E']
HYDROPHOBIC = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
UNCHARGE = ['S', 'T', 'N', 'Q']
SPECIAL = ['C', 'U', 'G', 'P']

SHORT_TO_LONG = {
    'R': 'ARG',
    'H': 'HIS',
    'K': 'LYS',
    'D': 'ASP',
    'E': 'GLU',
    'S': 'SER',
    'T': 'THR',
    'N': 'ASN',
    'Q': 'GLN',
    'C': 'CYS',
    'U': 'SEC',
    'G': 'GLY',
    'P': 'PRO',
    'A': 'ALA',
    'V': 'VAL',
    'I': 'ILE',
    'L': 'LEU',
    'M': 'MET',
    'F': 'PHE',
    'Y': 'TYR',
    'W': 'TRP'
}
# for accurate simulating
N_CA = 1.459
CA_C = 1.525
C_N = 1.33
N_CA_C = 180 - 111.56
CA_C_N = 180 - 116.53
C_N_CA = 180 - 121.23


def structure_data(pdb, chain):
    """
    for one chain, calculate distances and angles
    :param pdb:
    :param chain:
    :return:
    """
    critical_atom_df = read_pdb(pdb, elements=['N', 'CA', 'C'], chain=chain)
    atom_amount = len(critical_atom_df.index)
    assert atom_amount % 3 == 0, 'some aa does not have 3 critical atoms'
    aa_amount = atom_amount / 3
    distances = []  # amount of numbers = atom_amount - 1 (distance between two adjacent atom)
    vv_angles = []  # amount of numbers = atom_amount - 2 (angle for each atom)
    pv_angles = []  # amount of numbers = aa_amount - 1
    pp_angles = []  # amount of numbers = aa_amount - 1
    last_x = critical_atom_df.iloc[0].at['X']
    last_y = critical_atom_df.iloc[0].at['Y']
    last_z = critical_atom_df.iloc[0].at['Z']
    for i in np.arange(0, atom_amount):
        cur_row = critical_atom_df.iloc[i]
        #  calculate distance
        if i == 0:
            distance = np.NaN
        else:
            distance = math.sqrt((cur_row['X'] - last_x) ** 2 + (cur_row['Y'] - last_y)**2 + (cur_row['Z'] - last_z)**2)
        #  calculate vv_angle
        if i == 0 or i == atom_amount - 1:
            vv_angle = np.NaN
        else:
            last_row = critical_atom_df.iloc[i-1]
            next_row = critical_atom_df.iloc[i+1]
            vector1 = (
                cur_row['X'] - last_row['X'],
                cur_row['Y'] - last_row['Y'],
                cur_row['Z'] - last_row['Z']
            )
            vector2 = (
                next_row['X'] - cur_row['X'],
                next_row['Y'] - cur_row['Y'],
                next_row['Z'] - cur_row['Z']
            )
            vv_angle = vector_angle_vector(vector1, vector2)
        #  calculate pv_angle between AAs:
        if cur_row.at['AA sequence'] > 1 and cur_row.at['element long'] == 'N':
            last_aa_v1 = (
                critical_atom_df.iloc[i - 2].at['X'] - critical_atom_df.iloc[i - 3].at['X'],
                critical_atom_df.iloc[i - 2].at['Y'] - critical_atom_df.iloc[i - 3].at['Y'],
                critical_atom_df.iloc[i - 2].at['Z'] - critical_atom_df.iloc[i - 3].at['Z']
            )
            last_aa_v2 = (
                critical_atom_df.iloc[i - 1].at['X'] - critical_atom_df.iloc[i - 2].at['X'],
                critical_atom_df.iloc[i - 1].at['Y'] - critical_atom_df.iloc[i - 2].at['Y'],
                critical_atom_df.iloc[i - 1].at['Z'] - critical_atom_df.iloc[i - 2].at['Z']
            )
            vector3 = (
                critical_atom_df.iloc[i].at['X'] - critical_atom_df.iloc[i - 1].at['X'],
                critical_atom_df.iloc[i].at['Y'] - critical_atom_df.iloc[i - 1].at['Y'],
                critical_atom_df.iloc[i].at['Z'] - critical_atom_df.iloc[i - 1].at['Z']
            )
            pv_angle = plane_angle_vector(last_aa_v1, last_aa_v2, vector3)
        else:
            pv_angle = np.NaN
        #  calculate pp_angle between AAs:
        if cur_row.at['AA sequence'] > 1 and cur_row.at['element long'] == 'C':
            last_aa_v1 = (
                critical_atom_df.iloc[i - 4].at['X'] - critical_atom_df.iloc[i - 5].at['X'],
                critical_atom_df.iloc[i - 4].at['Y'] - critical_atom_df.iloc[i - 5].at['Y'],
                critical_atom_df.iloc[i - 4].at['Z'] - critical_atom_df.iloc[i - 5].at['Z']
            )
            last_aa_v2 = (
                critical_atom_df.iloc[i - 3].at['X'] - critical_atom_df.iloc[i - 4].at['X'],
                critical_atom_df.iloc[i - 3].at['Y'] - critical_atom_df.iloc[i - 4].at['Y'],
                critical_atom_df.iloc[i - 3].at['Z'] - critical_atom_df.iloc[i - 4].at['Z']
            )
            cur_aa_v1 = (
                critical_atom_df.iloc[i - 1].at['X'] - critical_atom_df.iloc[i - 2].at['X'],
                critical_atom_df.iloc[i - 1].at['Y'] - critical_atom_df.iloc[i - 2].at['Y'],
                critical_atom_df.iloc[i - 1].at['Z'] - critical_atom_df.iloc[i - 2].at['Z']
            )
            cur_aa_v2 = (
                critical_atom_df.iloc[i].at['X'] - critical_atom_df.iloc[i - 1].at['X'],
                critical_atom_df.iloc[i].at['Y'] - critical_atom_df.iloc[i - 1].at['Y'],
                critical_atom_df.iloc[i].at['Z'] - critical_atom_df.iloc[i - 1].at['Z']
            )
            pp_angle = plane_angle_plane(last_aa_v1, last_aa_v2, cur_aa_v1, cur_aa_v2)
        else:
            pp_angle = np.NaN

        distances.append(distance)
        vv_angles.append(vv_angle)
        pv_angles.append(pv_angle)
        pp_angles.append(pp_angle)
        last_x = cur_row.at['X']
        last_y = cur_row.at['Y']
        last_z = cur_row.at['Z']
    critical_atom_df['distance'] = distances
    critical_atom_df['atom angle'] = vv_angles
    critical_atom_df['CN bond angle'] = pv_angles
    critical_atom_df['AA angle'] = pp_angles
    return critical_atom_df


def read_db():
    cursor.execute(
        """
        select first_aa_name, second_aa_name, avg(CN_a) as CN_a, avg(CN_b) as CN_b, 
        avg(CN_c) as CN_c,  avg(NCA_a) as NCA_a, avg(NCA_b) as NCA_b, avg(NCA_c) as NCA_c, 
        avg(CAC_a) as CAC_a, avg(CAC_b) as CAC_b, avg(CAC_c) as CAC_c from two_adjacent_AAs
        group by first_aa_name, second_aa_name
        """
    )
    results = cursor.fetchall()
    with open('parameter.txt', 'w') as outfile:
        json.dump(outfile, results)


def calculate_parameters():
    cursor.execute(
        """
        SELECT DISTINCT first_aa_name, second_aa_name 
        FROM two_adjacent_AAs
        ORDER BY first_aa_name, second_aa_name
        """
    )
    results = cursor.fetchall()
    data = {}
    for combination in results:
        first = combination['first_aa_name']
        last = combination['second_aa_name']
        cursor.execute(
            """
            SELECT avg(CN_a) AS CN_a,
                   avg(CN_b) AS CN_b,
                   avg(CN_c) AS CN_c,
                   avg(NCA_a) AS NCA_a,
                   avg(NCA_b) AS NCA_b,
                   avg(NCA_c) AS NCA_c,
                   avg(CAC_a) AS CAC_a,
                   avg(CAC_b) AS CAC_b,
                   avg(CAC_c) AS CAC_c
            FROM two_adjacent_AAs
            WHERE first_aa_name = '{first}'
              AND second_aa_name = '{last}'
            """.format(
                first=first,
                last=last
            )
        )
        result = cursor.fetchone()
        data[first + last] = {
            'CN_a': result['CN_a'],
            'CN_b': result['CN_b'],
            'CN_c': result['CN_c'],
            'NCA_a': result['NCA_a'],
            'NCA_b': result['NCA_b'],
            'NCA_c': result['NCA_c'],
            'CAC_a': result['CAC_a'],
            'CAC_b': result['CAC_b'],
            'CAC_c': result['CAC_c'],
        }
    with open('parameter.json', 'w') as outfile:
        json.dump(data, outfile)


def generate_coordinates(aas):
    df = pd.DataFrame(
        [[1, 'N', 1, 0, 0, 0],
        [2, 'CA', 1, N_CA, 0, 0],
        [3, 'C', 1, N_CA + math.cos(N_CA_C) * CA_C, math.sin(N_CA_C), 0]],
        columns=['atom sequence', 'element long', 'AA sequence', 'X', 'Y', 'Z']
    )  # default start is N at (0, 0, 0)
    # df = pd.DataFrame(
    #     [[1, 'N', 1, -6.966, -25.812, 7.869],
    #      [2, 'CA', 1, -8.004, -25.537, 6.895],
    #      [3, 'C', 1, -7.837, -24.189, 6.222]],
    #     columns=['atom sequence', 'element long', 'AA sequence', 'X', 'Y', 'Z']
    # )
    # next is CA on x-axis
    # the 3rd atom is C on the x-y plane
    aa_amount = len(aas)
    for i in np.arange(1, aa_amount):
        last_aa = aas[i - 1]
        cur_aa = aas[i]
        # cursor.execute(
        #     """
        #     SELECT avg(CN_a) AS CN_a,
        #            avg(CN_b) AS CN_b,
        #            avg(CN_c) AS CN_c,
        #            avg(NCA_a) AS NCA_a,
        #            avg(NCA_b) AS NCA_b,
        #            avg(NCA_c) AS NCA_c,
        #            avg(CAC_a) AS CAC_a,
        #            avg(CAC_b) AS CAC_b,
        #            avg(CAC_c) AS CAC_c
        #     FROM two_adjacent_AAs
        #     WHERE first_aa_name = '{first}'
        #       AND second_aa_name = '{last}'
        #     """.format(
        #         first=SHORT_TO_LONG[last_aa],
        #         last=SHORT_TO_LONG[cur_aa]
        #     )
        # )
        json_parameter = json.load(open('parameter.json', 'r'))
        result = json_parameter[SHORT_TO_LONG[last_aa] + SHORT_TO_LONG[cur_aa]]
        last_n = df.iloc[-3]
        last_ca = df.iloc[-2]
        last_c = df.iloc[-1]
        last_v1 = np.array([
            last_ca.at['X'] - last_n.at['X'],
            last_ca.at['Y'] - last_n.at['Y'],
            last_ca.at['Z'] - last_n.at['Z']
        ])
        last_v2 = np.array([
            last_c.at['X'] - last_ca.at['X'],
            last_c.at['Y'] - last_ca.at['Y'],
            last_c.at['Z'] - last_ca.at['Z']
        ])
        new_cn = result['CN_a'] * last_v1 + result['CN_b'] * last_v2 + result['CN_c'] * np.cross(last_v1, last_v2)
        new_nca = result['NCA_a'] * last_v1 + result['NCA_b'] * last_v2 + result['NCA_c'] * np.cross(last_v1, last_v2)
        new_cac = result['CAC_a'] * last_v1 + result['CAC_b'] * last_v2 + result['CAC_c'] * np.cross(last_v1, last_v2)
        new_n = new_cn / np.linalg.norm(new_cn) * C_N + np.array([last_c.at['X'], last_c.at['Y'], last_c.at['Z']])
        new_ca = new_n + new_nca / np.linalg.norm(new_nca) * N_CA
        new_c = new_ca + new_cac / np.linalg.norm(new_cac) * CA_C
        # new_n = new_cn + np.array([last_c.at['X'], last_c.at['Y'], last_c.at['Z']])
        # new_ca = new_n + new_nca
        # new_c = new_ca + new_cac
        new_rows = pd.DataFrame(
            [[i * 3 + 1, 'N', i + 1, new_n[0], new_n[1], new_n[2]],
            [i * 3 + 2, 'CA', i + 1, new_ca[0], new_ca[1], new_ca[2]],
            [i * 3 + 3, 'C', i + 1, new_c[0], new_c[1], new_c[2]]],
            columns=['atom sequence', 'element long', 'AA sequence', 'X', 'Y', 'Z']
        )
        df = pd.concat([df, new_rows])
    return df


def rough_simulating(aas, x=0, y=0, z=0):
    df = pd.DataFrame(
        [1, x, y, z],
        columns=['AA sequence', 'X', 'Y', 'Z']
    )  # default start is N at (x, y, z)
    length = len(aas)
    for i in np.arange(0, length - 1):
        last_x = df.iloc[-1].at['X']
        last_y = df.iloc[-1].at['Y']
        last_z = df.iloc[-1].at['Z']
        if aas[i] in POSITIVE and aas[i+1] in POSITIVE:
            delta_x = -2.197  # 1AO6 10 - 9
            delta_y = -2.24
            delta_z = -2.066
        elif aas[i] in POSITIVE and aas[i+1] in NEGATIVE:
            delta_x = -3.344  # 1AO6 13 - 12
            delta_y = -1.251
            delta_z = -1.391
        elif aas[i] in NEGATIVE and aas[i + 1] in POSITIVE:
            delta_x = 3.344  # 1AO6 13 - 12
            delta_y = 1.251
            delta_z = 1.391
        elif aas[i] in NEGATIVE and aas[i+1] in NEGATIVE:
            delta_x = -3.344  # 1AO6 17 - 16
            delta_y = -1.251
            delta_z = -1.391

        new_row = {
            'AA sequence': i + 2,
            'X': last_x + delta_x
        }


def simulate(filepath, element='CA', save_as=None):
    fasta = open(filepath, 'r')
    fasta.readline()
    for line in fasta:
        aas = list(line)[:-1]
        break
    atoms_df = generate_coordinates(aas)
    atoms_df = atoms_df[atoms_df['element long'] == element]
    draw_backbone(atoms_df, save_as)


# calculate_parameters()
simulate('E:\Research\protein_turtle_ura-master\samples\\KAMP\\5KI0_A.fasta.txt',
         save_as='E:\Research\protein_turtle_ura-master\GIFs\\5KI0_A_Jul4.gif')