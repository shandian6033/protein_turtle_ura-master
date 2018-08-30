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


def calculate_parameters():
    """
    retrieve parameters (i.e. CN_a,CN_b ...) for simulation and store to parameter.json.
    This function is to group parameters by adjacent AAs ignoring order. (i.e. GLU LYS is same as LYS GLU)
    """
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
    """
    With parameters, generate simulated protein structure.
    :param aas: (str): the aa sequence of a protein chain
    :return: (pandas.DataFrame):
        a data frame containing simulated protein structure
    """
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


def simulate(filepath, element='CA', save_as=None):
    """
    given a protein chain, simulate its 3d structure
    :param filepath: (str):
        path of a fasta file containing aa sequence.
        The aa sequence should not be longer than one line.
    :param element: (str):
        one of following: 'CA', 'C', 'N'
        which atom in aa is used to draw
    :param save_as: (str):
        refer to draw_backbone()
    :return: null
    """
    fasta = open(filepath, 'r')
    fasta.readline()
    for line in fasta:
        aas = list(line)[:-1]
        break
    atoms_df = generate_coordinates(aas)
    atoms_df = atoms_df[atoms_df['element long'] == element]
    draw_backbone(atoms_df, save_as)


simulate('E:\Research\protein_turtle_ura-master\samples\\KAMP\\5KI0_A.fasta.txt',
         save_as='E:\Research\protein_turtle_ura-master\GIFs\\5KI0_A_Jul4.gif')