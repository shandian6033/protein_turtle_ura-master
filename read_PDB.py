"""
TODO: goal: draw a draft
"""
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import animation
import mpl_toolkits.mplot3d as plt3d
import os
import sqlite3

CWD = os.getcwd()
PDB_LOCATION = CWD + '\PDBs\\'
COLUMNS = ['item', 'atom_sequence', 'element_long', 'amino_acid', 'chain', 'aa_sequence', 'X', 'Y', 'Z', 'occupancy',
           'temperature factor', 'element']
NUMERIC_COLUMNS = ['atom_sequence', 'aa_sequence', 'X', 'Y', 'Z', 'occupancy', 'temperature factor']
CWD = os.getcwd()
fig = plt.figure()
ax = plt3d.Axes3D(fig)
VALID_CHAINS = [
    ('3I3Z', 'A'),
    ('3I3Z', 'B'),
    ('1AO6', 'A'),
    ('5SV4', 'A'),
    ('1FUW', 'A'),
    ('1IV7', 'A'),
    ('3VVJ', 'A'),
    ('1A4V', 'A'),
    ('2A3D', 'A'),
    ('1AKP', 'A'),
    ('1AP2', 'A'),
    ('1AP2', 'B'),
    ('5KI0', 'A'),
    ('5M4S', 'A'),
    ('1B9P', 'A')
]


def dict_factory(cur, row):
    d = {}
    for idx, col in enumerate(cur.description):
        d[col[0]] = row[idx]
    return d


connect = sqlite3.connect(CWD + '\protein.db')
connect.row_factory = dict_factory
cursor = connect.cursor()


def read_pdb(filename):
    """

    :param filename: (str)
    :return: (list of list of str)

    """
    file = open(PDB_LOCATION + filename)
    for line in file:
        if line[:4] == 'ATOM':
            words = [
                line[:8].strip(' '),
                int(line[8:12]),
                line[12:17].strip(' '),
                line[17:21].strip(' '),
                line[21:23].strip(' '),
                int(line[23:31]),
                float(line[31:38]),
                float(line[38:46]),
                float(line[46:54]),
                float(line[54:60]),
                float(line[60:66]),
                line[77]
            ]
            cursor.execute(
                """
                DELETE FROM protein_atoms
                WHERE
                    protein_id = '{protein}'
                    AND chain = '{chain}'
                    AND atom_sequence = {atom_sequence}
                """.format(
                    protein=filename[:4],
                    chain=words[4],
                    atom_sequence=words[1]
                )
            )
            cursor.execute(
                """
                INSERT INTO protein_atoms
                    (protein_id,
                    chain,
                    amino_acid,
                    aa_sequence,
                    element_long,
                    atom_sequence,
                    X,
                    Y,
                    Z)
                VALUES
                    ('{protein}',
                    '{chain}',
                    '{amino_acid}',
                    {aa_sequence},
                    '{element_long}',
                    {atom_sequence},
                    {X},
                    {Y},
                    {Z})
                """.format(
                    protein=filename[:4],
                    chain=words[4],
                    amino_acid=words[3],
                    aa_sequence=words[5],
                    element_long=words[2],
                    atom_sequence=words[1],
                    X=words[6],
                    Y=words[7],
                    Z=words[8]
                )
            )
    connect.commit()
    # df = pd.DataFrame(atoms, columns=COLUMNS)
    # df = df.sort_values(['atom_sequence']).filter(
    #     items=['chain', 'amino_acid', 'aa_sequence', 'element_long', 'atom_sequence', 'X', 'Y', 'Z']
    # )
    # df['protein_id'] = filename[:4]
    #
    # df.to_sql('protein_atoms', con=connect, if_exists='append', index=False)
    # return df


def rotate(angle):
    ax.view_init(azim=angle)


def draw_backbone(df, save_as=None):
    fig.set_size_inches(10, 10)
    # ax = fig.add_subplot(111, projection='3d', aspect='equal')
    # ax.view_init(azim=120)

    is_first = True
    for index, row in df.iterrows():
        if is_first:
            is_first = False
            # plot the nodes
            x = row['X']
            y = row['Y']
            z = row['Z']
            ax.scatter(x, y, z, color='red', marker='s')
            continue
        else:
            # plot the line
            xs = x, row['X']
            ys = y, row['Y']
            zs = z, row['Z']
            line = plt3d.art3d.Line3D(xs, ys, zs)
            ax.add_line(line)
            # plot the nodes
            x = row['X']
            y = row['Y']
            z = row['Z']
            ax.scatter(x, y, z, color='black', marker='s')
    # plt.show()
    if save_as:
        rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 360, 3), interval=150)
        rot_animation.save(save_as, dpi=80, writer='imagemagick')


def calculate_distance(df):
    first = True
    df['distance'] = np.NaN
    for index, row in df.iterrows():
        if not first:
            distance = math.sqrt((row['X'] - last_x) ** 2 + (row['Y'] - last_y)**2 + (row['Z'] - last_z)**2)
            # df.loc[:, ('distance', index)] = distance
            df[index, 'distance'] = distance
        else:
            first = False
        last_x = row['X']
        last_y = row['Y']
        last_z = row['Z']
    return df


def vector_angle_vector(v1, v2):
    """

    :param v1: tuple
    :param v2: tuple
    :return:
    angle in radian
    0 <= angle <= 180
    """
    cosang = np.dot(v1, v2)
    sinang = np.linalg.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang) / np.pi * 180


def plane_angle_vector(plane_v1, plane_v2, v3):
    """
    0 <= angle <= 90
    :param plane_v1:
    :param plane_v2:
    :param v3:
    :return:
    """
    perpendicular_vector = np.cross(plane_v1, plane_v2)
    return abs(90 - vector_angle_vector(perpendicular_vector, v3))


def plane_angle_plane(p1_v1, p1_v2, p2_v1, p2_v2):
    """
    0 <= angle <= 90
    :param p1_v1:
    :param p1_v2:
    :param p2_v1:
    :param p2_v2:
    :return:
    """
    p1_pvector = np.cross(p1_v1, p1_v2)  # perpendicular vector of plane 1
    p2_pvector = np.cross(p2_v1, p2_v2)  # perpendicular vector of plane 2
    pvectors_angle = vector_angle_vector(p1_pvector, p2_pvector)
    if pvectors_angle >= 90:
        return 180 - pvectors_angle
    else:
        return pvectors_angle


def calculate_angle(df):
    length = len(df.index)
    angles = [np.NaN]
    for i in np.arange(1, length-1):
        last_row = df[i-1]
        cur_row = df[i]
        next_row = df[i+1]
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
        angles.append(vector_angle_vector(vector1, vector2))
    angles.append(np.NaN)
    df['radian'] = angles
    return df


def is_constant(df, column, reference_atom):
    """
    based on the distance bewteen two atoms, see if they are constant.
    Args:
        df:
        column:
            which column to check. (one of following)
            'distance': check if the distance between 2 atoms are constant
            'radian': check if the angle at an atom is constant
        reference_atom:
            'CA': distance between CA and N OR the angle of N-CA-C
            'C': CA & C OR the angle of CA-C-N
            'N': C & N(for next atom) OR the angle of C-N-CA
    :return:
    """
    df = df[df['element long'] == reference_atom]
    to_check = df[column]
    return {
        'checking': column,
        'atom': reference_atom,
        'mean': to_check.mean(),
        'standard deviation': to_check.std()
    }


def solve_functions(v1, v2, v3):
    """
    a * v1 + b * v2 + c * v1 dot v2 = v3
    :param v1:
    :param v2:
    :param v3:
    :return:
    """
    v12 = np.cross(v1, v2)
    d = np.array([
        [v1[0], v2[0], v12[0]],
        [v1[1], v2[1], v12[1]],
        [v1[2], v2[2], v12[2]]
    ])
    v3 = np.array(v3)
    return np.matmul(np.linalg.inv(d), v3)


def structure_data(protein, chain):
    """
    for one chain, calculate angles between AAs and store to db
    :param pdb:
    :param chain:
    :return:
    """
    cursor.execute(
        """
        SELECT 
            protein_id, chain, amino_acid, aa_sequence, element_long, atom_sequence, X, Y, Z
        FROM
            protein_atoms
        WHERE
            protein_id = '{0}'
            AND chain = '{1}'
            AND element_long IN ('CA', 'C', 'N', 'CA A', 'C  A', 'N  A')
        ORDER BY
            atom_sequence
        """.format(protein, chain)
    )
    results = cursor.fetchall()
    atom_amount = len(results)
    assert atom_amount % 3 == 0, 'some aa does not have 3 critical atoms'
    aa_amount = int(atom_amount / 3)
    for i in list(range(1, aa_amount)):
        #  calculate pv_angle between AAs:
        last_N_CA = np.array([
            results[i * 3 - 2]['X'] - results[i * 3 - 3]['X'],
            results[i * 3 - 2]['Y'] - results[i * 3 - 3]['Y'],
            results[i * 3 - 2]['Z'] - results[i * 3 - 3]['Z']
        ])
        last_CA_C = np.array([
            results[i * 3 - 1]['X'] - results[i * 3 - 2]['X'],
            results[i * 3 - 1]['Y'] - results[i * 3 - 2]['Y'],
            results[i * 3 - 1]['Z'] - results[i * 3 - 2]['Z']
        ])
        CN_bond = np.array([
            results[i * 3]['X'] - results[i * 3 - 1]['X'],
            results[i * 3]['Y'] - results[i * 3 - 1]['Y'],
            results[i * 3]['Z'] - results[i * 3 - 1]['Z']
        ])
        pv_angle = plane_angle_vector(last_N_CA, last_CA_C, CN_bond)
        #  calculate pp_angle between AAs:
        cur_N_CA = np.array([
            results[i * 3 + 1]['X'] - results[i * 3]['X'],
            results[i * 3 + 1]['Y'] - results[i * 3]['Y'],
            results[i * 3 + 1]['Z'] - results[i * 3]['Z']
        ])
        cur_CA_C = np.array([
            results[i * 3 + 2]['X'] - results[i * 3 + 1]['X'],
            results[i * 3 + 2]['Y'] - results[i * 3 + 1]['Y'],
            results[i * 3 + 2]['Z'] - results[i * 3 + 1]['Z']
        ])
        pp_angle = plane_angle_plane(last_N_CA, last_CA_C, cur_N_CA, cur_CA_C)
        #  calculate C-N angle between AAs:
        c_angle = vector_angle_vector(last_CA_C, CN_bond)
        n_angle = vector_angle_vector(CN_bond, cur_N_CA)

        #  set vector to unit vector
        last_N_CA = last_N_CA / np.linalg.norm(last_N_CA)
        last_CA_C = last_CA_C / np.linalg.norm(last_CA_C)
        CN_bond = CN_bond / np.linalg.norm(CN_bond)
        cur_N_CA = cur_N_CA / np.linalg.norm(cur_N_CA)
        cur_CA_C = cur_CA_C / np.linalg.norm(cur_CA_C)
        #
        draw_cn = solve_functions(last_N_CA, last_CA_C, CN_bond)
        draw_nca = solve_functions(last_N_CA, last_CA_C, cur_N_CA)
        draw_cac = solve_functions(last_N_CA, last_CA_C, cur_CA_C)

        last_aa = results[i * 3 -1]['amino_acid']
        cur_aa = results[i * 3]['amino_acid']
        cursor.execute(
            """
            DELETE FROM two_adjacent_AAs
            WHERE 
                first_aa_sequence = '{first_aa_sequence}'
                AND protein_id = '{protein}'
                AND chain = '{chain}'
            """.format(
                first_aa_sequence=results[i * 3 - 1]['aa_sequence'],
                protein=protein,
                chain=chain
            )
        )
        cursor.execute(
            """
            INSERT INTO two_adjacent_AAs
                        (first_aa_name,
                         second_aa_name,
                         combined_aas,
                         ca_cn_angle,
                         cn_ca_angle,
                         p_cn_angle,
                         pp_angle,
                         protein_id,
                         chain,
                         first_aa_sequence,
                         CN_a,
                         CN_b,
                         CN_c,
                         NCA_a,
                         NCA_b,
                         NCA_c,
                         CAC_a,
                         CAC_b,
                         CAC_c)
            VALUES      ('{last_aa}',
                         '{cur_aa}',
                         '{combine}',
                         {CA_CN_angle},
                         {CN_CA_angle},
                         {p_CN_angle},
                         {pp_angle},
                         '{protein_id}',
                         '{chain}',
                         {first_aa_sequence},
                         {CN_a},
                         {CN_b},
                         {CN_c},
                         {NCA_a},
                         {NCA_b},
                         {NCA_c},
                         {CAC_a},
                         {CAC_b},
                         {CAC_c}) 
            """.format(
                last_aa=last_aa,
                cur_aa=cur_aa,
                combine=last_aa + '_' + cur_aa if last_aa < cur_aa else cur_aa + '_' + last_aa,
                CA_CN_angle=c_angle,
                CN_CA_angle=n_angle,
                p_CN_angle=pv_angle,
                pp_angle=pp_angle,
                protein_id=protein,
                chain=chain,
                first_aa_sequence=results[i * 3 - 1]['aa_sequence'],
                CN_a=draw_cn[0],
                CN_b=draw_cn[1],
                CN_c=draw_cn[2],
                NCA_a=draw_nca[0],
                NCA_b=draw_nca[1],
                NCA_c=draw_nca[2],
                CAC_a=draw_cac[0],
                CAC_b=draw_cac[1],
                CAC_c=draw_cac[2]
            )
        )
        connect.commit()


def draw_chain(protein_id, chain, element='CA'):
    cursor.execute(
        """
        SELECT X, Y, Z
        FROM protein_atoms
        WHERE protein_id = '{protein}'
            AND chain = '{chain}'
            AND element_long = '{element}'
        ORDER BY atom_sequence
        """.format(
            protein=protein_id,
            chain=chain,
            element=element
        )
    )
    results = cursor.fetchall()
    df = pd.DataFrame(results)
    draw_backbone(df, '{0}\GIFs\\real_{1}_{2}.gif'.format(CWD, protein_id, chain))


def calculate_valid_parameters():
    for pair in VALID_CHAINS:
        structure_data(pair[0], pair[1])


def read_all_pdbs():
    files = [f for f in os.listdir(PDB_LOCATION) if os.path.isfile(os.path.join(PDB_LOCATION, f))]
    for file in files:
        if file[-3:] == 'pdb':
            read_pdb(file)


if __name__ == '__main__':
    # read_all_pdbs()
    # calculate_valid_parameters()
    # structure_data('1IV7', 'A')
    # structure_data('3VVJ', 'A')
    draw_chain('5KI0', 'A')
