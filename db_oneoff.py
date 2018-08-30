"""
This file is used to create a database in the beginning. This file should be run only once.
"""
import sqlite3

DB_FILE = 'E:\Research\protein_turtle_ura-master\protein.db'

connect = sqlite3.connect(DB_FILE)
cursor = connect.cursor()


def create_protein_structure():
    cursor.execute(
        """
        CREATE TABLE `protein_atoms` (
            `protein_id`	TEXT,
            `chain`	TEXT,
            `amino_acid`	TEXT,
            `aa_sequence`	INTEGER,
            `element_long`	TEXT,
            `atom_sequence`	INTEGER,
            `X`	REAL,
            `Y`	REAL,
            `Z`	REAL,
            PRIMARY KEY(`protein_id`,`chain`,`atom_sequence`)
        );
        """
    )
    connect.commit()


def create_two_adjacent_aas():
    cursor.execute(
        """
        CREATE TABLE `two_adjacent_AAs` (
            `first_aa_name`	TEXT,
            `second_aa_name`	TEXT,
            `combined_aas`	TEXT,
            `CA_CN_angle`	REAL,
            'CN_CA_angle'   REAL,
            `p_CN_angle`	REAL,
            `pp_angle`	REAL,
            `protein_id`	TEXT,
            `chain`	TEXT,
            `first_aa_sequence`	INTEGER,
            `CN_a` REAL,
            `CN_b` REAL,
            `CN_c` REAL,
            `NCA_a` REAL,
            `NCA_b` REAL,
            `NCA_c` REAL,
            `CAC_a` REAL,
            `CAC_b` REAL,
            `CAC_c` REAL
        );
        """
    )
    connect.commit()


def create_two_aas_analysis():
    """
    vector1 = N-CA in first protein in absolute value
    vector2 = CA-C in first protein in absolute value
    vector3 = (0, 0, 1)
    C_N = a * vector1 + b * vector2 + c * vector3
    :return:

    """


def create_one_aa():
    cursor.execute(
        """
        
        """
    )


def create_four_AAs():
    pass

# create_protein_structure()
create_two_adjacent_aas()