from matplotlib import colors as col
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import random
import os.path

# instead of tracing back on itself, turn by a factor of "turn" (potentially try other values?)
turn = 1

# display 3D plot base
fig = plt.figure()
ax = plt.axes(projection='3d')

# aa_num: a[#] = amino acid letter (# of amino acid -> letter of amino acid)
aa_num = []
# num_loc: a[#] = location of amino acid face (# of amino acid -> location of amino acid face)
num_loc = []

# Dictionary of Colours REMAINS CONSTANT FOR ALL RUNS
aa_col = {}
# taking in face colours for each amino acid (remains constant at all times)
f_col = open("colours_locations", "r"); next(f_col)
for i in range(0, 20):
    col_amino = f_col.readline().split(' ')
    aa_col[col_amino[0]] = col_amino[1]  # take in colour assignment from file
f_col.close()

# names for relative fasta files
# FASTA_FILES = {"a3c8": "samples/a3c8/5SV4_A.fasta.txt",
#                "a3d": "samples/a3d/2A3D_A.fasta.txt",
#                "albumin": "samples/albumin/1AO6_A.fasta.txt",
#                "apo": "samples/apo/1AKP_A.fasta.txt",
#                "fv": "samples/fv/1AP2_A.fasta.txt",
#                "insulin": "samples/insulin/2LWZ_A.fasta.txt",
#                "KAMP": "samples/KAMP/5KI0_A.fasta.txt",
#                "lactal": "samples/lactal/1A4V_A.fasta.txt",
#                "monelin": "samples/monelin/1IV7_A.fasta.txt",
#                "ovalbumin": "samples/ovalbumin/3VVJ_A.fasta.txt",
#                "scm": "samples/scm/1FUW_A.fasta.txt",
#                "tfiia": "samples/tfiia/5M4S_A.fasta.txt",
#                "xiv": "samples/xiv/1B9P_A.fasta.txt"}
FASTA_FILES = {"human insulin A": "samples/3I3Z/3I3Z_A.fasta.txt",
               "human insulin B": "samples/3I3Z/3I3Z_B.fasta.txt"}

# 3D vector class
class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        return Vector(self.x + other.x,
                      self.y + other.y,
                      self.z + other.z)

# aminoassign: reads in amino assignment from file in /assignments/assign_input
# updates global aa_num and num_loc
def aminoassign(assign_file):
    global aa_num, num_loc

    # reset to 0's
    aa_num = [0] * 20
    num_loc = [0] * 20

    # face assignments
    if (assign_file == ""):
        f_assign = open("assignments/1", "r")  # Default Random
    else:
        f_assign = open("assignments/" + assign_file, "r")

    # face locations
    f_loc = open("colours_locations", "r")
    for i in range(23):
        next(f_loc)

    #  read in aa_loc values from files
    for i in range(0, 20):
        cur_amino = f_assign.readline().split()[0]  # take in aa assignment from file
        loc = [float(x) for x in f_loc.readline().split()]  # take in face location from file
        loc_v = Vector(loc[0], loc[1], loc[2])  # make vector

        aa_num[i] = cur_amino
        num_loc[i] = loc_v

    # close files (just to be clean)
    f_assign.close()
    f_loc.close()


# plot_sequence: plots the amino acid sequence using icosahedron approach
#   seq: some iterable structure of strings of valid amino acids according to face_assignments.txt
#   centre: the centre plotting should begin at, usually Vector(0,0,0) for first run
#   speed: speed to draw at, =0 if want to skip draw
#
#   if the last amino acid face number + the current amino acid = 19, they are on opposite faces of the shape
#       this would cause a weird shape, therefore preforms a slight turn
#
#   returns: the position of the last amino acid
def plot_sequence(seq, centre, prev_face, speed):

    for amino in seq:
        try:
            face = aa_num.index(amino)
        except ValueError:    # error catching and skipping invalid amino acids
            print("**********ERROR '" + amino + "' is not a valid amino acid and is skipped")
            continue

        # if opposite amino acids, make a dip before turning back around
        if (face + prev_face) == 19:
            turn_vec = Vector(centre.x, centre.y + turn, centre.z)

            plt.plot([centre.x, turn_vec.x],
                     [centre.y, turn_vec.y],
                     [centre.z, turn_vec.z],
                     aa_col[amino])

            centre = turn_vec

        new = centre + num_loc[face]

        plt.plot([centre.x, new.x],
                 [centre.y, new.y],
                 [centre.z, new.z],
                 aa_col[amino])  # plot previous centre to new amino acid, according to amino acid colour

        # pause speed if needed (not 0)
        if speed:
            draw()  # draw the line
            pause(speed)  # at this speed? or for this amount of time? something, causes delay. smaller = faster

        centre = new
        prev_face = face

    return [centre, prev_face]


# plotform: takes in the format and the speed and runs one instance of the turtle algortihm
#
#       form: can be 'c', 'r' or else (file name)
#           'c': command line input of sequence, for quick tests
#           'r': generates random sequence and plots it
#           else: preforms on sequence given in fasta file name: ** NOTE ** skips first line for description of file
#
#       plots a green diamond at the initial position and a red diamond at the final position
def plotform(form, speed):

    # Initial Position
    plt.plot([0], [0], [0], 'gD')  # plot start position of sequence: GREEN DIAMOND

    prev_face_default = -1
    centre_default = Vector(0, 0, 0)

    #  Command Line Input
    if form == 'c' or form == 'C':
        sequence = input("Enter protein sequence\n")
        plt.title(sequence)
        final = plot_sequence(sequence, centre_default, prev_face_default, speed)[0]

    #  Random Input
    elif form == 'r' or form == 'R':
        n = input("Enter number of random amino acids\n")
        sequence = ""

        for i in range(0, int(n)):
            rand_amino = random.randrange(20)
            sequence += aa_num[rand_amino]

        print("Sequence: " + sequence)
        plt.title("random input")
        final = plot_sequence(sequence, centre_default, prev_face_default, speed)[0]

    #  Filename Input
    else:
        try:
            sequence = open(form, "r")
        except FileNotFoundError:
            print("**********ERROR: Filename incorrect")
            exit()

        sequence.readline()  # skip first line
        result = [centre_default, prev_face_default]

        # plot each line in the file
        for line in sequence:
            if line[0] == '>':
                continue
            result = plot_sequence(line.strip('\n'), centre_default, prev_face_default, speed)
            centre_default = result[0]
            prev_face_default = result[1]

        final = result[0]

    # Final position
    plt.plot([final.x], [final.y], [final.z], 'rD')  # plot final position of sequence: RED DIAMOND


# for the rotated gif animation
def rotate(angle):
    ax.view_init(azim=angle)


#  turtle_main: runs the main turtle generation algorithm
#
#        command_input: 's', 'a', or ''
#          s: does all the assignments (1-9 in /assignments) for one sample,
#                  according to the name given in sample_input
#                  this won't display the images, and will create rotating gif images in samples/sample_input
#          a: does all the samples (in the /samples folder,
#                  according to global "fasta_files" dictionary on the assignment given in assign_input
#                  this won't display the images, and will create a folder assign_input/ to put the rotating gif images
#          '': to run one sample and display the graph in python,
#                  needs the assignment file (assign_input) with default 1,
#                  needs the pause speed with default 0,
#                  needs the fasta_file path in sample_input with default "sequence"
#              if enter 'c' for the fasta file path, can enter command line input for quick tests
#                                          still needs the assign_input and speed_input fields
#              if enter 'r' for fasta file path, generates random amino acid sequences
#                                          still needs the assign_input and speed_input fields
#
#       assign_input: name of assignment file, inside the "assignments" folder. (only needed for 'a' and '')
#
#       speed_input: value of pause speed, to watch the amino acid be drawn. 0 to not pause at all
#
#       sample_input: either the name of the sample for 's' or FASTA file name path for ''
def turtle_main(command_input, assign_input, speed_input, sample_input):

    if command_input == "s":
        # run for each assignment
        for i in range(1,10):
            ind = str(i)
            print("working on number "+ind)

            assign_input = ind

            # create path to save at, skip if already exists (no need to replace)
            file_name = ("samples/" + sample_input + "/" + assign_input + "_" + sample_input + ".gif")
            if(os.path.exists(file_name)):
                print("already done number " + ind)
                continue

            # create assignment, run plot
            aminoassign(assign_input)
            plotform(FASTA_FILES[sample_input], speed_input)

            # save 3d rotation of plot in gif format
            rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0,362,2), interval=100)
            rot_animation.save(file_name, dpi=80, writer='imagemagick')

            print("done number "+ind)

            ax.clear()

    elif command_input == "a":

        aminoassign(assign_input)

        for chain, file in FASTA_FILES.items():
            print("working on " + chain)

            # create path to save at, skip if already exists (no need to replace)
            file_name = (file[:file.find('.')] + '_' + assign_input + '.gif')
            if (os.path.exists(file_name)):
                print("already done " + chain)
                continue

            # create assignment, run plot
            plotform(file, speed_input)

            # save 3d rotation of plot in gif format
            rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
            rot_animation.save(file_name, dpi=80, writer='imagemagick')
            print("done " + chain)

            ax.clear()

    else: # nothing, C or R. simple: just one run

        # RUN PROGRAM
        aminoassign(assign_input)  # generate amino assignment
        plotform(sample_input, speed_input)  # plot
        plt.show()  # show
        ax.clear()


# taking in user input, only if this is the main (so runtest can run).
# description of the input can be found above turtle_main(...)
if __name__ == '__main__':
    # command = input('For user control, input nothing\n'
    #         'To do all samples for one assignment, enter "a", '
    #         'To do all assignments for one sample, enter "s"\n')
    #
    # if command == "a" or command != "s":    # (a,'')
    #     assign = input("Enter amino acid assignment file in assignments/ folder:\n")
    # else:
    #     assign = ""
    #
    # speed = input("Enter pause speed of plot (recommended 0.1-0.01) (if needed):\n")
    # if speed == "":
    #     speed = 0
    # else:
    #     speed = float(speed)
    #
    # if command == "s":
    #     sample = input('Enter the sample name (s)\n')
    # elif command != "a":
    #     sample = input('Enter the FASTA file name path ('')\n')
    # else:
    #     sample = ""
    #
    # turtle_main(command, assign, speed, sample)
    turtle_main('a', '1', 0.1, '')
