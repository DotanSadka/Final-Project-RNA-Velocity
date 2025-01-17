import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

CELL_TYPE = 'T-cell_CD8A'
FILTER = 0
T = 0
IS_FILTER = True


def plot_points_arrows(points_start, points_end, pc1, pc2):

    salmon_color = '#CCCCFF'  # CD8A
    golden = '#FFD700'  # CD3D
    turquoise = '#40E0D0'  # CD3E
    light_green = '#90EE90'  # CD4

    # Plot the points
    for point in points_start:
        plt.plot(point[0], point[1], marker='o', markersize=5, color=salmon_color,
                 markeredgecolor='black')

    # Plot the arrows
    for start, end in zip(points_start, points_end):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        plt.arrow(start[0], start[1], dx, dy, color='black', width=0.008, head_width=0.15, length_includes_head=True)

    # Set the aspect ratio
    plt.axis('equal')


    # Set the titles
    #if IS_FILTER:
    #    title = 'With Filter'
    #else:
    #    title = 'No Filter'

    #plt.xlim(-5, 5)
    #plt.ylim(-6.5, 5)

    plt.title("{} : t={}".format(CELL_TYPE, T))
    plt.xlabel(f"{pc1}")
    plt.ylabel(f"{pc2}")

    # # Show the plot
    # plt.show()
    plt.savefig(f'{CELL_TYPE}_scale/{pc1} {pc2}/{pc1} {pc2} - t={T}_ZO.png')
    plt.clf()


def load_s0_points(file_name):
    # Load CSV file into a dataframe
    df = pd.read_csv(file_name)

    # Filter the dataframe based on the 'label_c' column
    filtered_df = df[df['labels_c'] == CELL_TYPE]

    # Access the lists of data by column name todo: check if needs to cast into float
    pc1 = filtered_df['pc1']
    pc2 = filtered_df['pc2']

    coordinates = list(zip(pc1, pc2))

    return coordinates


def load_pcs(file_name, t):
    # Load CSV file into a dataframe
    df = pd.read_csv(file_name)

    # Filter the dataframe based on the 'label_c' column
    filtered_df = df[df['labels_c'] == CELL_TYPE]
    #if t:
    #filtered_df.to_csv(f'{CELL_TYPE}_files/{CELL_TYPE}_{t}.csv')

    # Access the lists of data by column name
    pc1 = filtered_df['pc1']
    pc2 = filtered_df['pc2']
    pc3 = filtered_df['pc3']

    return pc1, pc2, pc3

def load_vec_sum(file_name):
    # Load CSV file into a dataframe
    df = pd.read_csv(file_name)

    # Access the lists of data by column name
    sum1 = df['t=1']
    sum2 = df['t=2']
    sum3 = df['t=3']
    sum4 = df['t=4']

    return [sum1, sum2, sum3, sum4]


if __name__ == '__main__':
    #points_start = load_s0_points("cell_pc_df.csv")
    # vec_sums = load_vec_sum(f"FILTERATION-0.01/S_VEC_SUM.csv")
    pc1_start, pc2_start, pc3_start = load_pcs("cell_pc_df.csv", T)

    #for filter in [0.01]:
        #FILTER = filter
        #for bool in [False, True]:
            #IS_FILTER = bool
    #for t in np.arange(0.2, 5.2, 0.2):
    for t in [1.4]:
        #T = round(t,1)
        T=t
        #if IS_FILTER:
        file_name = f"st_filter/st_t={T}_filter.csv"
        #else:
        #file_name = f"st_t{T}_nofilter.csv"

        pc1_end, pc2_end, pc3_end, = load_pcs(file_name, T)
        plot_points_arrows(list(zip(pc1_start, pc2_start)), list(zip(pc1_end, pc2_end)), 'PC1', 'PC2')
        #plot_points_arrows(list(zip(pc1_start, pc3_start)), list(zip(pc1_end, pc3_end)), 'PC1', 'PC3')
        #plot_points_arrows(list(zip(pc2_start, pc3_start)), list(zip(pc2_end, pc3_end)), 'PC2', 'PC3')


