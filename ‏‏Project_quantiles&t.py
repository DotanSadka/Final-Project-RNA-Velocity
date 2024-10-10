# for each cell we have an initial state of its spliced RNA for a gene (s(t) = s0 + vt s.t. t = 0 => s = s0)
# plot embedding - PCA
# calculate all gammas for each gene using knn pooling (gamma = u/s) (seems like gammas are calculated from s0 and u0)
# calculate v for each gene according ds/dt = v = u(t) - gamma*s(t)
# for each cell we calculate the new s(t) and then get the trajectory

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression


# HYPER PARAMETERS

FIRST_RUN = False # if you dont have "neighbors.csv" yet- True. else - False
K = 200
t = 0
PERCENTAGE = 0.1
df_normal_st = None
df_filtered_st = None
GROUND_TRUTH = ['ACTA2', 'ACTR3B', 'ANLN', 'AP1M2', 'APOE', 'AR', 'BCL2', 'BGN', 'BRAF', 'C10orf54', 'C1QA', 'CCL4', 'CCNB1', 'CCND1', 'CCNE1', 'CD14', 'CD24', 'CD36',
                'CD38', 'CD3E', 'CD3G', 'CD5', 'CD63', 'CD69', 'CD7', 'CD79B', 'CD8A', 'CDC20', 'CDC6', 'CDH3', 'CDK4', 'CDK6', 'CDKN2A', 'CENPF', 'COL1A1', 'COL3A1', 
                'COL4A1', 'COL4A5', 'CRABP2', 'CSRP2', 'CSTB', 'CTSL', 'CXXC5', 'DCN', 'EFNA5', 'EIF3E', 'ELF5', 'ERBB2', 'ERBB3', 'ESR1', 'FABP7', 'FAP', 'FASN', 'FAT1', 
                'FCGR3A', 'FCN1', 'FCRL5', 'FGFR2', 'FGFR4', 'FN1', 'FOS', 'FOXC1', 'FOXP3', 'FTL', 'GATA3', 'GNLY', 'GPNMB', 'GRB7', 'GSN', 'GZMB', 'HIF1A', 'HLA-C', 
                'HLA-DRA', 'HLA-DRB5', 'HLA-E', 'HSPB1', 'HSPG2', 'ICOS', 'IFITM3', 'IGF1R', 'IGFBP5', 'IGHG1', 'IGHG4', 'IGHM', 'IGKC', 'IL32', 'IL7R', 'ISG20', 'JUN', 
                'KIF23', 'KIT', 'KRAS', 'KRT14', 'KRT15', 'KRT18', 'KRT19', 'KRT5', 'KRT8', 'LAMA1', 'LDB2', 'LGALS2', 'LGMN', 'LRP2', 'LST1', 'LTF', 'LYZ', 'MAPK13', 
                'MAPK3', 'MDM2', 'MELK', 'MLPH', 'MMP12', 'MMP9', 'MS4A1', 'MSR1', 'MT2A', 'MYB', 'MYLK', 'MYO10', 'MYO5B', 'MZB1', 'NFKBIA', 'NKG7', 'NOTCH1', 'NR3C1', 
                'NUF2', 'ORC6', 'PABPC1', 'PDPN', 'PGR', 'PHGDH', 'PLVAP', 'POU2AF1', 'PRLR', 'PTPRB', 'PTPRC', 'PTTG1', 'RPSA', 'S100A14', 'S100A9', 'SFRP1', 'SIAH2', 
                'SKAP1', 'SLC2A1', 'SLC39A6', 'SLPI', 'SOX10', 'SOX18', 'SOX4', 'SPDEF', 'SULF1', 'TCF4', 'TCL1A', 'TFF1', 'TFF3', 'TFRC', 'THEMIS', 'THY1', 'TMEM45B', 
                'TMSB10', 'TMSB4X', 'TPM2', 'TSPAN1', 'UBE2T', 'XCL1'] #165

def save_to_csv(file_name, header, row_names, values):
    df = pd.DataFrame(values, columns=header, index=row_names)
    
    # Save the data frame to CSV
    df.to_csv(file_name, index=True)

def plot_precision_recall():
    # CSV file path
    csv_file = "Hyperparameters/gene/F-score/precision_recall_unspliced.csv"

    # Read data from CSV file
    expression_levels = np.arange(0, 0.041, 0.001)
    

    with open(csv_file, mode='r') as file:
        lines = file.readlines()
        header = lines[0].strip().split(",")  # Read header and split into column names

        # for line in lines[1:]:
        line_p = lines[1]
        line_r = lines[2]

        values_p = line_p.strip().split(",")[1:]  # Split each line into values
        values_r = line_r.strip().split(",")[1:]  # Split each line into values

        precision_values = [float(x) for x in values_p]
        recall_values = [float(x) for x in values_r]

    # Calculate F-score
    f_scores = [(2 * precision * recall) / (precision + recall) for precision, recall in
                zip(precision_values, recall_values)]

    max_f_score = max(f_scores)
    max_f_score_index = f_scores.index(max_f_score)
    max_f_score_expression_level = expression_levels[max_f_score_index]

    # Plot precision, recall, and F-score
    plt.figure(figsize=(8, 6))
    plt.plot(expression_levels, precision_values, label='Precision')
    plt.plot(expression_levels, recall_values, label='Recall')
    plt.plot(expression_levels, f_scores, label='F-score')

    plt.scatter(max_f_score_expression_level, max_f_score, color='red')  # הוספת הנקודה האדומה
    plt.text(max_f_score_expression_level, max_f_score, 
             f'Max F-score: {max_f_score:.2f}\n@Expression: {max_f_score_expression_level:.4f}',
             fontsize=10, ha='right', color='red') 

    plt.xlabel('Expression Levels')
    plt.ylabel('Precision / Recall / F-score')
    plt.title('Precision, Recall, and F-score vs. Expression Levels')
    plt.legend()
    plt.grid(True)
    plt.savefig('Hyperparameters/gene/F-score/Precision_Recall_F-score_plot.png')
    plt.show()

def precision_and_recall(cells, l):
#def precision_and_recall(cells):
    """
    Precision: what predicted True Positive divided by all Positive Predictions - we want the lowest fp meaning catch less not related genes in PG
    Recall: what predicted True Positive divided by all Positive samples - we want the lowest fn meaning catch more GT genes
    """

    tp = 0
    fp = 0

    # PG - predicted genes
    predicted_genes = [gene.name for gene in cells.genes if gene.filter == False]
    for gene in predicted_genes:
        if gene in GROUND_TRUTH:
            print(f'{tp} gene:\t{gene}')
            tp += 1
        else:
            fp += 1
    print(f'fp: {fp}')
    fn = len(GROUND_TRUTH) - tp

    # Calculate precision and recall
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    vec = [precision, recall]
    cells.precision_recall_df[l] = vec 

    print("Precision: {}\tRecall: {}".format(precision, recall))
    #exit(0)


class CellsDict(dict):
    """
    Dict of cell type as keys and list of cells of that cell type as values
    """

    def __init__(self, df):
        super(CellsDict, self).__init__()
        self.df = df
        self.number_of_cells = self.df.shape[0]
        spliced_df = pd.read_csv("conbined_counts_0_271Genes.csv")
        unspliced_df = pd.read_csv("conbined_counts_1_271Genes.csv")
        if FIRST_RUN:
            self.neighbors_df = None
        else:
            self.neighbors_df = pd.read_csv('neighbors.csv')

        self.genes = []
        self.genes_name = []

        self.S = self.initiate_mRNA_expression(spliced_df, True)
        self.U = self.initiate_mRNA_expression(unspliced_df)
        
        self.expression_levels = np.arange(0, 0.041, 0.001)

        self.precision_recall_df = pd.DataFrame(columns=list(self.expression_levels))
        row_names = ['precision', 'recall']
        self.precision_recall_df = self.precision_recall_df.reindex(row_names)

    def initiate_mRNA_expression(self, df, is_first=False):
        """
        
        Initiate all expression values of all genes in all cells.
        :param df: dataframe which contains expression values
        :return: dict that holds all cells' names as keys and each cell initiate expression level of all its genes
                 as values, such that values are a dict of a gene name - key and it's expression val - value.
        """
        cells_name = list(df["cell_name"])
        genes_name = list(df.columns)[1:]
        values = df.values

        expression_dict = {}
        for i in range(self.number_of_cells):
            #dict: key- cell's name. valueis dict: key- gene's name. value- gene expression.
            expression_dict.update({cells_name[i]: {gene: value for gene, value in zip(genes_name, values[i, 1:])}}) 

        if is_first:
            self.init_genes(genes_name)

        return expression_dict

    def init_genes(self, genes_name):
        self.genes_name = genes_name
        for name in genes_name:
            self.genes.append(Gene(name))

    def set_dict(self):
        """
        Set the cells to be clustered according to their cell type
        """

        cells_name = list(self.df['Row.names'])
        cells_type = list(self.df["labels_c"])

        self.cell_types = list(set(cells_type))

        pc1 = self.df["pc1"]
        pc2 = self.df["pc2"]
        pc3 = self.df["pc3"]
        for i in range(self.number_of_cells):
            cell_name = cells_name[i]
            self.update({cells_name[i]: Cell(cells_name[i], self.S[cell_name], self.U[cell_name],
                                                 cells_type[i], pca=np.array([pc1[i], pc2[i], pc3[i]]))})

    def save_distances_to_csv(self):
        cells_neighbors = {}
        for current_cell in self.values():
            distances = []
            for i, cell in enumerate(self.values()):
                # Calculate Euclidean distance to all points in PCA space
                distance = np.linalg.norm(current_cell.pca - cell.pca, axis=0)
                distances.append((distance, cell))

            # sort the list based on the first value in each tuple
            sorted_list = sorted(distances, key=lambda x: x[0])

            # get all names in the order that is correlated to their distance
            sorted_cells_names = [x[1].name for x in sorted_list]
            cells_neighbors.update({current_cell.name: sorted_cells_names})

        df = pd.DataFrame(cells_neighbors)
        df.to_csv('neighbors.csv', index=False)
        self.neighbors_df = pd.read_csv('neighbors.csv')

    def calculate_gammas(self):
    #def calculate_gammas(self, l): #for the Fscore calculation
        """
        function calculates gamma for each gene and update all genes accordingly
        """
        gammas = []
        for gene in self.genes:
            #if gene.name in GROUND_TRUTH:
                points = gene.cells_coordinates.values()  # get all (x=s_knn, y=u_knn) of all cells
                x_values = [point[0] for point in points]  # get all x values from points
                y_values = [point[1] for point in points]  # get all y values from points

                """
                FILTRATION RULE
                1. brfore F-score calculation- no filter
                2. for F-score calculation - if max(x_values) < l:/ if max(y_values) < l:
                3. after F-score calculation - put the vale in the if statment.
                """
                
                if max(x_values) < 0.004 or max(y_values) < 0.005 or gene.name not in GROUND_TRUTH:
                #if max(x_values) < l: 
                    gene.filter = True
                
                ## Find gamma ##

                # Find the extreme points on the top right and bottom left
                max_x, max_y = max(x_values), max(y_values)
                min_x, min_y = 0, 0

                # Get the number of points for specific percentage of the total number of points
                num_points = len(points)  # supposed to be 4530
                num_points_percent = int(num_points * PERCENTAGE)

                # Sort the points by distance from the top right point
                sorted_points = sorted(points, key=lambda point: ((point[0] - max_x) ** 2 + (point[1] - max_y) ** 2) ** 0.5)

                # Get the top percentage of the sorted points
                top_points = sorted_points[:num_points_percent]

                # Sort the points by distance from the bottom left point
                sorted_points = sorted(points, key=lambda point: ((point[0] - min_x) ** 2 + (point[1] - min_y) ** 2) ** 0.5)

                # Get the bottom 5% of the sorted points 

                bottom_points = sorted_points[:num_points_percent]

                points = top_points + bottom_points
                x_values = np.array([point[0] for point in points]).reshape(-1, 1)  # get all x values from points
                y_values = np.array([point[1] for point in points])  # get all y values from points

                # calculate the slope for the gamma
                reg = LinearRegression(fit_intercept=False)
                reg.fit(x_values, y_values)
                slope = reg.coef_

                # update gamma
                gene.set_gamma(slope[0])
                gammas.append(slope[0])

        return gammas

    def set_gammas_for_cells(self, gammas):
        for cell in self.values():
            cell.set_gammas(gammas)


class Gene:
    """
    Gene class.
    Each gene has 2 dict of all the spliced and unspliced RNA expression (as value) in a current cell (as key)
    :param self.filter - flags if should filter in v calculation
    """

    def __init__(self, name):
        self.name = name
        self.gamma = None
        self.cells_coordinates = {}
        self.filter = False

    def set_gamma(self, slope):
        self.gamma = slope

    def get_gamma(self):
        return self.gamma

    def set_coordinates(self, cell_name, coordinates):
        """
        Set the cell's coordinates for that specific gene - spliced is x-axis and unspliced is y-axis
        :param cell_name: the cell that gene coordinates were calculated for
        :param coordinates: the coordinates from the s_knn and u_knn of the cell
        """
        self.cells_coordinates.update({cell_name: coordinates})


def create_boolean_vector(filter_indices, vector):  # todo: change func name
    # Initialize the boolean vector with all zeros
    # boolean_vector = [1] * length

    # Set the corresponding positions to 1
    for index in filter_indices:
        vector[index] = 0

    return vector


class Cell:
    """
    Cell class.
    Each cell has a dict of all its genes such that gene name is the key and the value is a tuple of gene expression (s, u)
    """

    def __init__(self, name, S0, U0, cell_type, pca=None):
        """
        :param name: The cell name.
        :param S0: A dict contains cell's initial amount of spliced mRNA (values) for each gene (gene objects - keys)
        :param s0_norm: normalized s0 dict - divide s0 by all spliced count of current cell - use for velocity calculation
        :param total_S_counts: The amount of Spliced counts of current cell
        :param s_knn: A dict like S0, but normalized according to K nearest neighbors (divide by sum of neighbors spliced counts) - use for gamma calculation
        :param U0: A dict contains cell's initial amount of unspliced mRNA (values) for each gene (gene objects - keys)
        :param u0_norm: normalized u0 dict - divide u0 by all unspliced count of current cell - use for velocity calculation
        :param total_U_counts: The amount of Unspliced counts of current cell
        :param u_knn: A dict like U0, but normalized according to K nearest neighbors (divide by sum of neighbors unspliced counts) - use for gamma calculation
        :param cell_type: What cluster the cell belongs to.
        :param neighbors: List of all cell's neighbors (according to the definition of neighbors that is used such as
                                                        celltype, PCA embedding, etc)
        :param pca: vector (pc1, pc2, pc3) of current cell.
        """
        self.name = name
        self.S0 = S0
        self.total_S_counts = self.get_total_counts(self.S0)
        self.s0_norm = self.normalization(self.S0, self.total_S_counts)
        self.s_knn = {}
        self.U0 = U0
        self.total_U_counts = self.get_total_counts(self.U0)
        self.u0_norm = self.normalization(self.U0, self.total_U_counts)
        self.u_knn = {}
        self.cell_type = cell_type
        self.neighbors = None

        self.pca = pca
        self.gammas = None

    def set_gammas(self, gammas):
        self.gammas = gammas

    def get_st(self, cells, t):
        """
        :param t: The time we want to predict the S in the cell.
        :return: A vector with all the spliced amount of each gene of the cell in time t.
        """
        v = self.get_v_calculation(cells)

        normal_st = np.array(list(self.s0_norm.values())) + v[0] * t
        filtered_st = np.array(list(self.s0_norm.values())) + v[1] * t
        return normal_st, filtered_st


    def get_v_calculation(self, cells):
        """
        Calculate V according to V=ds/dt which is the constant u - gamma * s - o.
        :return: vector size 297 which contains all the velocities for each gene in cell.
        """
        # if self.gammas is None:
        #     raise ValueError("Gammas attribute is not set.")

        normal_v = np.array(list(self.u0_norm.values())) - self.gammas * np.array(list(self.s0_norm.values()))

        filter_indices = [index for index, gene in enumerate(cells.genes) if gene.filter == True] #יוצר מערך של האינדקסים של הגנים המפולטרים
        filtered_v = create_boolean_vector(filter_indices, normal_v.copy())
        #print(normal_v, filtered_v)
        return normal_v, filtered_v


    def initialize_neighbors(self, cells):
        """
         Initialize the  cell's neighbors according to PCA method.
        """
        self.neighbors = self.get_PCA_neighbors(cells)


    def get_total_counts(self, expression_dict):
        """
        Calculate the total count of mRNA in current cell
        :param expression_dict: S0 or U0
        :return: the total count
        """
        total = 0
        for x in expression_dict.values():
            total += x

        return total

    def normalization(self, expression_dict, total_counts):
        """
        normalize S0 and U0 for knn pooling
        :param expression_dict: S0 or U0
        :param total_counts: counts of spliced or unspliced.
        :return: s0_norm or u0_norm
        """
        dict = {gene: expression / total_counts for gene, expression in expression_dict.items()}
        return dict

    def substitute_S_by_knn(self):
        """
        Calculate new S0 according to cell's neighbors
        :return: the new S0 according to the knn pooling
        """
        sum_vec = 0  # holds the sum of all neighbors spliced counts for each gene
        total_neighbor_S_counts = 0  # holds the total splice count of all neighbors
        for neighbor in self.neighbors:
            sum_vec += np.array(list(neighbor.S0.values()))
            total_neighbor_S_counts += neighbor.total_S_counts

        # sum of vector is 1 - normalization for the vector
        sum_vec = sum_vec / total_neighbor_S_counts

        dict = {gene: sum_vec[i] for i, gene in enumerate(self.S0.keys())}
        #dict = {gene: sum_vec[i] if gene in GROUND_TRUTH else 0 for i, gene in  enumerate(self.S0.keys())}
        self.s_knn = dict

    def substitute_U_by_knn(self):
        """
        Calculate new U0 according to cell's neighbors
        :return: the new U0 according to the knn pooling
        """
        sum_vec = 0  # holds the sum of all neighbors unspliced counts for each gene
        total_neighbor_U_counts = 0  # holds the total ubsplice count of all neighbors
        for neighbor in self.neighbors:
            sum_vec += np.array(list(neighbor.U0.values()))
            total_neighbor_U_counts += neighbor.total_U_counts

        # sum of vector is 1 - normalization for the vector
        sum_vec = sum_vec / total_neighbor_U_counts

        dict = {gene: sum_vec[i] for i, gene in enumerate(self.U0.keys())}
        #dict = {gene: sum_vec[i] if gene in GROUND_TRUTH else 0 for i, gene in enumerate(self.U0.keys())}
        self.u_knn = dict

    def get_PCA_neighbors(self, cells):
        """
        According to PCA calculation - find k nearest neighbors using euclidean distance metric.
        :param cells: data structure which holds all cells
        :return: k neighbors of current cell
        """
        global K
        names = list(cells.neighbors_df[str(self.name)].head(K))
        neighbors = []
        for i in range(K):
            neighbors.append(cells[names[i]])

        return neighbors




def main():
    global PERCENTAGE, t

    #k = [70, 100, 200, 250 ,300,350,400,450] #run this to check the quantiels (loop on this in line 511)
    k = [200]
    gammas_lists = []
    for i in k:
        global K
        K = i
        cells = CellsDict(pd.read_csv("cell_pc_df.csv"))
        cells.set_dict()

        # todo: if first time running code - run this to calculate neighbors of each cell
        if FIRST_RUN:
            cells.save_distances_to_csv()
        # return

        for cell in cells.values():
            # initialize all neighbors and calculate new s_knn and u_knn accordingly
            cell.initialize_neighbors(cells)
            cell.substitute_S_by_knn()
            cell.substitute_U_by_knn()

            # set the knn_normalization (the coordinates) of each gene for each cell
            for j, (x, y) in enumerate(zip(cell.s_knn.values(), cell.u_knn.values())):
                cells.genes[j].set_coordinates(cell, (x, y))

        """
        Calculate gamma for all genes and set gammas for each cell in order to calculate velocity and s(t)
        """


        #for l in cells.expression_levels:
        #    gammas = cells.calculate_gammas(l)
        #    precision_and_recall(cells, l)  # todo: used for hyperparameter fine-tuning
        gammas = cells.calculate_gammas()
        #precision_and_recall(cells) 
        #cells.precision_recall_df.to_csv(f'precision_recall_unspliced.csv', index=True)
        #plot_precision_recall()
        gammas = np.array(gammas)
        gammas_lists.append(gammas)
    
    #genes = [g.name for g in cells.genes if g.filter == False] #gene that passed filtration for validation check
    #print(genes)

    #genes = GROUND_TRUTH
    #save_to_csv(f'quantiles_70-500_{PERCENTAGE}.csv', header=genes, row_names=k, values=gammas_lists)
    #stop here for the quantiles sensitivity check

    
       
    # First row
    #T = np.arange(0.2, 5.2, 0.2)
    T =[1.2]
    for r in T:
        global t
        t = r


        global df_normal_st, df_filtered_st
        df_normal_st = pd.DataFrame({'cells_name': cells.genes_name})


        df_filtered_st = pd.DataFrame({'cells_name': cells.genes_name})
        df = pd.DataFrame({'cells name': ['s0', 's(t)', 's(t)_filtered']})


        for cell in cells.values():
            cell.set_gammas(gammas)
            s0_norm = np.array(list(cell.s0_norm.values()))
            #print(f's0: {s0_norm}\ts(t): {cell.get_st(1)}\n')

            st_values = cell.get_st(cells, t)  # first value: normal s(t), second value: filtered s(t)
            sum_st_normal = np.sum(st_values[0])
            sum_st_filtered = np.sum(st_values[1])
            sum_s0_norm = np.sum(s0_norm)

            # print(f'cell name: {cell.name}\ts0_sum: {sum_s0_norm}\ts(t)_sum: {sum_st_normal}\t'
            #       f's(t)_filtered_sum: {sum_st_filtered}\n')
            df = pd.concat([df, pd.DataFrame({cell.name: [sum_s0_norm, sum_st_normal, sum_st_filtered]})], axis=1)
            df_normal_st = pd.concat([df_normal_st, pd.DataFrame({cell.name: st_values[0]})], axis=1)
            df_filtered_st = pd.concat([df_filtered_st, pd.DataFrame({cell.name: st_values[1]})], axis=1)

        df.to_csv(f"t={t}_check.csv", index= False)  # csv file to check sum of new vector s_t
        #df_normal_st.transpose().to_csv(f"FILTERATION-0.1/csv_for_plots/t={t}_Gene_expression.csv")  # the new s_t matrix


        #df_filtered_st.transpose().to_csv(f"Hyperparameters/t/geneExpAndT/t={t}_Gene_expression_filtered.csv")
        df_filtered_st.transpose().to_csv(f"t={t}_Gene_expression_filtered_check.csv")



def get_top_st_sum():
    df = pd.read_csv(f"FILTERATION/t={t}.csv")
    row = df.iloc[1][2:]
    sorted_row = row.sort_values(ascending=False)

    print(f'Top value: {sorted_row.head(10)}\tBottom value: {sorted_row.tail(10)}\n')
    print(f'max:{max(row)} min:{min(row)}')


def calculate_mean_st_sum(path):
    df = pd.read_csv(path)
    row1 = df.iloc[1][2:]
    row2 = df.iloc[2][2:]

    #print(f'ST mean sum value: {np.mean(row1)}\tST Filtered mean sum value: {np.mean(row2)}\n')
    #print(f'max:{max(row1)} min:{min(row1)}\n')
    #print(f'max:{max(row2)} min:{min(row2)}\n')
    return np.mean(row2)


def t_iteration_to_find_mean():
    T = np.arange(0.2, 5.2, 0.2)
    rows =[]
    for r in T:
        global t
        t = round(r, 1)
        path = f"Hyperparameters/t/geneExpAndT/t={t}.csv"
        mean_st_sum = calculate_mean_st_sum(path)
        rows.append({"t": t, "mean_s(t)_sum": mean_st_sum})
    return rows

if __name__ == "__main__":

    main()
    """
    # Run s(t) calculation

    df = t_iteration_to_find_mean()
    df = pd.DataFrame(df)
    df.to_csv(f"Hyperparameters/t/geneExpAndT/mean_s(t)_sum.csv", index=False)
    
    # used to check which "t" is best - got convergence around t=3
    for i in range(1, 21):
        t = i
        path = f"FILTERATION-0.1/csv_for_plots/t={t}.csv"
        get_top_st_sum()
        df = pd.DataFrame(columns=["t", "mean_s(t)_sum"])
        mean_st_sum = calculate_mean_st_sum(path)
        rows =[]
        rows.append({"t": i, "mean_s(t)_sum": mean_st_sum})
        df = pd.DataFrame(rows)
    """   

    
    
