import argparse
from data_io import read_and_process_AML, read_from_folder, read_and_process_BMMC
import phenograph
from sklearn.neighbors import kneighbors_graph
import networkx as nx
import community
import markov_clustering as mc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import Normalize
from sklearn import metrics
import os
from sklearn.manifold import t_sne
import time
import DBCV
from scipy.spatial.distance import euclidean

def build_parser():
    """
    Make object that parses command like arguments
    :return: argparse.ArgumentParser object
    """
    parser = argparse.ArgumentParser(description=
                                     'Run Louvain community detection, Markov clustering or phenograph on dataset.'
                                     'Picture of clusters as embedding is saved at results/tsne_mode_k.png',
                                     prog="main.py")
    parser.add_argument('filepath', type=str,
                        help='Filepath of data, expect csv or folder with csvs seperated in labels and samples')
    parser.add_argument('mode', type=str, help='Which algorithm to run, "lou", "mcl" or "pheno"',
                        choices={"lou", "mcl", "pheno"})
    parser.add_argument('k', type=int, help='Number of neighbours for knn graph')
    parser.add_argument('--f', type=str, help='File storage type, either a single csv or '
                                              'a folder with samples and labels subfolders containing csvs',
                        default="AML", nargs='?', choices={"AML", "PANORAMA", "BMMC"})
    parser.add_argument('--e', type=str, help='Location of t-sne embeddings if precomputed',
                        default="embeddings.csv", nargs='?')
    parser.add_argument('--i', type=float,
                        help='Inflation parameter for Markov clustering', nargs='?', default=None)
    return parser


def read_data(filepath, filetype):
    """

    :param filepath: string pointing to file
    :param filetype: string specifying which file reader function to call
    :return: data, dataframe with column for each protein, datalabels, dataframe with cell_type column
    """
    if filetype == 'AML':
        data, datalabels = read_and_process_AML(filepath)
        data = data.reset_index().drop(columns='index')
        datalabels = datalabels.reset_index().drop(columns='index')
    elif filetype == 'PANORAMA':
        data, datalabels = read_from_folder(filepath)
    elif filetype == 'BMMC':
        data, datalabels = read_and_process_BMMC(filepath)
    return data, datalabels


def run_louvain(data, k):
    """
    Run louvain clustering on data
    :param data: dataframe
    :param k: int, k in knn
    :return: list, cluster assignments
    """
    sparse_graph = build_sparse_matrix(data, k)
    graph = nx.Graph(sparse_graph)
    partition = community.best_partition(graph, randomize=False)
    return list(partition.values())


def build_sparse_matrix(data, k, return_distance=True):
    """
    :param data: dataframe
    :param k: k in knn
    :param return_distance: return similarity if False, distance if True
    :return: scipy sparse matrix in coo format
    """
    knn = kneighbors_graph(data.values, n_neighbors=k, mode='distance')
    if not return_distance:
        knn.data = 1 / knn.data
    return knn.tocoo()


def run_mcl(data, k, inflation):
    """
    Markov clustering on data
    :param data: dataframe
    :param k: int, k in knn
    :param inflation: float, algorithm parameter
    :return: list of cluster assignments
    """
    knn = build_sparse_matrix(data, k)
    graph = nx.to_scipy_sparse_matrix(nx.Graph(knn), format='coo')
    result = mc.run_mcl(graph, inflation=inflation)
    clusters = mc.get_clusters(result)
    clustersize = sum([len(cluster) for cluster in clusters])
    clusters_flat = np.zeros(clustersize)
    for cluster_num, cluster in enumerate(clusters):
        for cell in cluster:
            clusters_flat[cell] = cluster_num
    return clusters_flat


def run_tsne(args, data):
    embeddings = t_sne.TSNE(n_components=2).fit_transform(data.values)
    np.savetxt(args.f + "embeddings.csv", embeddings, delimiter=",")
    return embeddings

def read_embeddings(args, data):
    '''Read TSNE embeddings from args.e
    :param data: dataframe
    :param args:argparser object
    :return: flat float64 numpy arrays x and y
    '''
    try:
        embeddings = pd.read_csv(args.f + args.e, header=None)

    except FileNotFoundError:
        embeddings = run_tsne(args, data)
        embeddings = embeddings.transpose()

    x = embeddings[0]
    y = embeddings[1]
    return x, y


def cluster_to_label(clusters, labels):
    '''Count occurence of label per cluster to assign each cluster a label by maximum occurence
    :param clusters: list, cluster assignments
    :param labels: dataframe, ground truth
    :return: list of most prevalent label per cluster
    '''
    louvaindf = pd.DataFrame(clusters)
    louvaindf.columns = ['community']
    grouped = louvaindf.groupby('community').groups
    clusterlabels = []
    for group in grouped:
        labs, counts = np.unique(labels.iloc[grouped[group]], return_counts=True)
        maximum = max(counts)
        index = list(counts).index(maximum)
        clusterlabels.append(labs[index])
    return clusterlabels


def draw_embeddings(clusters, data, labels, args):
    '''Saves TSNE embeddings colored with clustering results to ./results/
    Don't look at this horrible code :@
    :param clusters: list of cluster assignments
    :param data: dataframe
    :param labels:dataframe
    :param args: argparser object
    '''
    labelset = list(set(list(clusters)))
    x, y = read_embeddings(args, data)

    x = x[:len(clusters)]
    y = y[:len(clusters)]
    print("len x %i" % (len(x)))
    print("len y %i" % (len(y)))
    print("Len clusters %i" % (len(clusters)))
    norm = Normalize(vmin=0, vmax=len(labelset) - 1)
    cmap = plt.cm.get_cmap('nipy_spectral')
    map_colors = []
    for i in range(len(labelset)):
        map_colors.append(cmap(norm(i)))
    # Gather clusternumber's color for each cluster assignment
    plot_colors = [map_colors[labelset.index(cell)] for cell in clusters]
    # Plot TSNE embeddings with clusternumber as color
    plt.scatter(x, y, c=plot_colors, s=2, alpha=0.5)
    plt.xlabel('TSNE 1')
    plt.ylabel('TSNE 2')
    # Get cluster counts for legend
    _, counts = np.unique(clusters, return_counts=True)
    clusterlabels = cluster_to_label(clusters, labels)
    legend_elements = [Patch(label=str(ctype) + ': ' + str(counts[labelset.index(ctype)]) +
                                   ' ' + str(clusterlabels[labelset.index(ctype)])
                             , facecolor=map_colors[int(labelset.index(ctype))]) for ctype in labelset]
    plt.legend(handles=legend_elements, prop={'size': 4}, loc='lower right')
    plt.title('Clustering in TSNE space for k is %i'
              ' \nNumber of clusters is %i'
              '\nSilhouette score is %f'
              % (args.k, len(set(clusters)), metrics.silhouette_score(data.values, clusters)))
    plt.savefig('./results/tsne_' + args.mode + '_' + args.f + str(args.k) + '.png', format='png', dpi=300)


def gather_statistics(data, datalabels, clusters, args):
    """
    :param data: dataframe
    :param datalabels: dataframe
    :param clusters: list
    :param args: argparser object
    :return: dictionary
    """
    dfdict = {"K": [args.k],
              "Clusters": [len(set(clusters))],
              'Completeness': [calc_completeness(datalabels, clusters)],
              'Homogeneity': [calc_homogeneity(datalabels, clusters)],
              'Adjusted Mutual Information Score': [calc_MIS(datalabels, clusters)],
              'Adjusted Rand Index': [calc_rand(datalabels, clusters)],
              'Silhouette': [calc_silhouette(clusters, data)],
              'Davies Bouldain': [calc_davies(data, clusters)],
              'Calinski Harabaz': [calc_calhar(data, clusters)],
              'mode': [args.mode],
              'inflation': [args.i],
              'dataset': [args.f]
#              'Density based evaluation' : [calc_DBCV(data, clusters)]
              }
    statistics = pd.DataFrame(dfdict)
    return statistics

def calc_DBCV(data, clusters):
    tic = time.time()
    dbscore = DBCV.DBCV(data.values, clusters, dist_function=euclidean)
    toc = time.time()
    print(toc-tic)
    return dbscore

def calc_completeness(labels, clusters):
    return metrics.completeness_score(labels.values.ravel(), clusters)


def calc_homogeneity(labels, clusters):
    return metrics.homogeneity_score(labels.values.ravel(), clusters)


def calc_MIS(labels, clusters):
    return metrics.adjusted_mutual_info_score(labels.values.ravel(), clusters)


def calc_rand(labels, clusters):
    return metrics.adjusted_rand_score(labels.values.ravel(), clusters)


def calc_silhouette(clusters, data):
    return metrics.silhouette_score(data.values, clusters, metric='euclidean')


def calc_davies(data, clusters):
    return metrics.davies_bouldin_score(data.values, clusters)


def calc_calhar(data, clusters):
    return metrics.calinski_harabaz_score(data.values, clusters)


def main():
    # Parse command line arguments
    parser = build_parser()
    args = parser.parse_args()
    k = args.k
    filepath = args.filepath
    mode = args.mode.lower()
    filetype = args.f
    inflation = args.i
    # Inflation is only used as a parameter for MCL
    if hasattr(args, 'inflation') and args.mode != 'mcl':
        parser.error("Inflation is only a valid argument in combination with mcl algorithm")
    # Dataframes for data and datalabels
    data, datalabels = read_data(filepath, filetype)
    # Start specified method
    if mode == 'pheno':
        clusters, _, _ = phenograph.cluster(data.values, n_jobs=1, k=k)
    elif mode == 'lou':
        clusters = run_louvain(data, k)
    elif mode == 'mcl':
        clusters = run_mcl(data, k, inflation)
    # Draw clusters on embeddings
    #draw_embeddings(clusters, data, datalabels, args)
    # Save stats to new file
    stats_df = gather_statistics(data, datalabels, clusters, args)
    if not os.path.isfile('./results/benchmark_results.csv'):
        stats_df.to_csv('./results/benchmark_results.csv', header='column_names')
    else:  # else it exists so append without writing the header
        stats_df.to_csv('./results/benchmark_results.csv', mode='a', header=False)


if __name__ == '__main__':
    main()
