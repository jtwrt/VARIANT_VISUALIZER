import argparse
from variant_visualizer import clusters
from pathos.multiprocessing import ProcessingPool

def main(process_number):
    """Run instance of ClusterGenerator"""
    clusters.ClusterGenerator(
        process_number=process_number,
        n_processes=n_processes,
        selected_clusters=selected_clusters)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-n','--n_processes',
                        dest='n_processes',
                        help='<Required> Number of parallel processes that will be launched.', 
                        required=True,
                        type=int)
    parser.add_argument('-c','--cluster_ids', 
                        dest='selected_clusters',
                        nargs='+', # parse as list
                        help='<Optional> Clusters that will be generated. Generates all clusters if not provided.',
                        required=False,
                        type=int)
    args = parser.parse_args()

    global n_processes
    n_processes = args.n_processes
    global selected_clusters
    selected_clusters = args.selected_clusters

    print(f'Generating Clusters using {n_processes} parallel processes ...')
    with ProcessingPool(n_processes) as pool:
        pool.map(main, [i for i in range(0,n_processes)])
    print('Finished generating clusters!')