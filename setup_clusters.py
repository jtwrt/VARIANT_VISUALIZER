import sys
from variant_visualizer import clusters
from pathos.multiprocessing import ProcessingPool

def main(process_number):
    """Run instance of ClusterGenerator"""
    clusters.ClusterGenerator(
        process_number=process_number,
        n_processes=n_processes)

if __name__ == '__main__':
    global n_processes
    n_processes = int(sys.argv[1])
    print(f'Generating Clusters using {n_processes} parallel processes ...')
    with ProcessingPool(n_processes) as pool:
        pool.map(main, [i for i in range(0,n_processes)])
    print('Finished generating clusters!')