# phenograph_analysis
Comparison of Phenograph clustering to native Louvain method on same graph
# Benchmark
Comparison of Phenograph vs Louvain on BMMC AML and PANORAMA CyTOF benchmark sets, Markov clustering is an option but too slow
for datasets of this size.


    usage: main.py [-h] [--f [{BMMC,AML,PANORAMA}]] [--e [E]] [--i [I]]
                   filepath {pheno,lou,mcl} k
    Run Louvain community detection, Markov clustering or phenograph on
    dataset.Picture of clusters as embedding is saved at results/tsne_mode_k.png

    positional arguments:
      filepath              Filepath of data, expect csv or folder with csvs
                            seperated in labels and samples
      {lou,mcl,pheno}       Which algorithm to run, "lou", "mcl" or "pheno"
      k                     Number of neighbours for knn graph

    optional arguments:
      -h, --help            show this help message and exit
      --f [{AML,PANORAMA,BMMC}]
                            File storage type, either a single csv or a folder
                            with samples and labels subfolders containing csvs
      --e [E]               Location of t-sne embeddings if precomputed
      --i [I]               Inflation parameter for Markov clustering
