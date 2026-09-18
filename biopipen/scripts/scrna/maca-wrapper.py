"""Run MACA on an AnnData and save the per-cell labels.

Used by CellTypeAnnotation-maca.R (through
biopipen.utils::RunCellTypeAnnotation()), the same way as the other
python-based annotation tools.

Prerequisites
-------------
MACA is on PyPI (`MACA-Python`), but that release pins `scanpy==1.6.0` and
`anndata==0.7.5`, which cannot coexist with a modern scanpy. Use the
modernized fork (it wants `scanpy>=1.10`/`anndata>=0.10`) builtin below
using the python that runs this wrapper (the `envs.maca.python` of the
process). MACA is imported lazily and a clear error is raised naming that
environment when the import fails.

`maca.singleMACA()` scores the cells against the marker sets (`cell_markers`,
a dict of cell type to marker genes, built here from the 2-column marker file),
clusters the score matrix at every combination of `res` and `n_neis`, and maps
the clusters to cell types by the majority of their cells. The per-cell labels
are aligned with the rows of the AnnData.

MACA keeps only the markers that are in the AnnData and then drops every cell
type left with fewer than 3 (or more than 300) of them, so a marker table whose
genes are largely absent from the object -- or that lists too many of them for
a cell type -- leaves nothing to annotate and MACA's own `np.argmax` fails on an
empty sequence. That case is reported with the per-cell-type marker counts
instead of the raw ValueError.
"""
from argparse import ArgumentParser
import sys

"""
Created on Tue Jul 21 16:31:17 2020

@author: Yang Xu
Patched by pwwang 2026
"""
import umap
import scipy
#import random
import anndata
import collections
import numpy as np
import pandas as pd
import scanpy as sc
import multiprocessing

from sklearn.decomposition import PCA
from sklearn.metrics import confusion_matrix
from sklearn.feature_extraction.text import TfidfTransformer

import warnings
warnings.filterwarnings("ignore")

##-----------------------------------------------------------------------------
def ensemble_labels(multi_labels=None):
    ensemble = []
    for i in range(multi_labels.shape[0]):
        ks=[]
        vs=[]
        for k,v in collections.Counter(multi_labels[i,:].tolist()).items():
            ks.append(k)
            vs.append(v)
        ensemble.append(ks[vs.index(max(vs))])
    return ensemble

def singleMACA(ad=None, cell_markers=None,n_pcs=None,res=[1,2,3],n_neis = [5,10],
               freq=0.5,use_weight=False):
    ##TF-IDF transformation
    X = ad.X.copy()
    if scipy.sparse.issparse(X):
        X = X.toarray()

    tf_transformer = TfidfTransformer(use_idf=True).fit(X)
    X= tf_transformer.transform(X).toarray()

    labels = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    labels.columns = cell_markers.keys()
    exprsed = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    exprsed.columns = cell_markers.keys()
    celltype_size = {}

    ##create artifical labels for each cell
    if use_weight == True:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            marker_index = -1
            for i in v:
                marker_index += 1
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0

                    ##consider marker weight
                    l = l*(1-marker_index/(len(v)*2))##default 2

                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])

    else:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            for i in v:
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])

    assess1 = np.argmax((labels*exprsed).values,axis=1)
    vals1 = 0
    for k,v in collections.Counter(assess1).items():
        if v >= 5:
            vals1 += 1

    assess1 = vals1

    assess2 = np.argmax((labels).values,axis=1)
    vals2 = 0
    for k,v in collections.Counter(assess2).items():
        if v >= 5:
            vals2 += 1

    assess2 = vals2

    ass = [assess1,assess2]

    new_labels = [labels*exprsed,labels][ass.index(max(ass))]
    #new_labels = labels*exprsed##consider the number of expressed marker of each cell-type for each cell
    #new_labels = labels#*exprsed

    celltype_size = pd.DataFrame.from_dict(celltype_size,orient='index')

    ##remove cell types with > 300 and < 5 marker genes
    labels = labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    new_labels = new_labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    celltype_size = celltype_size.iloc[np.logical_and(celltype_size.values<=300,celltype_size.values>=3),:]
    print(labels.shape)

    ##create Label1
    labels1 = np.argmax(new_labels.values,axis=1)
    ad.obs['Label1'] = labels1

    ##UMAP visualization
    embedding = umap.UMAP(n_neighbors=15,min_dist=0.2,n_components=2,
                          metric='cosine').fit_transform(labels.values)#(l2)#
    embedding = pd.DataFrame(embedding)
    embedding.columns=['UMAP1','UMAP2']
    ad.obsm['X_umap'] = embedding.iloc[:,:2].values

    ##create Label2
    if n_pcs is not None:
        pca = PCA(n_components=n_pcs)
        scores = pca.fit_transform(labels.values)
        ad.obsm['Score']=scores
    else:
        ad.obsm['Score']=labels.values

    ##scanpy reuses a pre-existing adata.uns['neighbors'] instead of building the
    ##graph from use_rep; a foreign one (e.g. carried over from a Seurat
    ##conversion) has no params['n_neighbors'] and would raise a KeyError here
    ad.uns.pop('neighbors', None)

    label_list = np.zeros((ad.X.shape[0],len(res)*len(n_neis))).astype('str')
    indexs = 0
    for r in res:
        for nei in n_neis:

            sc.pp.neighbors(ad, use_rep="Score", n_neighbors=nei,metric='cosine')
            ##sc.tl.louvain was removed from scanpy; leiden with the igraph flavor,
            ##2 iterations and an undirected graph is scanpy's Louvain-equivalent
            sc.tl.leiden(ad, resolution=r, key_added = 'louvain',
                         flavor='igraph', n_iterations=2, directed=False)

            ##Map Label2 to Label1; remove all non-candiate cell types
            labels1_2 = list(collections.Counter(labels1.tolist()).keys())
            labels_2 = labels.copy().iloc[:,labels1_2]
            new_labels_2 = new_labels.copy().iloc[:,labels1_2]
            #celltype_size_2 = celltype_size.iloc[labels1_2,:].copy()
            celltype_names = new_labels_2.columns

            cm = confusion_matrix(ad.obs['louvain'].values.astype(int),
                                  np.argmax(new_labels_2.values,axis=1).astype(int))
            cm = cm[:np.max(ad.obs['louvain'].values.astype(int)),:]

            normed_cm = cm.copy().T
            normed_cm = normed_cm/np.sum(normed_cm,axis=0)
            normed_cm = np.nan_to_num(normed_cm)
            normed_cm = normed_cm.T
            mapping={}

            mapping = np.argmax(normed_cm,axis=1)
            mapmax = np.max(normed_cm,axis=1)
            mapmax = np.nan_to_num(mapmax)

            clustering = ad.obs['louvain'].values.astype(int)
            new_cluster = np.zeros((len(clustering)))
            for i in range(len(mapping)):
                tof = clustering==i
                if mapmax[i]>=freq:
                    new_cluster[tof]=mapping[i]
                else:
                    sub_labels = new_labels_2.values[tof,:]
                    if sub_labels.shape[0]>0:
                        freqs = []
                        for j in range(sub_labels.shape[0]):
                            zscore=scipy.stats.zscore(sub_labels[j,:])
                            orderi = np.array([g for g in range(sub_labels.shape[1])])[zscore>3].tolist()
                            orderj = np.argsort(sub_labels[j,:])[::-1][:3].tolist()#default 3
                            a = [len(orderi),len(orderj)]
                            freqs += [orderi,orderj][a.index(max(a))]
                        vals = 0
                        ks = 0
                        for k,v in collections.Counter(freqs).items():
                            if v >= vals:
                                ks = k
                                vals = v

                        if vals/sub_labels.shape[0]>=freq:
                            new_cluster[tof]=ks
                        else:
                            new_cluster[tof]=-1

            mapped =new_cluster.astype('int')
            mapped=mapped.astype('str')
            ad.obs['Mapped'] = mapped

            cell_dict = {}
            for k,v in collections.Counter(new_cluster.tolist()).items():
                if int(k)>=0:
                    cell_dict[int(k)]=labels_2.columns[int(k)]
                else:
                    cell_dict[int(k)]="unassigned"
            label_list[:,indexs]=ad.obs['Mapped'].values
            indexs+=1

    ensemble = ensemble_labels(label_list)
    ensemble = np.array(ensemble)

    annotations=[]
    for e in ensemble:
        if int(e)>=0:
            annotations.append(celltype_names[int(e)])
        else:
            annotations.append("unassigned")

    return ad, np.array(annotations)

def gene2cell(ad=None, cell_markers=None,use_weight=False):
    ##TF-IDF transformation
    X = ad.X.copy()
    if scipy.sparse.issparse(X):
        X = X.toarray()

    tf_transformer = TfidfTransformer(use_idf=True).fit(X)
    X= tf_transformer.transform(X).toarray()

    labels = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    labels.columns = cell_markers.keys()
    exprsed = pd.DataFrame(np.zeros((X.shape[0],len(cell_markers))))
    exprsed.columns = cell_markers.keys()
    celltype_size = {}

    ##create artifical labels for each cell
    if use_weight == True:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            marker_index = -1
            for i in v:
                marker_index += 1
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0

                    ##consider marker weight
                    l = l*(1-marker_index/(len(v)*2))##default 2

                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])

    else:
        for k, v in cell_markers.items():
            celltype_size[k]=0
            sums=0
            n = np.zeros((X.shape[0]))
            for i in v:
                if i in ad.var.index:
                    expr95 = np.percentile(X[:,ad.var.index == i],95)
                    thresh = 0.25 * expr95
                    l = np.array(X[:,ad.var.index == i])
                    l[X[:,ad.var.index == i]<=thresh]=0
                    n[np.array(l>0).reshape(X.shape[0])] += 1
                    sums += 1
                    labels[k] += l.reshape(X.shape[0])
            n = n/sums
            celltype_size[k]=sums
            exprsed[k] = n.reshape(X.shape[0])

    assess1 = np.argmax((labels*exprsed).values,axis=1)
    vals1 = 0
    for k,v in collections.Counter(assess1).items():
        if v >= 5:
            vals1 += 1

    assess1 = vals1

    assess2 = np.argmax((labels).values,axis=1)
    vals2 = 0
    for k,v in collections.Counter(assess2).items():
        if v >= 5:
            vals2 += 1

    assess2 = vals2

    ass = [assess1,assess2]

    new_labels = [labels*exprsed,labels][ass.index(max(ass))]
    #new_labels = labels*exprsed##consider the number of expressed marker of each cell-type for each cell
    #new_labels = labels#*exprsed

    celltype_size = pd.DataFrame.from_dict(celltype_size,orient='index')

    ##remove cell types with > 300 and < 5 marker genes
    labels = labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    new_labels = new_labels.iloc[:,np.logical_and(celltype_size.values<=300,celltype_size.values>=3)]
    celltype_size = celltype_size.iloc[np.logical_and(celltype_size.values<=300,celltype_size.values>=3),:]
    print(labels.shape)

    return labels, new_labels

def multiMACA(labels=None,new_labels=None):

    ad = anndata.AnnData(X=labels)
    ##create Label1
    labels1 = np.argmax(new_labels.values,axis=1)
    ad.obs['Label1'] = labels1

    ##create Label2
    #if n_pcs is not None:
    #    pca = PCA(n_components=n_pcs)
    #    scores = pca.fit_transform(labels.values)
    #    ad.obsm['Score']=scores
    #else:
    #    ad.obsm['Score']=labels.values

    ##l2 normalize
    #transformer = Normalizer().fit(new_labels)
    #l2 = transformer.transform(new_labels)

    ad.obsm['Score']=labels.values#l2#
    sc.pp.neighbors(ad, use_rep="Score", n_neighbors=5,metric='cosine')
    sc.tl.leiden(ad, resolution=2, key_added = 'louvain',
                 flavor='igraph', n_iterations=2, directed=False)##default 1

    ##Map Label2 to Label1; remove all non-candiate cell types
    labels1_2 = list(collections.Counter(labels1.tolist()).keys())
    #labels_2 = labels.copy().iloc[:,labels1_2]
    new_labels_2 = new_labels.copy().iloc[:,labels1_2]
    #celltype_size_2 = celltype_size.iloc[labels1_2,:].copy()
    celltype_names = new_labels_2.columns

    cm = confusion_matrix(ad.obs['louvain'].values.astype(int),
                          np.argmax(new_labels_2.values,axis=1).astype(int))
    cm = cm[:np.max(ad.obs['louvain'].values.astype(int)),:]

    normed_cm = cm.copy().T
    normed_cm = normed_cm/np.sum(normed_cm,axis=0)
    normed_cm = np.nan_to_num(normed_cm)
    normed_cm = normed_cm.T
    mapping={}

    mapping = np.argmax(normed_cm,axis=1)
    mapmax = np.max(normed_cm,axis=1)
    mapmax = np.nan_to_num(mapmax)

    clustering = ad.obs['louvain'].values.astype(int)
    new_cluster = np.zeros((len(clustering)))
    for i in range(len(mapping)):
        tof = clustering==i
        if mapmax[i]>=0.5:
            new_cluster[tof]=mapping[i]
        else:
            sub_labels = new_labels_2.values[tof,:]
            if sub_labels.shape[0]>0:
                freqs = []
                for j in range(sub_labels.shape[0]):
                    zscore=scipy.stats.zscore(sub_labels[j,:])
                    orderi = np.array([g for g in range(sub_labels.shape[1])])[zscore>3].tolist()
                    orderj = np.argsort(sub_labels[j,:])[::-1][:3].tolist()##default 3
                    a = [len(orderi),len(orderj)]
                    freqs += [orderi,orderj][a.index(max(a))]
                    vals = 0
                    ks = 0
                    for k,v in collections.Counter(freqs).items():
                        if v >= vals:
                            ks = k
                            vals = v

                    if vals/sub_labels.shape[0]>=0.5:
                        new_cluster[tof]=ks
                    else:
                        new_cluster[tof]=-1

        mapped=[]
        for e in new_cluster:
            if int(e)>=0:
                mapped.append(celltype_names[int(e)])
            else:
                mapped.append("unassigned")

    return np.array(mapped)

def parallel(scores=None,labels=None,batch_size=20000,repeats=9,n_core=10):
    index = np.array([i for i in range(scores.shape[0])])
    label_list = np.zeros((scores.shape[0],repeats)).astype('str')

    for i in range(repeats):
        r = np.random.permutation(scores.shape[0])
        r_index = index[r]
        r_scores = scores.iloc[r,:]
        r_label1 = labels.iloc[r,:]
        scores_list= []
        for j in range(scores.shape[0]//batch_size+1):
            scores_list.append((r_scores.iloc[j*batch_size:(j+1)*batch_size,:],
                                r_label1.iloc[j*batch_size:(j+1)*batch_size,:]))

        pool = multiprocessing.Pool(processes=n_core)
        mapped = pool.starmap(multiMACA, scores_list)
        merged = []
        for m in mapped:
            merged += m.tolist()
        merged = np.array(merged)
        merged = merged[r_index.argsort()]
        label_list[:,i]=merged

    return label_list



def csv_list(value, cast):
    """Parse a comma/space separated list, e.g. `1,2,3`."""
    if value is None:
        return None
    return [cast(item) for item in value.replace(",", " ").split()]


def main():
    parser = ArgumentParser(description="Run MACA")
    parser.add_argument(
        "-i", "--input", required=True, help="Input H5AD file (AnnData)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file (barcode, cell type)"
    )
    parser.add_argument(
        "-m", "--marker", required=True,
        help="Marker file: a 2-column headerless TSV of cell type and gene"
    )
    parser.add_argument(
        "--n-pcs", type=int, default=None,
        help="Number of principal components of the marker scores"
    )
    parser.add_argument(
        "--res", default=None,
        help="Comma separated Louvain resolutions (default: MACA's [1, 2, 3])"
    )
    parser.add_argument(
        "--n-neis", default=None,
        help="Comma separated numbers of neighbors (default: MACA's [5, 10])"
    )
    parser.add_argument(
        "--freq", type=float, default=None,
        help="Frequency threshold of the cluster mapping (default: 0.5)"
    )
    parser.add_argument(
        "--use-weight", action="store_true",
        help="Weight the markers by their order in the marker file"
    )
    args = parser.parse_args()

    adata = sc.read_h5ad(args.input)
    markers = pd.read_csv(
        args.marker, sep="\t", header=None, names=["cell_type", "gene"]
    )
    cell_markers = {
        str(cell_type): genes.astype(str).tolist()
        for cell_type, genes in markers.groupby("cell_type")["gene"]
    }
    print(
        f"MACA: annotating {adata.n_obs} cells with "
        f"{len(cell_markers)} cell type(s)"
    )

    # None is MACA's own default for these: leave them out so the defaults stay
    # in one place (and `res`/`n_neis` must be lists)
    kwargs = {
        "n_pcs": args.n_pcs,
        "res": csv_list(args.res, float),
        "n_neis": csv_list(args.n_neis, int),
        "freq": args.freq,
        "use_weight": args.use_weight,
    }
    kwargs = {key: value for key, value in kwargs.items() if value is not None}

    try:
        _, labels = singleMACA(
            ad=adata, cell_markers=cell_markers, **kwargs
        )
    except ValueError as exc:
        if "argmax" not in str(exc):
            raise
        # `celltype_size` counts the markers of each cell type that MACA found
        # among the object's features; with every cell type outside MACA's
        # 3..300 range, `labels` is narrowed to 0 columns and the argmax above
        # has nothing to pick from
        counts = "\n".join(
            "  {}: {} of {} marker(s) in the object".format(
                cell_type,
                sum(gene in adata.var_names for gene in genes),
                len(genes),
            )
            for cell_type, genes in cell_markers.items()
        )
        sys.exit(
            f"MACA cannot annotate anything ({exc}). It keeps only the "
            f"markers that are in the {adata.n_vars} features of the h5ad and "
            "drops every cell type left with fewer than 3 or more than 300 "
            "of them:\n"
            f"{counts}\n"
            "Check the marker table's genes against the object's features "
            "(a Seurat object is converted with only its variable features), "
            "and install the modernized fork with:\n"
            "  pip install -e ~/github/MACA"
        )
    pd.DataFrame(
        {"cell": adata.obs_names, "maca_celltype": labels}
    ).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
