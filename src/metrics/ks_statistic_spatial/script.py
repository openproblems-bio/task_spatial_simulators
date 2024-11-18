import anndata as ad
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
import squidpy as sq
from scipy.stats import ks_2samp

## VIASH START
par = {
  'input_spatial_dataset': 'temp_prostate.srtsim.ks_statistic_spatial/_viash_par/input_spatial_dataset_1/output_sp.h5ad',
  # 'input_singlecell_dataset': 'resources_test/spatialsimbench_mobnew/MOBNEW_sc.rds',
  'input_simulated_dataset': 'temp_prostate.srtsim.ks_statistic_spatial/_viash_par/input_simulated_dataset_1/gastrulation.srtsim.generate_sim_spatialcluster.output_sp.h5ad',
  'output': 'output.h5ad'
}
meta = {
  'name': 'ks_statistic_spatial'
}

## VIASH END

print('Reading input files', flush=True)
input_spatial_dataset = ad.read_h5ad(par['input_spatial_dataset'])
# input_singlecell_dataset = ad.read_h5ad(par['input_singlecell_dataset'])
input_simulated_dataset = ad.read_h5ad(par['input_simulated_dataset'])

def get_spatial_network(num_sample=None, spatial=None, radius=None, coord_type="grid", n_rings=2, set_diag=False):
    spatial_adata = ad.AnnData(np.empty((num_sample, 1), dtype="float32"))
    spatial_adata.obsm["spatial"] = spatial
    # sq.gr.spatial_neighbors(spatial_adata, n_rings=n_rings, coord_type=coord_type, n_neighs=n_neighs, radius=radius,set_diag =set_diag)
    sq.gr.spatial_neighbors(spatial_adata, n_rings=n_rings, coord_type=coord_type, radius=radius, set_diag=set_diag,
                            delaunay=True)
    sn = spatial_adata.obsp["spatial_connectivities"]

    return sn


def get_onehot_ct(init_assign=None):
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(init_assign)
    onehot_encoder = OneHotEncoder(sparse_output=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_ct = onehot_encoder.fit_transform(integer_encoded)
    return onehot_ct.astype(np.float32)


# @numba.jit("float32[:, ::1](float32[:, ::1], float32[:, ::1])")
def get_nb_freq(nb_count=None, onehot_ct=None):
    #     nb_freq = onehot_ct.T @ nb_count
    nb_freq = np.dot(onehot_ct.T, nb_count)
    res = nb_freq / nb_freq.sum(axis=1).reshape(onehot_ct.shape[1], -1)
    return res

def get_trans(adata=None, ct=None):
    sn = get_spatial_network(num_sample=adata.obs.shape[0],
                             spatial=adata.obsm["spatial"], coord_type="generic")
    onehot_ct = get_onehot_ct(init_assign=ct)
    nb_count = np.array(sn * onehot_ct, dtype=np.float32)
    target_trans = get_nb_freq(nb_count=nb_count, onehot_ct=onehot_ct)
    return target_trans


input_spatial_dataset.obsm["spatial"] = np.array(input_spatial_dataset.obs[['col', 'row']].values.tolist())
input_spatial_dataset.obs["celltype"] = input_spatial_dataset.obs["spatial_cluster"]
input_spatial_dataset.obs["celltype"] = input_spatial_dataset.obs["celltype"].astype('category')

sq.gr.spatial_neighbors(input_spatial_dataset, coord_type="generic", set_diag=False, delaunay=True)
# neighborhood enrichment matrix
sq.gr.nhood_enrichment(input_spatial_dataset, cluster_key="celltype")
# centrality scores matrix
sq.gr.centrality_scores(input_spatial_dataset, cluster_key="celltype")

input_simulated_dataset.obsm["spatial"] = np.array(input_simulated_dataset.obs[['col', 'row']].values.tolist())
input_simulated_dataset.obs["celltype"] = input_simulated_dataset.obs["spatial_cluster"]
input_simulated_dataset.obs["celltype"] = input_simulated_dataset.obs["celltype"].astype('category')

# remove NaN before compute the matrix
valid_indices = input_simulated_dataset.obs["celltype"].notnull()
input_simulated_dataset = input_simulated_dataset[valid_indices]
input_spatial_dataset = input_spatial_dataset[valid_indices]

sq.gr.spatial_neighbors(input_simulated_dataset, coord_type="generic", set_diag=False, delaunay=True)
# centrality scores matrix
sq.gr.centrality_scores(input_simulated_dataset, cluster_key="celltype")
# neighborhood enrichment matrix
input_simulated_dataset = input_simulated_dataset[input_simulated_dataset.obs["celltype"].notnull()]
sq.gr.nhood_enrichment(input_simulated_dataset, cluster_key="celltype")


target_enrich_real = input_spatial_dataset.uns["celltype_nhood_enrichment"]["zscore"]
target_enrich_scale_real = target_enrich_real/np.max(target_enrich_real)
target_enrich_sim = input_simulated_dataset.uns["celltype_nhood_enrichment"]["zscore"]
target_enrich_scale_sim = target_enrich_sim/np.max(target_enrich_sim)

#error_enrich = np.linalg.norm(target_enrich_sim - target_enrich_real)
#error_enrich_scale = np.linalg.norm(target_enrich_scale_sim - target_enrich_scale_real)
    
target_enrich_real_ds = target_enrich_real.flatten()
target_enrich_sim_ds = target_enrich_sim.flatten()
ks_enrich, p_value = ks_2samp(target_enrich_real_ds, target_enrich_sim_ds)

# KS central

real_central_real = np.array(input_spatial_dataset.uns["celltype_centrality_scores"])
real_central_sim = np.array(input_simulated_dataset.uns["celltype_centrality_scores"])
    
real_central_real_ds = real_central_real.flatten()
real_central_sim_ds = real_central_sim.flatten()
ks_central, p_value = ks_2samp(real_central_real_ds, real_central_sim_ds)

# transition matrix
real = np.array(input_spatial_dataset.obs['spatial_cluster'].values.tolist())
sim = np.array(input_simulated_dataset.obs['spatial_cluster'].values.tolist())


transition_matrix_real = get_trans(adata=input_spatial_dataset, ct=real)
transition_matrix_sim = get_trans(adata=input_simulated_dataset, ct=sim)

# error = np.linalg.norm(transition_matrix_sim - transition_matrix_real)
transition_matrix_real_ds = transition_matrix_real.flatten()
transition_matrix_sim_ds = transition_matrix_sim.flatten()
ks_stat_error, p_value = ks_2samp(transition_matrix_real_ds, transition_matrix_sim_ds)


uns_metric_ids = [
  "ks_statistic_transition_matrix",
  "ks_statistic_central_score",
  "ks_statistic_enrichment"
]

uns_metric_values = [
  ks_stat_error,
  ks_central,
  ks_enrich
]

print("Write output AnnData to file", flush=True)
output = ad.AnnData(
  uns={
    'dataset_id': input_simulated_dataset.uns['dataset_id'],
    'method_id': input_simulated_dataset.uns['method_id'],
    'metric_ids': uns_metric_ids,
    'metric_values': uns_metric_values
  }
)
output.write_h5ad(par['output'], compression='gzip')
