import numpy as np
import matplotlib.pyplot as plt

from astropy.table import Table
from sklearn.decomposition import PCA

#=====================================

orbs = Table.read('../queorbita_new.out', data_start=50, header_start=49, format='ascii')

orbs_matrix = np.vstack(np.array(orbs[str(cols)]) for cols in orbs.columns).T
orbs_matrix_just_floats = np.delete(orbs_matrix, [21, 23], 1)

pca = PCA(n_components=5)
pca.fit(orbs_matrix_just_floats)

#=====================================
# 
# def draw_vector(v0, v1, ax=None):
#     ax = ax or plt.gca()
#     arrowprops=dict(arrowstyle='->',
#                     linewidth=2,
#                     shrinkA=0, shrinkB=0)
#     ax.annotate('', v1, v0, arrowprops=arrowprops)
#
# # plot data
# plt.scatter(orbs_matrix_just_floats[:, 0], orbs_matrix_just_floats[:, 1], alpha=0.2)
# for length, vector in zip(pca.explained_variance_[0:2], pca.components_):
#     v = vector * 3 * np.sqrt(length)
#     draw_vector(pca.mean_[0:2], pca.mean_[0:2] + v)
# plt.axis('equal');
