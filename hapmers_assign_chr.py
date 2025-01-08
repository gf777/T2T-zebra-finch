import itertools
import numpy as np
import pandas as pd
import sys
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import scipy.stats as stats
from scipy import linalg

# Usage: eg. python hapmers_assign_chr.py hapmers.out
color_iter = itertools.cycle(["navy", "c", "cornflowerblue"])

def plot_results(X, Y_, means, covariances, index, title):
    splot = plt.subplot()
    for i, (mean, covar, color) in enumerate(zip(means, covariances, color_iter)):
        v, w = linalg.eigh(covar)
        v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
        u = w[0] / linalg.norm(w[0])
        # as the DP will not use every component it has access to
        # unless it needs it, we shouldn't plot the redundant
        # components.
        if not np.any(Y_ == i):
            continue
        plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 0.8, color=color)
        # Plot an ellipse to show the Gaussian component
        angle = np.arctan(u[1] / u[0])
        angle = 180.0 * angle / np.pi  # convert to degrees
        ell = mpl.patches.Ellipse(mean, v[0], v[1], angle=180.0 + angle, color=color)
        ell.set_clip_box(splot.bbox)
        ell.set_alpha(0.5)
        splot.add_artist(ell)

# Load the tables
df = pd.read_csv(sys.argv[1], sep='\t')
df.columns = ['utg', 'total_kmers', 'total_maternal', 'maternal_in_utg', 'total_paternal', 'paternal_in_utg']
df['utg'] = df['utg'].str.removesuffix('_path')

X_train = df[['maternal_in_utg','paternal_in_utg']].div(df.total_kmers, axis=0).copy(deep=False)
X_train.fillna(0, inplace=True)

gm = GaussianMixture(n_components=3)
gmm = gm.fit(X_train)
mat_index = np.argmax(gmm.means_, axis=0)[0]
pat_index = np.argmax(gmm.means_, axis=0)[1]
hom_index = [0,1,2]
hom_index.remove(mat_index)
hom_index.remove(pat_index)
cluster_labels = gmm.predict(X_train)
cluster_probabilities = gmm.predict_proba(X_train)
df['cluster_labels'] = cluster_labels
df[['probability_maternal','probability_paternal','probability_homozygous']] = cluster_probabilities
pd.set_option('display.max_rows', None)
print(df)

plt.figure(0)
sns.displot(df, x="probability_maternal", bins=50)
plt.savefig("probability_maternal.png", dpi=300)
plt.close()

plt.figure(1)
sns.displot(df, x="probability_paternal", bins=50)
plt.savefig("probability_paternal.png", dpi=300)
plt.close()

plt.figure(2)
sns.displot(X_train, x="maternal_in_utg", bins=50)
plt.savefig("maternal_in_utg.png", dpi=300)
plt.close()

plt.figure(3)
sns.displot(X_train, x="paternal_in_utg", bins=50)
plt.savefig("paternal_in_utg.png", dpi=300)
plt.close()

plt.figure(4)
plot_results(np.array(X_train), cluster_labels, gmm.means_, gmm.covariances_, 0, "Gaussian Mixture")
plt.savefig("scatter_linear.png", dpi=300)
plt.close()

plt.figure(5)
plt.xscale('symlog', linthresh=0.00001)
plt.yscale('symlog', linthresh=0.00001)
colors = {mat_index: "#FF8888", pat_index: "#8888FF", hom_index[0]: "#88FF88"}
assigned_colors = [colors[0],colors[1],colors[2]]
cmap = mpl.colors.ListedColormap(assigned_colors)
plt.scatter(X_train['maternal_in_utg'], X_train['paternal_in_utg'], s=np.log(df['total_kmers']), c=cluster_labels, cmap=cmap)
centers = np.zeros(shape=(gmm.n_components, X_train.shape[1]))
for i in range(gmm.n_components):
    density = stats.multivariate_normal(cov=gmm.covariances_[i], mean=gmm.means_[i]).logpdf(X_train)
    centers[i, :] = X_train.loc[np.argmax(density)]

plt.scatter(centers[:, 0], centers[:, 1], s=20, color='black', marker='x')
plt.scatter(gmm.means_[:, 0], gmm.means_[:, 1], s=20, color='green', marker='x')

plt.savefig("scatter_log.png", dpi=300)

classified = df[['utg', 'maternal_in_utg', 'paternal_in_utg']]
classified['parental_color'] = [assigned_colors[i] for i in cluster_labels]
classified = classified.assign(maternal_percent=(X_train['maternal_in_utg']*100).round(3), paternal_percent=(X_train['paternal_in_utg']*100).round(3))
classified = classified[['utg', 'parental_color', 'maternal_in_utg', 'paternal_in_utg', 'maternal_percent', 'paternal_percent']]
classified.to_csv('hapmers.csv', index=False)
