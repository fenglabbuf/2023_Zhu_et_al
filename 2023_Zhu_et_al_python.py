import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import rcParams
import numpy as np
from scipy.cluster import hierarchy
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

plt.rcParams.update(plt.rcParamsDefault)
rcParams['font.family'] = 'arial'
rcParams.update({'font.size': 13})

### heatmaps for Amnp samples ###
gaba = ['DLX1', 'DLX2', 'DLX5','DLX6','LHX6','GBX1', 'GBX2', 'OTX2',
        'LAMP5','VIP','GAD1','GAD2','SLC6A1','SLC17A8','SLC32A1']
fib = ['DCN','FN1','FBLN1','FBLN2','DMKN','HSPB1','PTX3','TGFBI','COL1A2','COL13A1',
       'S100A6','S100A4','LGALS1','LOXL1','KRT5']

df_amnp = pd.read_csv('./../results/Amnp_day3_day14.csv').set_index('Row.names').filter(regex='Amnp')
gabam = np.log2(df_amnp.loc[gaba].to_numpy().reshape((len(gaba)*3,3)).mean(axis=1).reshape((len(gaba),3))+1)
fibm = np.log2(df_amnp.loc[fib].to_numpy().reshape((len(fib)*3,3)).mean(axis=1).reshape((len(fib),3))+1)

# GABA markers
fig = sns.clustermap(gabam,yticklabels=df_amnp.loc[gaba].index, col_cluster=False, row_cluster=False,
                    xticklabels=['Day3', 'Day7', 'Day14'],z_score=0,vmax=1.2, vmin=-1.2,
                     cmap='magma', figsize=(2,5.1), method='complete',linewidth=0.2)
fig.ax_heatmap.tick_params(bottom=False, right=False)
fig.ax_row_dendrogram.set_visible(False)

# fibroblast markers
fig = sns.clustermap(fibm,yticklabels=df_amnp.loc[fib].index, col_cluster=False, row_cluster=False,
                    xticklabels=['Day3', 'Day7', 'Day14'],z_score=0,vmax=1.2, vmin=-1.2,
                     cmap='magma', figsize=(2,5.1), method='complete',linewidth=0.2)
fig.ax_heatmap.tick_params(bottom=False, right=False)
fig.ax_row_dendrogram.set_visible(False)


### line plots for longitudinal gene expression ###
# read RNA binding protein (RBP) list
rbps = pd.read_csv('./../data/RBPgenes.csv')
# read RBP gene expression in samples
df_mtx = pd.read_csv('./../results/Amnp_day3_day14.csv').set_index('Row.names')
df_mtx_rbp = df_mtx[df_mtx.index.isin(rbps['Approved symbol'])]

anp = df_mtx_rbp.filter(regex='Anp')
amp = df_mtx_rbp.filter(regex='Amp')
amnp = df_mtx_rbp.filter(regex='Amnp')

# highlighted RBP list
lstsp = ['RBFOX1', 'RBFOX3','KHDRBS2','LIN28B','ELAVL2', 
        'ELAVL4','CELF4','CELF5','CPEB3','NOVA2']

# background RBP expression
anp_g = anp.drop(lstsp)
amnp_g = amnp.drop(lstsp)
amp_g = amp.drop(lstsp)

mtxanp = np.log2(anp_g.to_numpy().reshape((len(anp_g)*3,3)).mean(axis=1).reshape((len(anp_g),3))+1)
mtxamnp = np.log2(amnp_g.to_numpy().reshape((len(amnp_g)*3,3)).mean(axis=1).reshape((len(amnp_g),3))+1)
mtxamp = np.log2(amp_g.to_numpy().reshape((len(amp_g)*3,3)).mean(axis=1).reshape((len(amp_g),3))+1)

# RBP expression in Anp
fig = plt.figure(figsize=(3,6))
ax = fig.add_subplot(1, 1, 1)
for n in range(len(anp_g)):
    plt.plot(np.arange(3), mtxanp[n,:], color='grey', alpha=0.3)
a = 0.9
plt.plot(np.arange(3), np.log2(anp.loc['CELF4'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:purple',label='CELF4', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['CELF5'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:olive',label='CELF5', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['CPEB3'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='magenta',label='CPEB3', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['ELAVL2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:green',label='ELAVL2', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['ELAVL4'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:red',label='ELAVL4', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['KHDRBS2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:orange',label='KHDRBS2', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['LIN28B'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='darkred',label='LIN28B', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['NOVA2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='darkblue',label='NOVA2', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['RBFOX1'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:blue', label='RBFOX1', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(anp.loc['RBFOX3'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:cyan',label='RBFOX3', linewidth=3.7, linestyle='--', alpha=a)
plt.ylim(bottom=0, top=14)
plt.xticks([0,1,2], ['3', '7', '14'], fontsize=16)
plt.ylabel('Log2 relative expression', fontsize=16)
plt.xlabel('Day', fontsize=16)
plt.title('Anp', fontsize=16)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

# RBP expression in Amp
fig = plt.figure(figsize=(3,6))
ax = fig.add_subplot(1, 1, 1)
for n in range(len(amp_g)):
    plt.plot(np.arange(3), mtxamp[n,:], color='grey', alpha=0.3)
a = 0.9
plt.plot(np.arange(3), np.log2(amp.loc['CELF4'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:purple',label='CELF4', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['CELF5'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:olive',label='CELF5', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['CPEB3'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='magenta',label='CPEB3', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['ELAVL2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:green',label='ELAVL2', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['ELAVL4'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:red',label='ELAVL4', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['KHDRBS2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:orange',label='KHDRBS2', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['LIN28B'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='darkred',label='LIN28B', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['NOVA2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='darkblue',label='NOVA2', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['RBFOX1'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:blue', label='RBFOX1', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amp.loc['RBFOX3'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:cyan',label='RBFOX3', linewidth=3.7, alpha=a)
plt.ylim(bottom=0, top=14)
plt.xticks([0,1,2], ['3', '7', '14'], fontsize=16)
plt.ylabel('Log2 relative expression', fontsize=16)
plt.xlabel('Day', fontsize=16)
plt.title('AMp', fontsize=16)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')

# RBP expression in Amnp
fig = plt.figure(figsize=(3,6))
ax = fig.add_subplot(1, 1, 1)
for n in range(len(amnp_g)):
    plt.plot(np.arange(3), mtxamnp[n,:], color='grey', alpha=0.3)
a=0.9
plt.plot(np.arange(3), np.log2(amnp.loc['CELF4'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:purple',label='CELF4', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['CELF5'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:olive',label='CELF5', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['CPEB3'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='magenta',label='CPEB3', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['ELAVL2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:green',label='ELAVL2', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['ELAVL4'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:red',label='ELAVL4', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['KHDRBS2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:orange',label='KHDRBS2', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['LIN28B'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='darkred',label='LIN28B', linewidth=3.7, alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['NOVA2'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='darkblue',label='NOVA2', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['RBFOX1'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:blue', label='RBFOX1', linewidth=3.7, linestyle='--', alpha=a)
plt.plot(np.arange(3), np.log2(amnp.loc['RBFOX3'].to_numpy().reshape((3,3)).mean(axis=1).reshape(3)+1), 
         color='tab:cyan',label='RBFOX3', linewidth=3.7, alpha=a)
plt.ylim(bottom=0, top=14)
plt.legend(loc='upper right', bbox_to_anchor=(1.7,1.0), frameon=False)
plt.xticks([0,1,2], ['3', '7', '14'], fontsize=16)
plt.ylabel('Log2 relative expression', fontsize=16)
plt.xlabel('Day', fontsize=16)
plt.title('AMnp', fontsize=16)
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')


### heatmaps for RBP in GTEx ###
gtex = pd.read_csv('./../data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct', sep='\t', skiprows=2)
rbps = pd.read_csv('./../data/RBPgenes.csv')

map = rbps.merge(gtex, left_on='Approved symbol', 
                 right_on='Description').drop(columns=['Name', 'Approved symbol']).set_index('Description')

dcolor = {'Adipose':'#ff0000','Artery':'#ff8900','Brain':'#fcd900',
          'Cerebellum':'#be8600','Cells':'#bafc03','Cervix':'#03fc1c',
          'Colon':'#16f09c','Esophagus':'#00c6ff','Heart':'#003eff',
          'Kidney':'#4b00ff','Skin':'#d400ff'}

xlabel = []
xcolor = []
for n in map.columns:
    if ' - ' in n:
        if n.split(' - ')[0] in ['Breast','Muscle', 'Nerve', 'Small Intestine']:
            xlabel.append(n.split(' - ')[0])
            xcolor.append('#ffffff')
        elif ' (' in n:
            xlabel.append(n.split(' - ')[1].split(' (')[0])
            xcolor.append(dcolor[n.split(' - ')[0]])
        else:
            xlabel.append(n.split(' - ')[1])
            if 'Cerebell' in n:
                xcolor.append(dcolor['Cerebellum'])
            else:
                xcolor.append(dcolor[n.split(' - ')[0]])
    else:
        xlabel.append(n)
        xcolor.append('#ffffff')

fig = sns.clustermap(np.log2(map+0.5).to_numpy(), yticklabels=[], 
               xticklabels=xlabel, cmap='magma', figsize=(10,20), col_colors=xcolor)


### dendrogram of samples based on differentially spliced junctions (DSJs) of RNA ###
def readcluster(prefix, samples=['AMp', 'AMnp']):
    sg = pd.read_csv(prefix+'_cluster_significance.txt', sep='\t')
    ef = pd.read_csv(prefix+'_effect_sizes.txt', sep='\t')
    sg['cluster_n'] = [n.split('_')[1] for n in list(sg['cluster'])]
    ef['cluster_n'] = [n.split('_')[1] for n in list(ef['intron'])]
    efsg = ef.merge(sg, left_on='cluster_n', right_on='cluster_n', how='left')
    efsg['loc'] = [n.split('_')[0][:-4] for n in list(efsg['intron'])]
    return efsg[['loc', 'deltapsi']+samples+['p.adjust', 'cluster']]

cd3 = readcluster('./../results/STAR/leafcutter/SJ_Amp_Amnp_3')
cd7 = readcluster('./../results/STAR/leafcutter/SJ_Amp_Amnp_7')
cd14 = readcluster('./../results/STAR/leafcutter/SJ_Amp_Amnp_14')
cdall = cd3.merge(cd7, left_on='loc', right_on='loc').merge(cd14, left_on='loc', right_on='loc')

gtesb = readcluster('./../results/STAR/leafcutter/GTEx/fibroblast_brain', samples=['Skin', 'Brain'])
gtesb = gtesb[(gtesb['p.adjust']<0.05)&(gtesb['deltapsi'].abs()>0.1)].sort_values(by='p.adjust')\
        .sort_values(by='deltapsi', key=abs, ascending=False).rename(columns={'deltapsi':'gtex_deltapsi'})

names = pd.read_csv('./../data/gene_names.csv')

SJmap = cdall.merge(gtesb, left_on='loc', right_on='loc')\
        .sort_values(by='gtex_deltapsi', ascending=False)\
        .merge(names, left_on='clusterID', right_on='clusterID', how='left').set_index('loc')\
        [['Skin','AMp_x', 'AMnp_x', 'AMp_y', 'AMnp_y', 'AMp', 'AMnp', 'Brain']]

fig = sns.clustermap(SJmap.to_numpy(),yticklabels=[],
                    xticklabels=['Fibroblast', 'Day3 AMp', 'Day3 AMnp', 'Day7 AMp', 'Day7 AMnp', 'Day14 AMp', 'Day14 AMnp', 'Brain'], 
                     cmap='magma', figsize=(4,10),vmin=0,vmax=1, method='complete',
                     tree_kws=dict(linewidth=1.5))

rcParams.update({'font.size': 16})
rcParams['lines.linewidth'] = 2.5

fig1 = plt.figure(figsize=(4, 5))
ax = fig1.add_subplot(1, 1, 1)
hierarchy.set_link_color_palette(['r', 'b'])
dn = hierarchy.dendrogram(fig.dendrogram_col.linkage, orientation='left', above_threshold_color='black',
                         labels = ['Fibroblast', 'Day3 AMp', 'Day3 AMnp', 'Day7 AMp', 
                                   'Day7 AMnp', 'Day14 AMp', 'Day14 AMnp', 'Brain'])
plt.xlabel('Height')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['left'].set_color('none')
plt.yticks(fontsize=16);


### PCA plot of samples based on DSJs ###
pc_scale = StandardScaler().fit_transform(SJmap.to_numpy().T)
pca = PCA(n_components=2)
pca.fit(pc_scale)
pcav = pca.transform(pc_scale)
pc_r = pca.explained_variance_ratio_*100

fig = plt.figure(figsize=(5, 5))
ax = fig.add_subplot(1, 1, 1)

plt.scatter(pcav[0,0], pcav[0,1], s=90, alpha=0.6, marker='o', color='r', label='Fibroblast')
plt.scatter(pcav[-1,0], pcav[-1,1], s=90, alpha=0.6, marker='o', color='b', label='Brain')
plt.scatter(pcav[1,0], pcav[1,1], s=100, alpha=0.7, marker='^', color='#ff2e7e', label='Day3 AMp')
plt.scatter(pcav[2,0], pcav[2,1], s=100, alpha=0.7, marker='P', color='#ff2e7e', label='Day3 AMnp')
plt.scatter(pcav[3,0], pcav[3,1], s=100, alpha=0.7, marker='^', color='#a21ddb', label='Day7 AMp')
plt.scatter(pcav[4,0], pcav[4,1], s=100, alpha=0.7, marker='P', color='#a21ddb', label='Day7 AMnp')
plt.scatter(pcav[5,0], pcav[5,1], s=100, alpha=0.7, marker='^', color='#4f34eb', label='Day14 AMp')
plt.scatter(pcav[6,0], pcav[6,1], s=100, alpha=0.7, marker='P', color='#4f34eb', label='Day14 AMnp')
plt.xlabel('PC1 (' + str(pc_r[0])[:4] + '%)', fontsize=16)
plt.ylabel('PC2 (' + str(pc_r[1])[:4] + '%)', fontsize=16)
plt.grid(linewidth=2, alpha=0.2, color='grey', linestyle=':')
plt.legend(bbox_to_anchor=(1.01,1), prop={'size':14})