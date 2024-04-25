#!/usr/bin/python 3

import sklearn as sk
import numpy as np
import skbio
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn import preprocessing
from sklearn.manifold import TSNE
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from skbio.stats.ordination import pcoa 
from skbio import DistanceMatrix
from sklearn.metrics import pairwise_distances

#from PHIbase mapping 
raw_data = pd.read_csv('../taxonomy_totalread.csv')

seq_info = pd.read_csv('../sequencing_details_samp_comp.csv')

# Create a column of the read counts divided by the read counts for that barcode - normalised
raw_data['read_norm_by_total'] = raw_data.read_count/raw_data.total_read

#Extracting the genus information
# creating a table of just genus & barcode
genus_data = raw_data[['genus', 'barcode']]

#because this is a binary presence/absence matrix i only want unique genera for each barcode
#subset means we are dropping rows that have duplicate values in both barcode and genus col
genus_data = genus_data.drop_duplicates(subset=['barcode', 'genus'])

#reading in MEGAN data
#meg_count_data = pd.read_csv('../MEGAN_taxonNametoCount_sum_SampComp.csv')
meg_count_data = pd.read_csv('../MEGAN_taxonNametoCount_assi_SampComp.csv')
#renaming the samples
new_columns = [col[-2:] for col in meg_count_data.columns[1:]]
# Replace the column names with the extracted numbers
meg_count_data.columns = ["barcode"] + new_columns


#Creating the matrix, for a t-SNE it's a case of presence or absence of each taxa
#genera as columns and barcodes and rows
#aggfunc - assigns a value of 1 to each cell in the binary matrix
#fill_value - any missing values in the matrix set to 0
#binary_matrix = pd.pivot_table(genus_data, index='barcode', columns='genus', aggfunc=lambda x: 1, fill_value=0)

#For MEGAN data it's slightly diff, this isn't just 0 and 1 but I'm going to try it
# Pivot the DataFrame to create a matrix
binary_matrix = pd.pivot_table(meg_count_data, index=None, columns='barcode')




# Perform T-SNE on your binary matrix---------------------------------
#n_components - data is reduced to 2 components for x and y
tsne = TSNE(n_components=2, random_state=0)
# #T-SNE co-ordinates are stored in t-sne_results
tsne_results = tsne.fit_transform(binary_matrix)

# Extract T-SNE coordinates for x and y axes
tsne_x = tsne_results[:, 0]
tsne_y = tsne_results[:, 1]

#colouring by location
#barcode 1 - 20 NHM and 25 - 44 Church farm
your_color_data = ['#008080' if i < 20 else '#FF7F50' for i in range(40)]

# Add labels to the points - just barcode numbers, to check it's being plotted correctly
#doesn't work for megan data as it's not indexed
#for i, barcode in enumerate(binary_matrix.index):
#    plt.text(tsne_x[i], tsne_y[i], barcode)


#sampler data
sampler_barcode = seq_info[['Sampler', 'barcode']]

# Create a dictionary to map sampler labels to marker symbols
marker_data = {
    'Bobcat': 'o',  # circle marker
    'Cub': '+',     # Use '+' for plus marker
    'SASS': 'x',     # Use 'x' for x marker (doesn't quite match r aesthetics)
    'Coriolis_Micro': 's', # Use 's' for square marker
    'Coriolis_Compact': '^' # Use '^' for upward triangle marker
}

# Create custom legend handles and labels
legend_handles = []
legend_labels = []

#Creating the legend to be black symbols based on my marker_data dict
for sampler in marker_data.keys():
    marker = Line2D([0], [0], marker=marker_data[sampler], color='black', markersize=10, label=sampler)
    legend_handles.append(marker)
    legend_labels.append(sampler)

# Plotting scatter t-sne
# Add marker symbols for samplers using sampler_barcode list & marker_data dictionary
for i in range(tsne_results.shape[0]):
    barcode = i + 1  # Barcode is assumed to be an index from 1 to n
    
    if barcode in sampler_barcode['barcode'].values:
        #barcode from sampler_barcode compared to tsne_results - creates a boolean array
        #boolean array is used to index the sampler_barcode dataframe
        sampler = sampler_barcode[sampler_barcode['barcode'] == barcode]['Sampler'].values[0]
        
        #Using above loop to specify the shape of each point 
        if sampler in marker_data:
            # Use modulo operator to loop through marker colors - specifed above
            color = your_color_data[(barcode - 1) % len(your_color_data)]
            #make sure marker corresponds to symbol obtained from marker data
            marker = marker_data[sampler]
            plt.scatter(tsne_x[i], tsne_y[i], marker=marker, color=color, s=75, alpha=0.75)
            # add barcode number to the point
            #plt.text(tsne_x[i], tsne_y[i], barcode)
            


#Setting aesthetics
# Create a legend with custom handles and labels
# plt.legend(handles=legend_handles, labels=legend_labels, loc='upper left', bbox_to_anchor=(1.01, 1), handlelength=0)
# plt.xlabel('T-SNE Dimension 1')
# plt.ylabel('T-SNE Dimension 2')
# plt.title('t-SNE of Genus-level taxonomic diversity (MEGAN assigned counts)')
# #plt.show()

# plt.savefig('/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Sampler comparison/images/t-SNE_Genus_MEGAN_SampComp_assi.pdf', bbox_inches='tight')


#Perform a PCoA on the binary matrix----------------------------------

# Compute dissimilarity matrix using Jaccard distance - expects a numpy array
#trying both methods
#binary_matrix = np.array(binary_matrix)
dissimilarity_matrix = pairwise_distances(binary_matrix, metric='jaccard')

# Perform PCoA on dissimilarity matrix
pcoa_results = pcoa(dissimilarity_matrix)

# Access the PCoA results
eigvals = pcoa_results.eigvals
samples = pcoa_results.samples

# Scale the samples for better visualization
scale_factor = 10
scaled_samples = samples * scale_factor

# Access the PCoA scores from the OrdinationResults object
pcoa_scores = pcoa_results.samples.values

# Plot the first two principal coordinates with colours and shapes
for i in range(pcoa_scores.shape[0]):
    barcode = i + 1  # Barcode is assumed to be an index from 1 to n
    
    if barcode in sampler_barcode['barcode'].values:
        #barcode from sampler_barcode compared to tsne_results - creates a boolean array
        #boolean array is used to index the sampler_barcode dataframe
        sampler = sampler_barcode[sampler_barcode['barcode'] == barcode]['Sampler'].values[0]
        
        #Using above loop to specify the shape of each point 
        if sampler in marker_data:
            # Use modulo operator to loop through marker colors - specifed above
            color = your_color_data[(barcode - 1) % len(your_color_data)]
            #make sure marker corresponds to symbol obtained from marker data
            marker = marker_data[sampler]
            plt.scatter(pcoa_scores[i, 0], pcoa_scores[i, 1], marker=marker, color=color, s=75, alpha=0.75)

# plt.legend(handles=legend_handles, labels=legend_labels, loc='upper left', bbox_to_anchor=(1.01, 1), handlelength=0) 
# plt.xlabel("PC1")
# plt.ylabel("PC2")
# plt.title("PCoA of Dissimilarity Matrix using Jaccard distance (points scaled x10, Genus, MEGAN, Summarised counts)")
# #plt.savefig('/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Sampler comparison/images/PCoA_MEGAN_Genus_SampComp_Jaccard_sum.pdf', bbox_inches='tight')
# plt.show()


#Bray curtis distance using a different matrix--------------
#Using my raw_data from the top of the script and keeping the norm read counts

# Group by genus and sum the read counts
#genus_df = raw_data.groupby(['genus', 'barcode']).agg({'read_norm_by_total': 'sum'})

# Pivot the grouped DataFrame to convert it into a matrix
# abundance_matrix = genus_df.pivot_table(index='genus', columns='barcode', values='read_norm_by_total')

# Fill missing values with 0
# abundance_matrix = abundance_matrix.fillna(0)

# Calculate pairwise Bray-Curtis dissimilarity
# bc_matrix = pairwise_distances(abundance_matrix, metric='braycurtis')
#same but can use MEGAN binary matrix
bc_matrix = pairwise_distances(binary_matrix, metric='braycurtis')

# Perform PCoA on dissimilarity matrix
bc_pcoa_results = pcoa(bc_matrix)

# Access the PCoA results
bc_eigvals = bc_pcoa_results.eigvals
bc_samples = bc_pcoa_results.samples

# Access the PCoA scores from the OrdinationResults object
bc_pcoa_scores = bc_pcoa_results.samples.values

# Plot the first two principal coordinates with colours and shapes
for i in range(bc_pcoa_scores.shape[0]):
    barcode = i + 1  # Barcode is assumed to be an index from 1 to n
    
    if barcode in sampler_barcode['barcode'].values:
        #barcode from sampler_barcode compared to tsne_results - creates a boolean array
        #boolean array is used to index the sampler_barcode dataframe
        sampler = sampler_barcode[sampler_barcode['barcode'] == barcode]['Sampler'].values[0]
        
        #Using above loop to specify the shape of each point 
        if sampler in marker_data:
            # Use modulo operator to loop through marker colors - specifed above
            color = your_color_data[(barcode - 1) % len(your_color_data)]
            #make sure marker corresponds to symbol obtained from marker data
            marker = marker_data[sampler]
            plt.scatter(bc_pcoa_scores[i, 0], bc_pcoa_scores[i, 1], marker=marker, color=color, s=75, alpha=0.75)

plt.legend(handles=legend_handles, labels=legend_labels, loc='upper left', bbox_to_anchor=(1.01, 1), handlelength=0) 
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCoA of Genus-level taxonomic diversity (MEGAN data assigned counts)")
plt.savefig('/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Sampler comparison/images/PCoA_MEGAN_Genus_SampComp_BC_assi.pdf', bbox_inches='tight')
# plt.show()