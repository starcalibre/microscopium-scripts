
# coding: utf-8

# In[50]:

# dependencies
import os
import pandas as pd
import sys
from pymongo import MongoClient
from bson.binary import Binary
from skimage import io
from skimage.transform import resize
from microscopium.screens import cellomics
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler
import numpy as np


# In[51]:

screen_name = 'LNC-SP-10A'

# load the feature data frame
def string2tuple(string_tuple):
    string_values = string_tuple.split(', ')
    coords = (int(string_values[0][1:]), string_values[1][1:-2])
    return coords

def new_tuple(coord):
    return '-'.join((screen_name, str(coord[0]), str(coord[1])))

feature_data = pd.read_csv('10A-LNC-SP-stitch2.csv', index_col=0, converters={0: string2tuple})
feature_data.head()

# update the keys to use our new format
old_index = feature_data.index
new_index = map(new_tuple, feature_data.index)
feature_data.index = new_index
feature_data.head()


# In[52]:

# normalise data frame
feature_data_std = StandardScaler().fit_transform(feature_data.values)
feature_data_std = pd.DataFrame(feature_data_std, index=feature_data.index, columns=feature_data.columns)
feature_data_std.head()


# In[53]:

# need PCA co-ordinates now
pca = PCA(n_components=2)
pca_data = pca.fit_transform(feature_data_std)
pca_data = pd.DataFrame(pca_data, index=new_index)
pca_data.head()


# In[57]:

# now a nearest neighbours dictionary
neighbours = NearestNeighbors(n_neighbors=25).fit(feature_data_std)
neighbours_dict = {}
for sample in feature_data.index:
    sample_point = feature_data_std.loc[sample].values
    distances, indices = neighbours.kneighbors(sample_point)
    neighbours_dict[sample] = list(feature_data.index[indices[0]])
    
print neighbours_dict['LNC-SP-10A-140206180002-A01']


# In[32]:

# load the images
def get_by_ext(path, extension, full=True):
    fns = []
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in [f for f in filenames if f.endswith("." + extension)]:
            fns.append(os.path.join(dirpath, filename))
    return fns

image_dir = '/media/New Volume/231-LNC-SP-stitch2-TIF'
image_fns = get_by_ext(image_dir, 'TIF')
print len(image_fns)

# convert to JPG functions
def tif2jpg(fn, outdir, new_size=None, thumb=False):
    image = io.imread(fn)
    image_base, image_ext = os.path.basename(fn).split('.')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    if thumb is True:
        new_fn = image_base + '-thumb.jpg'
    else:
        new_fn = image_base + '.jpg'
        
    if new_size is None:
        io.imsave(os.path.join(outdir, new_fn), image)
    else:
        image = resize(image, new_size)
        io.imsave(os.path.join(outdir, new_fn), image)
    return os.path.join(outdir, new_fn)

def image2key(fn):
    fn = os.path.basename(fn)
    split = fn.split('-')
    plate = str(split[1])
    well = split[2][:3]
    return '-'.join([screen_name, plate, well])

key2file = dict()
for image_fn in image_fns:
    key = image2key(image_fn)
    key2file[key] = image_fn


# In[58]:

# connect to mongo, make two connections each db
client = MongoClient('localhost', 27017)
db_myo = client['myofusion']
db_micro = client['microscopium']
myo_wells = db_myo['wells']
micro_screens = db_micro['screens']
micro_samples = db_micro['samples']
micro_images = db_micro['images']


# In[34]:

for i in range(0, feature_data.shape[0]):
    new_index  = feature_data.index[i]
    feature_vector = list(feature_data.loc[new_index].values)
    feature_vector_std = list(feature_data_std.loc[new_index].values)
    pca_vector = list(pca_data.loc[new_index].values)
    neighbours = neighbours_dict[new_index]
    image_full_fn = key2file[new_index]
    split = os.path.basename(image_full_fn).split('-')
    plate = int(split[1])
    well = split[2]
    row = well[0]
    column = well[1:3]
    
    image_full = tif2jpg(image_full_fn, '/media/New Volume/web_images/', new_size=(512, 512))
    image_thumb = tif2jpg(image_full_fn, '/media/New Volume/web_images/', new_size=(150, 150), thumb=True)
    image_full_binary = Binary(open(image_full).read())
    image_thumb_binary = Binary(open(image_thumb).read())
    
    try:
        gene_name = doc['gene_name']
    except:
        gene_name = ""
    
    new_image = {"sample_id": new_index,
                 "image_full": image_full_binary,
                 "image_thumb": image_thumb_binary}
    
    image_id = micro_images.insert(new_image)

    new_sample = {"_id": new_index,
               "gene_name": gene_name,
               "screen": screen_name,
               "control_pos": False,
               "control_neg": False,
               "feature_vector": feature_vector,
               "feature_vector_std": feature_vector_std,
               "neighbours": neighbours,
               "column": column,
               "row": row,
               "pca": pca_vector,
               "image": image_id,
               "plate": plate}
    
    micro_samples.insert(new_sample)


# In[35]:

screen = {"_id": screen_name,
          "screen_desc": "Genome-wide protein coding siRNA SMARTpool screen using colon carcinoma cells",
          "number_samples": feature_data.shape[0],
          "screen_features": list(feature_data.columns)}


# In[36]:

micro_screens.insert(screen)


# In[59]:

for i in range(0, feature_data.shape[0]):
    new_index  = feature_data.index[i]
    pca_vector = list(pca_data.loc[new_index].values)
    micro_images.update({"_id":new_index},{"$set":{"pca":pca_vector}})


# In[ ]:



