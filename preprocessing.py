# full_gtex_pipeline.py
# This script resolves DRS URIs for GTEx V8 public data files, loads the gene expression matrix and annotations,
# subsets to selected tissues, generates images from the matrix, creates ludwig_input.csv with image_path and label,
# zips the JPG images, and uploads files to the bucket. It includes progress prints for feedback.

import pandas as pd
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from terra_notebook_utils import drs  # Ensure terra-notebook-utils is installed: pip install --upgrade terra-notebook-utils
import zipfile
import shutil
import sys  # For flushing prints
import time  # For timing prints

# Define DRS URIs from the file_inventory table
TPM_DRS = 'drs://drs.anv0:v2_13903d0b-9dc6-3e09-a6f9-0275f7e5d78f'  # GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
READS_DRS = 'drs://drs.anv0:v2_9694ab46-2c6f-3a0c-ab00-732699214c03'  # GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz (optional)
ANNOT_DRS = 'drs://drs.anv0:v2_cc9c26d1-2ba8-3a0a-85ac-d06c72a6d77a'  # GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

# Selected tissues for subsetting (adjust as needed)
SELECTED_TISSUES = ['Brain - Cortex', 'Heart - Left Ventricle', 'Liver', 'Lung', 'Muscle - Skeletal']
SAMPLES_PER_TISSUE = 200  # Number of samples per tissue

# Image parameters
NORMALIZATION_METHOD = 'log'
JPG_DIR = 'output_images'

# Your bucket path (replace with your actual bucket, e.g., 'gs://fc-secure-abc123/')
BUCKET = 'gs://your-bucket-name/'

def resolve_drs(drs_uri):
    """Resolve DRS URI to a signed URL."""
    start_time = time.time()
    print(f"Resolving DRS URI: {drs_uri}...")
    url = drs.access(drs_uri)
    elapsed = time.time() - start_time
    print(f"Resolved DRS URI in {elapsed:.2f} seconds.")
    return url

def load_gct(signed_url):
    """Load .gct.gz file into DataFrame and transpose."""
    start_time = time.time()
    print("Loading GCT file...")
    df = pd.read_csv(signed_url, sep='\t', skiprows=2, index_col=0, compression='gzip')
    df = df.iloc[:, 1:].transpose()  # Samples as rows, genes as columns
    elapsed = time.time() - start_time
    print(f"GCT file loaded and transposed in {elapsed:.2f} seconds.")
    return df

def load_annotations(signed_url):
    """Load annotations TXT file."""
    start_time = time.time()
    print("Loading annotations file...")
    df = pd.read_csv(signed_url, sep='\t')
    elapsed = time.time() - start_time
    print(f"Annotations file loaded in {elapsed:.2f} seconds.")
    return df

def preprocess_data():
    start_time = time.time()
    print("Starting preprocessing (Step 1/4)...")
    # Resolve URLs
    tpm_signed_url = resolve_drs(TPM_DRS)
    # reads_signed_url = resolve_drs(READS_DRS)  # Uncomment if using reads instead
    annot_signed_url = resolve_drs(ANNOT_DRS)

    # Load data (use reads_signed_url if preferring counts)
    df_tpm = load_gct(tpm_signed_url)
    print(f"Loaded TPM matrix with shape: {df_tpm.shape}")

    annot = load_annotations(annot_signed_url)
    print(f"Loaded annotations with columns: {annot.columns.tolist()}")

    # Subset to selected tissues
    subset_start = time.time()
    print(f"Subsetting to {len(SELECTED_TISSUES)} tissues with {SAMPLES_PER_TISSUE} samples each...")
    meta = annot[annot['SMTSD'].isin(SELECTED_TISSUES)].groupby('SMTSD').sample(SAMPLES_PER_TISSUE, random_state=42)
    # Filter meta to only include samples present in matrix index to avoid KeyError
    meta = meta[meta['SAMPID'].isin(df_tpm.index)]
    print(f"Filtered to {len(meta)} common samples after index check.")
    df_filtered = df_tpm.loc[meta['SAMPID']]
    subset_elapsed = time.time() - subset_start
    print(f"Subset complete in {subset_elapsed:.2f} seconds.")

    # Save outputs
    save_start = time.time()
    df_filtered.to_csv('expression_matrix.csv', index=True)
    meta[['SAMPID', 'SMTSD']].to_csv('metadata_base.csv', index=False, header=['sample_id', 'label'])
    save_elapsed = time.time() - save_start
    print(f"Saved expression_matrix.csv and metadata_base.csv in {save_elapsed:.2f} seconds.")
    
    total_elapsed = time.time() - start_time
    print(f"Preprocessing complete (Step 1/4 done) in {total_elapsed:.2f} seconds.")

def load_matrix(file_path, sep=","):
    return pd.read_csv(file_path, sep=sep, index_col=0)

def normalize_row(row, method="log"):
    if method == "log":
        return np.log1p(row)
    elif method == "minmax":
        return (row - np.min(row)) / (np.max(row) - np.min(row) + 1e-8)
    elif method == "zscore":
        return (row - np.mean(row)) / (np.std(row) + 1e-8)
    else:
        return row

def get_image_size(num_genes):
    return math.ceil(math.sqrt(num_genes))

def row_to_image(row, image_size, method="log"):
    norm_row = normalize_row(row, method)
    padded = np.zeros(image_size * image_size, dtype=np.float32)
    padded[:len(norm_row)] = norm_row
    return padded.reshape((image_size, image_size))

def convert_matrix_to_images(df, normalization="log"):
    num_samples, num_genes = df.shape
    image_size = get_image_size(num_genes)
    image_tensor = np.zeros((num_samples, image_size, image_size, 1), dtype=np.float32)
    for i, row in enumerate(df.values):
        image_tensor[i, ..., 0] = row_to_image(row, image_size, normalization)
        if (i + 1) % 100 == 0:  # Progress feedback every 100 samples
            print(f"Converted {i + 1}/{num_samples} rows to images ({((i + 1) / num_samples * 100):.2f}% done)...", flush=True)
    print("All rows converted to images.")
    return image_tensor, image_size

def save_images_as_jpg(tensor, image_dir, sample_names=None):
    os.makedirs(image_dir, exist_ok=True)
    num_images = tensor.shape[0]
    for i, img in enumerate(tensor):
        img_2d = img[:, :, 0]
        sample_id = sample_names[i] if sample_names is not None else f"sample_{i}"
        output_path = os.path.join(image_dir, f"{sample_id}.jpg")
        plt.imsave(output_path, img_2d, cmap="gray", format="jpg")
        if (i + 1) % 100 == 0:  # Progress feedback every 100 images
            print(f"Saved {i + 1}/{num_images} images ({((i + 1) / num_images * 100):.2f}% done)...", flush=True)
    print(f"Saved {num_images} JPEG images to {image_dir}")

def generate_images(matrix_path, jpg_dir, normalization):
    print("Starting image generation (Step 2/4)...")
    df = load_matrix(matrix_path)
    tensor, image_size = convert_matrix_to_images(df, normalization=normalization)
    save_images_as_jpg(tensor, jpg_dir, sample_names=df.index)
    print("Image generation complete (Step 2/4 done).")

def create_ludwig_csv(metadata_path, jpg_dir, bucket):
    print("Creating Ludwig input CSV...")
    meta = pd.read_csv(metadata_path)
    meta['image_path'] = [row['sample_id'] + '.jpg' for _, row in meta.iterrows()]  # Only image name, no full path
    meta[['image_path', 'label']].to_csv('ludwig_input.csv', index=False)
    print("Saved ludwig_input.csv")

def zip_images(jpg_dir, zip_path):
    print("Zipping images...")
    num_files = len([f for f in os.listdir(jpg_dir) if f.endswith('.jpg')])
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        for root, _, files in os.walk(jpg_dir):
            for j, file in enumerate(files):
                zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), jpg_dir))
                if (j + 1) % 100 == 0:  # Progress feedback every 100 files
                    print(f"Zipped {j + 1}/{num_files} files ({((j + 1) / num_files * 100):.2f}% done)...", flush=True)
    print(f"Created ZIP file: {zip_path}")

def upload_to_bucket(local_path, bucket_path):
    print(f"Uploading {local_path} to bucket...")
    os.system(f"gsutil cp -r {local_path} {bucket_path}")
    print(f"Uploaded {local_path} to {bucket_path}")

def full_pipeline():
    preprocess_data()  # Step 1: Preprocess and save matrix/metadata

    generate_images('expression_matrix.csv', JPG_DIR, NORMALIZATION_METHOD)  # Step 2: Generate images

    print("Starting Ludwig preparation and zipping (Step 3/4)...")
    # Upload images to bucket
    upload_to_bucket(JPG_DIR, BUCKET)

    # Create Ludwig CSV
    create_ludwig_csv('metadata_base.csv', JPG_DIR, BUCKET)

    # Upload Ludwig CSV to bucket
    upload_to_bucket('ludwig_input.csv', BUCKET)

    # Zip images
    zip_images(JPG_DIR, 'output_images.zip')

    # Optional: Upload ZIP to bucket
    upload_to_bucket('output_images.zip', BUCKET)
    print("Ludwig preparation and zipping complete (Step 3/4 done).")

    print("Starting cleanup (Step 4/4)...")
    # Cleanup local images (optional, to save space)
    shutil.rmtree(JPG_DIR)
    print("Cleanup complete. Pipeline finished (Step 4/4 done). All files ready for Ludwig.")

if __name__ == "__main__":
    full_pipeline()
