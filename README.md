# AnVIL_Galaxy to Run Experiments with GTEx Dataset and Others

This README provides a complete, detailed step-by-step guide to using AnVIL (via Terra.bio) and Galaxy to process the GTEx public dataset for gene expression analysis. We'll import the data, preprocess it in Jupyter Notebook (loading expression matrices and annotations, subsetting to tissues for classification), generate grayscale images from the expression data using a custom Python script, prepare input for a Ludwig CNN model (creating a CSV with image paths and labels, and zipping images), and finally run the Ludwig experiment in Galaxy for tissue classification (e.g., brain vs. lung vs. heart). This is based on public GTEx V8 data, which is open-access and doesn't require dbGaP approval.
The pipeline assumes you have a Terra account (sign up at https://app.terra.bio/ with a Google account) and have set up billing (required for compute, though data access is free; costs ~$0.05–$0.50/hour per VM). All steps are as of August 15, 2025—Terra/AnVIL interfaces may evolve, so check docs if issues arise.
Prerequisites

A Google account for Terra login.
Basic familiarity with Python (for Jupyter) and command-line (for gsutil).
Install nothing externally—use Terra's built-in environments.
If running locally for testing, adapt paths, but this guide focuses on Terra.

## Step 1: Access the Dataset via AnVIL Website

Open your browser and go to the AnVIL homepage: https://anvilproject.org/.
Scroll down the main page until you see the "Datasets" section (it's typically below "About" or "Consortia" sections—use Ctrl+F to search for "Datasets" if needed).
Click on the "Datasets" link or button. This redirects you to the AnVIL Data Explorer at https://explore.anvilproject.org/datasets, where you can browse available datasets.

## Step 2: Filter and Select the GTEx Public Dataset

In the AnVIL Data Explorer, look at the left sidebar (toolbar) for filters.
Expand the "Access" filter category by clicking the arrow or plus icon next to it.
Check the box for "Granted" (this shows datasets you have immediate access to, including open/public ones like GTEx; it filters out controlled-access data requiring approvals).
In the main search results area, scroll or use the search bar at the top to find "ANVIL_GTEx_Public_Data" (it appears as an archive containing GTEx V8/V9/V10 public files, such as RNA-seq expression matrices in .gct.gz format, eQTLs, histology images, and annotations like sample attributes for tissue labels).
Click on the dataset name to open its details page (this shows metadata like participant count, data types, and export options).

## Step 3: Export the Dataset to Terra

On the dataset details page, locate and click the "Export" button in the upper-right corner (it might be a dropdown or icon—hover for tooltips).
In the export dialog that appears, select the option "Analyze in Terra" (this prepares the data for import into Terra workspaces).
If prompted for cohort or filter options, mark the checkbox for "Unspecified" (this includes all public data without restrictions; if other filters appear, leave them default unless you want to subset).
Click "Request Link" (this generates a secure, temporary link for importing the data into Terra—copy it if not automatically redirected).

## Step 4: Import into Terra Workspace

Open a new browser tab and go to Terra: https://app.terra.bio/. Log in with your Google account (sign up if you haven't—it's free).
In the main Terra dashboard, click the "Workspaces" tab in the left sidebar (it lists your existing workspaces).
If you don't have a workspace, click "Create Workspace" (top-right button), give it a name like "My_GTEx_Pipeline", add a description (optional), and select a billing project (create one if needed under Profile > Billing).
Select your existing or new workspace by clicking its name.
In the workspace dashboard, look for an "Import" or "Add Data" field/button (often in the Data tab or dashboard—paste the export link from Step 3 here if prompted).
Follow the import flow: Paste the link, confirm settings, and click "Import".
Wait for the import to finish (monitor the progress bar or indicator in the upper-right corner; it may take a few minutes). Once complete, the workspace will have populated tables like:

anvil_dataset (1 row)
anvil_file (25,792 rows for public files)
anvil_project (1)
duos_dataset_registration (1)
file_inventory (25,792)
workspace_attributes (37)
These tables list files with DRS URIs (e.g., drs://...) for access.



## Step 5: Find Your Google Bucket Path
The bucket path (e.g., gs://fc-secure-your-bucket-name/) is needed for uploading images/CSV in code (replace "your-bucket-name" in scripts).

In your workspace dashboard, go to the "Data" tab.
Open a table like file_inventory or anvil_file (click its name).
Click on any file row or link to open "File Details"—this shows the full path starting with "gs://" (the bucket is the part after gs:// up to the first /, e.g., gs://fc-secure-abc123/).
Alternatively, on the dashboard, scroll to the "Google Bucket" section (bottom-left) and click "Open in browser"—this opens Google Cloud Console; the URL bar shows the gs:// path.
Copy the full bucket prefix (e.g., gs://fc-secure-abc123/)—use it in code for gsutil commands or paths.

## Step 6: Launch Jupyter Notebook and Preprocess Data

In the workspace dashboard, click the cloud icon (top-right, looks like a cloud with gear) > Select "Jupyter" environment > Choose a standard VM (4 CPUs, 15 GB RAM) > Click "Create" or "Start" (wait ~5–10 minutes for setup).
Once running, click "Open" to enter Jupyter > Create a new notebook (File > New > Notebook > Python 3 kernel).
Copy-paste and run the following preprocessing script in cells (or save as preprocessing_gtex.py and run %run preprocessing_gtex.py)

Run the script cell-by-cell or as a whole (Shift+Enter). It resolves DRS, loads data, subsets (~1,000 samples), and saves expression_matrix.csv (expression values) and metadata_base.csv (sample IDs and tissue labels). Verify prints: Matrix shape ~ (1000, 56200), columns include 'SAMPID' and 'SMTSD'.

## Step 7: Generate Images and Prepare Ludwig Input in Jupyter

In the same notebook, copy-paste and run the image generation and Ludwig prep script (or save as generate_images_and_ludwig.py and run %run generate_images_and_ludwig.py

Run the script. It generates JPG images in output_images/ (grayscale, log-normalized, padded to square based on ~56,200 genes), uploads to bucket, creates ludwig_input.csv (with gs:// image paths and tissue labels), zips images into output_images.zip, and uploads everything. Verify prints for success.

## Step 8: Run Ludwig Experiment in Galaxy

Pause Jupyter (cloud icon > Pause to save costs).
Switch to Galaxy: Cloud icon > Select "Galaxy" > Create/Start (~8–10 minutes) > Click "Open Galaxy" (launches the Galaxy interface).
Import data into Galaxy:

Click "Get Data" in the left toolbar > "Upload File" or "Paste/Fetch data".
Upload ludwig_input.csv from your local (or fetch from bucket URL: gs://your-bucket-name/ludwig_input.csv).
For images: Click Tools > Search for "Create Collection" > Add dataset collection from bucket paths (gs://your-bucket-name/output_images/*.jpg) or upload the ZIP and unzip in Galaxy (Tools > "Unzip").


Run the Ludwig experiment:

In the Tools sidebar, search for "Ludwig" (it should be available in Terra Galaxy for ML workflows; if not, check Galaxy tools or install via admin if possible).
Select "Ludwig Experiment" or "Train Ludwig Model".
Inputs:

CSV metadata file: Select ludwig_input.csv from your history.
Image collection: Select the uploaded images collection.


Config: Paste this YAML in the config input field (for CNN on images, classifying tissues)

Optional: Set train/test split (80/20), enable GPU (if VM upgraded), or add hyperparameters like dropout:0.2 in trainer.

Click "Execute" or "Submit" > Monitor the job in the History panel (right sidebar; green when done).
View results: In History, click the eye icon on outputs—check metrics like accuracy (>85% expected for tissues), confusion matrix, and predictions. Download model/predictions if needed.

## Step 9: Cleanup and Next Steps

Pause Galaxy (cloud icon > Pause) to stop billing.
Delete temporary files in bucket (via gsutil rm in terminal) if not needed.
Iterate: Rerun with different tissues (edit SELECTED_TISSUES) or normalization (e.g., 'zscore'). For other datasets (e.g., TCGA), repeat Steps 1–4 with their names in Data Explorer.
Troubleshooting: If tools missing in Galaxy, search Galaxy docs or email help@lists.anvilproject.org. For code errors, check prints and reduce samples (e.g., SAMPLES_PER_TISSUE=50 for testing).

This pipeline runs end-to-end in Terra, leveraging public GTEx for biomarker discovery or tissue classification. For advanced use, explore Galaxy workflows to automate.
