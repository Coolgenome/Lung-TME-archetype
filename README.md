
# Spatial archetyping maps immune ecosystem transitions from lung precursor lesions to adenocarcinoma

<p align="center">
   <img width="360" height="660" alt="Cover" src="https://github.com/user-attachments/assets/3954db16-b1c1-4691-be33-85ff16cd10ab" />
</p>

## Spatial analysis in the study

1. Spatial TME cellular landscape across different disease stages during early LUAD development.
2. Spatial immune archetypes (IMPACT) presenting in LUAD and its precursor lesions.
3. Spatial architectures and validations of distinct IMPACT archetypes.
4. Immune cell functional states in different IMPACT archetypes.
5. Molecular profiles and immune-related cell-cell communications across different lesion stages and IMPACT archetypes.
6. Tumor genomic alterations and their associations with spatial TME architectures.
7. AI-based spatial decoding of lung TME architectures in histology images.
8. Analysis of TCGA cohort by spatial decoding of histology images, survival analysis.

<p align="center">
   <img width="600" height="460" alt="Study design" src="https://github.com/user-attachments/assets/83a07397-98b3-46af-b3f7-74c1a27a201d" />
</p>

## histo-IMPACT
We also provide the histo-IMPACT framework, which enables high-resolution prediction of lung lesion spatial architectures using H&E images.

<p align="center">
   <img width="500" height="450" alt="Artboard 1" src="https://github.com/user-attachments/assets/240bc805-8e56-45d6-b347-dd22d1273f37" />
</p>

This pipeline allows for identification of the tumor region, as well as detailed cell lineages (including tumor, non-tumor epithelium, and TME components). 
We provided high-resolution H&E image for demo, which can be downloaded here: https://doi.org/10.5281/zenodo.23483997
The results from demo run are included in ```bash results_histoimpact/ ```.

To run the pipeline:
```bash
cd ./histoIMPACT
bash setup_env.sh # this creates a conda environment named histoimpact
conda activate histoimpact

cd ./code
python run_histoimpact.py --image /path/to/image.tif --mpp 0.25 # set the image path and mpp here

```

## Notes
For any questions, please leave your comment in GitHub or contact Yibo Dai (ydai6@mdanderson.org). We will help address the issues as soon as possible.
