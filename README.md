
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

<p align="center">
   <img width="600" height="460" alt="Study design" src="https://github.com/user-attachments/assets/83a07397-98b3-46af-b3f7-74c1a27a201d" />
</p>

## LungTME-Net
We also provide the LungTME-Net framework, which enables high-resolution prediction of lung lesion spatial architectures using H&E images.

<p align="center">
   <img width="600" height="300" alt="LungTME-Net" src="https://github.com/user-attachments/assets/8b2a2970-3e47-4cb0-9ee0-655f28aa8187" />
</p>

This pipeline allows for identification of the tumor region, as well as detailed cell lineages (including tumor, non-tumor epithelium, and TME components). 
We provided three high-resolution H&E image patches for demo: 
  - ADC1 is a tumor-region-dominant patch with rich TME components, with abundant tumor cells, as well as lymphocytes, macrophages, and fibrous tissue; 
  - ADC2 is a tumor-region-dominant patch, with tumors along with vessels; 
  - ADC3 is the tumor-normal interface.
    
The three demo patches include multiple distinct tissue structures and demonstrate the performance of the model in different lung tissue contexts. The users can also run the pipeline using their own high-resolution lung H&E images. The spatial prediction results can be leveraged for various downstream tasks based on the users' needs, such as spatial pattern, niche, and archetype analyses and correlative analyses using clinical metadata.

To run the pipeline:
```bash
cd ./LungTME-Net
python3.12 -m venv .venv # create an environment for running the model, the model has been tested in python 3.12
source .venv/bin/activate # activate the environment
pip install -r requirements.txt # install essential packages

cd ./code
python3 demo_workflow.py # run the demo pipeline will give output for the three examples in the /results/ folder
```

## Notes
For any questions, please leave your comment in GitHub or contact Yibo Dai (ydai6@mdanderson.org). We will help address the issues as soon as possible.
