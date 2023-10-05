# Model Deployment App
The PCR-ribotype of uploaded Clostridioides Difficile genomes is predicted in real-time through this applcation.
https://github.com/dibsondivya/amr-prediction/assets/56643379/697d03d2-6f2d-4fbf-9a7d-d639ef61ca73

## Key Features
* Upload of only one of the paired end genomes is required
* Selection of k-mer size (17 or 21)
* L1 Logistic Regression used for prediction performed with weighted average F1-score of 93% for k-mer size 21 and UPLOAD% for k-mer size 17
* Real-time results for rapid analysis

## Installation
Pyhton 3.7's Flask for server side is required for the local version of this application.

## Contents
```bash
├── kmc
├── files # for temporary files  on upload, processing etc
│   ├── raw_upload
│   ├── processed_upload
│   └── dump
├── library # outputs generated from generate_columns.py
│   ├── k17_col_list.txt 
│   └── k21_col_list.txt 
├── models
│   ├── k17_all_binary_l1.joblib
│   └── k21_all_binary_l1.joblib
├── templates
│   └── index.html
├── app.py # Flask app
├── generate_columns.py  # generates names of library k-mers from uploaded .csv (where k-mers are columns)
├── requirements.txt
└── README.md
```

## Project Status
Project has yet to be deployed.

## Future Changes
* Allow upload of file of genomes and output of text file of predicted ribotypes

## Contact
If iin need of assistance, please do reach out at dds122@ic.ac.uk.
