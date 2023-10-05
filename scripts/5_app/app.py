from flask import Flask, render_template, request, redirect, url_for, flash
import subprocess
import joblib
import sys
import collections
import csv
import sklearn
import pandas as pd

app = Flask(__name__)
app.secret_key = 'super secret key'

# load required models and library files
model17 = joblib.load('models/k17_all_binary_l1.joblib')
model21 = joblib.load('models/k21_all_binary_l1.joblib')
file_content17 = open("library/k17_col_list.txt", "r") 
data17 = file_content17.read() 
kmers17 = data17.split("\n")
file_content21 = open("library/k21_col_list.txt", "r") 
data21 = file_content21.read() 
kmers21 = data21.split("\n")


@app.route('/', methods=["GET", "POST"])
def index():
    if request.method == 'POST':
        # form import
        k = request.form['inlineRadioOptions']
        file = request.files['file']
        filename=file.filename
        file.save('file/raw_upload/'+filename)
        # selecting required model and library files
        if k=='17':
            model = model17
            kmers_list = kmers17[:-1] # to remove last element that is empty ''
            kmers_dict = dict.fromkeys(kmers_list)
        elif k=='21':
            model = model21
            kmers_list = kmers21[:-1]
            kmers_dict = dict.fromkeys(kmers_list)
        else:
            message = 'There is an error with the model chosen.'

        naming = filename.strip('.fastq.gz')
        # run kmer count
        checking = False
        count = 0
        while checking == False and count < 3:
            try:
                count += 1
                subprocess.check_output('kmc/bin/kmc -k' + k + ' file/raw_upload/' + filename + ' file/processed_upload/' + naming + ' file/temporary', shell=True)
                subprocess.check_output('kmc/bin/kmc_tools transform file/processed_upload/'+naming+ ' dump file/dump/' + naming+'.txt', shell=True)
                checking = True
            except subprocess.CalledProcessError as error:
                message = 'There is an error with the file uploaded. ' + str(error)
                checking = False
        
        if checking == True:
            # prep for model input
            read_kmer_count = collections.defaultdict(int)
            output = open('file/dump/'+naming+'.txt', "r")
            Lines = output.readlines()
            for line in Lines:
                insert = line.strip().split('\t')
                if insert and len(insert) == 2:
                    kmer = insert[0]
                    if kmer in kmers_dict: # if kmer present in library
                        read_kmer_count[kmer] = 1 # if using count data, change to id_count[id_] += int(insert[1])

            # account for missing kmers
            missing_kmers = list(set(kmers_dict.keys()) - set(read_kmer_count.keys()))
            for missing_kmer in missing_kmers:
                read_kmer_count[missing_kmer] = 0

            read_kmer_count_df = pd.DataFrame.from_dict(read_kmer_count, orient='index').T
            read_kmer_count_df.reindex(columns = kmers_list)
            ##print(read_kmer_count_df)

            
            # run model
            pred_ribo = model.predict(read_kmer_count_df)

            # get ribotype
            ribo = str(pred_ribo[0])
            message = "For file "+filename+" with k-mer size "+k+", predicted ribotype is "+ribo+'.'

        flash(message)
        return render_template('index.html') 
    return render_template('index.html') 

if __name__ == '__main__':
    app.run(debug=True)
