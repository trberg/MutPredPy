import pandas as pd
import os

from . import fasta



def write_output(data):
    print (f"{os.getcwd()}/scores/MutPred2.tsv")
    #data.to_csv("scores/MutPred2.tsv", sep="\t", index=False)



def existing_scores():

    scores = pd.read_csv("scores/MutPred2.tsv", sep="\t")
    
    return scores



def collect_scores():

    merged_df = pd.DataFrame()

    score_dir = os.listdir("intermediates/scores/")
    cur_file_num = 1
    for filename in score_dir:

        print (f"Reading {cur_file_num} of {len(score_dir)} files.", end="\r")
        
        if filename == "scores":
            mutType = 'PASS'
        else:
            mutType = filename.split(".")[1].split("_")[0]
            
        if mutType == 'missense':
            #print (filename)
            data = pd.read_csv("intermediates/scores/" + filename)

            if len(data) > 0:
                #print (data)
                ID_column_parts = data["ID"][0].split("|")

                cols = []
                
                for col in ID_column_parts:
                    if col.startswith("ENSG"):
                        cols.append("gene_transcript_id")
                    elif col.startswith("ENSP"):
                        cols.append("ensembl_protein_id")
                    elif col.startswith("ENST"):
                        cols.append("ensembl_transcript_id")
                    elif not col.startswith("ENS"):
                        cols.append("gene_symbol")
                    else:
                        cols.append("unknown")
                        print ("column name unknown")
                        exit()

                if len(data) > 0:
                    #print (data)
                    data[cols] = data["ID"].str.split("|",expand=True)
                    data["hgvsp"] = data["ensembl_protein_id"] + ":" + data["Substitution"].apply(lambda x: fasta.mutation_mapping(x))
                    
                    try:
                        data["Molecular Mechanism"] = data["Molecular mechanisms with Pr >= 0.01 and P < 0.05"]
                    except KeyError:
                        data["Molecular Mechanism"] = data["Remarks"]
                    data["Molecular Mechanism"] = data["Molecular Mechanism"].fillna("-")

                    data["Posterior Probabilities"] = 0.00
                    data["P-values"] = 1.00
                    data.fillna(".", inplace=True)

        else:
            pass
        
        merged_df =  pd.concat([data, merged_df])
        cur_file_num += 1
        
    print ("")
    key_columns = ["ensembl_protein_id","ensembl_transcript_id","gene_symbol","hgvsp"]
    for col in key_columns:
        if col not in merged_df.columns:
            merged_df[col] = "-"

    merged_df = merged_df[["ensembl_protein_id","ensembl_transcript_id","gene_symbol","hgvsp","Substitution","MutPred2 score", "Molecular Mechanism", "Posterior Probabilities", "P-values"]]

    return merged_df



def merge():

    merged_df = collect_scores()
    scores = existing_scores()
    prev_scores = len(scores)

    merged_df = pd.concat([scores, merged_df])
    merged_df = merged_df.drop_duplicates(keep="first", subset="hgvsp")

    added_scores = len(merged_df) - prev_scores

    write_output(merged_df)
    print (f"{added_scores} scored mutations added.")