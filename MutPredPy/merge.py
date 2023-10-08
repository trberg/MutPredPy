import pandas as pd
import os
import re

from . import fasta
#import fasta


def write_output(data, dry_run):
    if dry_run:
        print (f"(Dry Run) Would write scores to > {os.getcwd()}/scores/MutPred2.tsv")
    else:
        print (f"Writing scores to > {os.getcwd()}/scores/MutPred2.tsv")
        data.to_csv("scores/MutPred2.tsv", sep="\t", index=False)



def existing_scores():
    if os.path.exists("scores/MutPred2.tsv"):   
        scores = pd.read_csv("scores/MutPred2.tsv", sep="\t")
    else:
        scores = pd.DataFrame()
    
    return scores



def collect_scores():

    projects = [project for project in os.listdir("intermediates/scores/")]
    score_dir = []
    for project in projects:
        score_dir += os.listdir(f"intermediates/scores/{project}")

    
    #score_dir = os.listdir("intermediates/scores/")
    score_pattern = re.compile(".*\.missense_output_\d+.txt$")

    scores = [s for s in score_dir if score_pattern.match(s)]

    cur_file_num = 1
    score_df = []

    for filename in scores:

        print (f"Reading {cur_file_num} of {len(scores)} files.", end="\r")
        
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
                    if "gene_symbol" in cols:
                        data = data[data["gene_symbol"]!="None"]
                    try:
                        data["hgvsp"] = data["ensembl_protein_id"] + ":p." + data["Substitution"].apply(lambda x: fasta.mutation_mapping(x))
                    except:
                        print (filename)
                        print (data)
                        exit()
                    
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
        
        score_df.append(data)
        cur_file_num += 1

    print ("")
    merged_df =  pd.concat(score_df)
    print ("Outputs merged")
    
    key_columns = ["ensembl_protein_id","ensembl_transcript_id","gene_symbol","hgvsp"]
    for col in key_columns:
        if col not in merged_df.columns:
            merged_df[col] = "-"

    merged_df = merged_df[["ensembl_protein_id","ensembl_transcript_id","gene_symbol","hgvsp","Substitution","MutPred2 score", "Molecular Mechanism", "Posterior Probabilities", "P-values"]]

    return merged_df



def merge(dry_run):

    merged_df = collect_scores()
    scores = existing_scores()
    prev_scores = len(scores)

    merged_df = pd.concat([scores, merged_df])
    merged_df = merged_df.drop_duplicates(keep="first", subset="hgvsp")

    added_scores = len(merged_df) - prev_scores

    write_output(merged_df, dry_run)
    print (f"{added_scores} scored mutations added.")



if __name__=="__main__":

    id_mappings = fasta.collect_fasta("/Users/bergqt01/Research/MutPredPy/resources/Homo_sapiens.GRCh38.combined.pep.all.fa",primary="Ensembl_transcriptid")
    id_mappings = id_mappings[["Ensembl_proteinid_v","Ensembl_transcriptid"]]
    id_mappings.rename(columns={"Ensembl_proteinid_v":"Ensembl_proteinid"}, inplace=True)


    def collect_value(x, val):
        cur_dict = {k.split("=")[0]:k.split("=")[1] for k in x.split(";")}
        
        if val in cur_dict.keys():
            protein = cur_dict[val]#.split(":")[0]
        else:
            protein = "."
        return protein
    

    f = "/Users/bergqt01/Research/MutPredPy/inputs/nilah/gnomad.final_set.vcf"

    header_rows = []
    header_info = ""
    for i, line in enumerate(open(f)):
        if line.startswith('##'):
            header_rows.append(i)
            header_info += line

    columns = []
    for line in header_info.split("\n"):
        if "ID=ANN" in line:
            columns = line.split("'")[1].split(" | ")
    #print (columns)

    inputs = pd.read_csv(f, sep="\t", skiprows = [i for i, line in enumerate(open(f)) if line.startswith('##')])
    inputs["INFO"] = inputs["INFO"].apply(lambda x: collect_value(x, "ANN"))
    #print (inputs)
    #exit()
    inputs["INFO"] = inputs["INFO"].str.split(",") 
    inputs = inputs.explode("INFO")
    inputs[columns] = inputs["INFO"].str.split("|", expand=True)
    
    inputs = inputs[['#CHROM','POS','ID','REF','ALT','Annotation', 'Annotation_Impact', 'Feature_ID', 'Gene_Name', 'Transcript_BioType', 'HGVS.p']]
    inputs = inputs[inputs["Annotation"]=="missense_variant"]
    inputs["Substitution"] = inputs["HGVS.p"].apply(lambda x: fasta.mutation_mapping(x.split(".")[-1]))
    inputs = inputs.merge(id_mappings, left_on="Feature_ID", right_on="Ensembl_transcriptid")
    print (inputs)

    scores = pd.read_csv("scores/MutPred2.tsv", sep="\t")
    print (scores)

    scored = inputs.merge(scores, left_on=["Ensembl_proteinid","Substitution"], right_on=["ensembl_protein_id","Substitution"], how="left")
    scored = scored[['#CHROM','POS','ID','REF','ALT','Ensembl_proteinid', 'Feature_ID','Substitution','MutPred2 score']]

    scored["MutPred2 score"].fillna('.', inplace=True)
    scored.rename(columns={"Feature_ID":"Ensembl_transcriptid"}, inplace=True)
    scored["Location"] = scored["#CHROM"]+":"+scored["POS"].astype(str)
    print (scored[scored.duplicated(subset=["Location"], keep=False)].sort_values(["#CHROM","POS"]))
    #print (scored[scored["MutPred2 score"]=="."])

    #scored.to_csv("/Users/bergqt01/Research/MutPredPy/inputs/nilah/gnomad.final_set.scored.txt", sep="\t", index=False)