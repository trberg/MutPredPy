import pandas as pd
import numpy as np
import os
import re


amino_acid_mapping = {
    "Ala":"A",
    "Arg":"R",
    "Asn":"N",
    "Asp":"D",
    "Cys":"C",
    "Glu":"E",
    "Gln":"Q",
    "Gly":"G",
    "His":"H",
    "Ile":"I",
    "Leu":"L",
    "Lys":"K",
    "Met":"M",
    "Phe":"F",
    "Pro":"P",
    "Ser":"S",
    "Sec":"U",  
    "Thr":"T",
    "Trp":"W",
    "Tyr":"Y",
    "Val":"V"
}
inv_amino_acid_mapping = {v: k for k, v in amino_acid_mapping.items()}


def check_sequences(row):

    mutations = row["mutation"].split(" ")
    mutation_positions  = [int(mut[1:-1]) for mut in mutations]
    ref_AA              = [mut[0] for mut in mutations]
    sequence            = row["sequence"]

    max_position = max(mutation_positions)

    sequence_length = len(sequence)

    PASS = True
    reasons = {
        "Sequence Errors":[],
        "Mutation Errors":[]
    }

    if sequence_length < max_position:
        PASS = False
        reasons["Sequence Errors"].append(f"Sequence {sequence_length} shorter than max position {max_position}")

    if len(sequence) < 30:
            PASS = False
            reasons["Sequence Errors"].append("Sequence < 30")

    for i in range(len(mutation_positions)):
        
        if not PASS or ref_AA[i] != sequence[mutation_positions[i]-1]:    
            PASS = False
            reasons["Mutation Errors"].append(f"{mutations[i]}")
            
        if "*" in sequence:
            PASS = False
            reasons["Mutation Errors"].append(f"{mutations[i]}")

        if "U" in sequence:    
            PASS = False
            reasons["Mutation Errors"].append(f"{mutations[i]}")

        

    return PASS, reasons["Sequence Errors"], reasons["Mutation Errors"]



def clean_FASTA_sequence(sequence):

    return sequence.replace("*","A")


def mutation_mapping(x):
    
    ref = x[:3]
    loc = x[3:-3]
    alt = x[-3:]

    loc = re.findall(r'\d+', x)[0]

    ref, alt = x.replace(str(loc),"|").split("|")
    
    
    if len(ref) == 1:
        m_ref = inv_amino_acid_mapping[ref]
    else:
        m_ref = amino_acid_mapping[ref]

    if len(alt) == 1:
        m_alt = inv_amino_acid_mapping[alt]
    else:
        m_alt = amino_acid_mapping[alt]

    return f"{m_ref}{loc}{m_alt}"



def dbnsfp_mutation_position(mut):
    pos = re.findall(r'\d+', mut)[0]
    return pos


def dbnsfp_aa_change(mut):
    pos = dbnsfp_mutation_position(mut)
    mut = mut.split(".")[-1].replace(str(pos),"|")
    return mut.split("|")


def collect_dbNSFP_mutations(x):

    if x["aaref"] != "." and x["aaalt"] != "." and x["aapos"] != ".":
        pos = x["aapos"]
        ref = x["aaref"]
        alt = x["aaalt"]
        return f"{ref}{pos}{alt}"

    elif x["HGVSp_ANNOVAR"] != ".":
        return x["HGVSp_ANNOVAR"].split(".")[-1]
    
    elif x["HGVSp_VEP"] != ".":
        pos = dbnsfp_mutation_position(x["HGVSp_VEP"])
        ref,alt = dbnsfp_aa_change(x["HGVSp_VEP"])
        return f"{ref}{pos}{alt}"
    
    else:
        return "."



def filter_non_missense(data, file_format):
    if file_format == "dbNSFP":
        
        print (f"Pre-filtered count: {len(data)}")
        data = data[((data["aaref"].isin(amino_acid_mapping.values())) & (data["aaalt"].isin(amino_acid_mapping.values()))) | (data["mutation"]!=".")]
        print (f"Post-filtered count: {len(data)}")

        return data
    
    elif file_format == "VEP":

        data = data[data["Consequence"].str.contains("missense_variant")]

        return data
    
    else:
        print (f"Pre-filtered count: {len(data)}")
        #print (data)
        muts = pd.DataFrame()
        muts[["ref","alt"]] = data["mutation"].apply(lambda x: pd.Series(dbnsfp_aa_change(x)))
        data = data[(muts["ref"].isin(amino_acid_mapping.values())) & (muts["alt"].isin(amino_acid_mapping.values()))]
        print (f"Post-filtered count: {len(data)}")

        return data
        

def memory_estimate_function():
    return np.poly1d(np.load("resources/memory_usage.npy"))

def time_estimate_function():
    return np.poly1d(np.load("resources/sequence_time.npy"))


def read_mutpred_input_fasta(faa_file):

    output = []

    temp = {
            "ID": "",
            "mutations":""
        }
    
    index = faa_file.split("/")[-1].split(".")[-2].split("_")[-1]
    
    with open(faa_file, "r") as fasta:
        
        for line in fasta:
            
            if line.startswith(">"):
                
                if temp["ID"] != "":
                    output.append(temp.copy())

                    
                header_info = line.replace(">","").replace("\n","").split(" ")

                key = header_info[0]
                mutations = header_info[1:]

                temp["ID"] = key
                temp["mutations"] = mutations
                
                
            else:
                pass
                
        else:
            output.append(temp)
                

    output = pd.DataFrame(output)

    output["num_mutations"] = output["mutations"].apply(len)

    output["index"] = index

    return output


def collect_fasta(pep_file):

    output = []

    temp = {
            "Ensembl_proteinid_v": "",
            "gene_symbol":"",
            "sequence": ""
        }
    
    with open(pep_file, "r") as fasta:
        
        for line in fasta:
            
            if line.startswith(">"):
                
                if temp["Ensembl_proteinid_v"] != "":
                    output.append(temp.copy())
                    
                header_info = line.replace(">","").split(" ")

                protein = header_info[0]
                gene_symbol = header_info[7].split(":")[1]

                temp["Ensembl_proteinid_v"] = protein
                temp["gene_symbol"] = gene_symbol
                temp["sequence"] = ""
                
            else:
                
                temp["sequence"] += line.replace("\n","")
        else:
            output.append(temp)
                

    output = pd.DataFrame(output)

    
    output["version"] = output["Ensembl_proteinid_v"].str.split(".").str[1]
    output["Ensembl_proteinid"] = output["Ensembl_proteinid_v"].str.split(".").str[0]
    output = output.sort_values("Ensembl_proteinid")

    output["length"] = output["sequence"].apply(lambda x: len(x))

    memory_function = memory_estimate_function()
    output["Memory Estimate (MB)"] = output["sequence"].apply(lambda x: memory_function(len(x)))

    time_function = time_estimate_function()
    output["Time per Mutation (hrs)"] = output["sequence"].apply(lambda x: time_function(len(x)))

    return output