import difflib
from pkg_resources import Requirement, resource_filename

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
    "Val":"V",
    "Ter":"X",
    "%3D":"="
}
inv_amino_acid_mapping = {v: k for k, v in amino_acid_mapping.items()}


def check_sequences(row):
    #print (row)

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
    
    if "*" in sequence:
        PASS = False
        reasons["Sequence Errors"].append("* in Sequence")

    if "U" in sequence:    
        PASS = False
        reasons["Sequence Errors"].append("U in Sequence")

    for i in range(len(mutation_positions)):
        
        if not PASS:
            pass
        elif ref_AA[i] != sequence[mutation_positions[i]-1]:    
            reasons["Mutation Errors"].append(f"{mutations[i]}: {ref_AA[i]} != {sequence[mutation_positions[i]-1]}")

    if len(reasons["Mutation Errors"]) > 0:
        PASS = False

    return PASS, reasons["Sequence Errors"], reasons["Mutation Errors"]



def clean_FASTA_sequence(sequence):

    sequence = sequence.replace("*","A")
    sequence = sequence.replace("U","A")

    if sequence.startswith("X"):
        sequence = sequence.replace("X", "", 1)

    return sequence


def mutation_mapping(x):
    #print (x)
    try:
        loc = re.findall(r'\d+', x)[0]
    except TypeError:
        print (x)
        exit()

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


def alignment_score(data):
    mutations = data["mutation"].split(" ")
    refs =      "".join([mut[0] for mut in mutations])
    positions = [int(dbnsfp_mutation_position(m))-1 for m in mutations]
    
    alignment = []
    for pos in range(len(positions)):
        try:
            if data["sequence"][positions[pos]]==refs[pos]:
                alignment.append(1)
            else:
                alignment.append(0)
        except IndexError:
            alignment.append(0)
    
    return float(sum(alignment))/float(len(positions))



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
    

def collect_dbNSFP_hgvsp(protein_id, mutation):
    #print (protein_id, mutation)

    if mutation != ".":
        return f"{protein_id}:p.{mutation_mapping(mutation)}"
    else:
        return "."



def filter_non_missense(data, file_format, annotation):
    if file_format == "dbNSFP":
        
        print (f"Pre-filtered count: {len(data)}")
        data = data[((data["aaref"].isin(amino_acid_mapping.values())) & (data["aaalt"].isin(amino_acid_mapping.values()))) | (data["aapos"] != '-1')]
        data = data[(data["aaref"] != "X") & (data["aaalt"] != "X")]
        print (f"Post-filtered count: {len(data)}")

        return data
    
    elif file_format == "VEP":

        data = data[data["Consequence"].str.contains("missense_variant")]
        data = data[~data["hgvsp"].str.contains("delins")]

        return data
    
    elif file_format=="VCF-info" and annotation=="SnpEff":

        print (f"Pre-filtered count: {len(data)}")
        data = data[data["Annotation"]=="missense_variant"]
        print (f"Post-filtered count: {len(data)}")

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
    #("/Users/bergqt01/Research/MutPredPy/resources/memory_usage.npy"))
    return np.poly1d(np.load(os.path.abspath(resource_filename(Requirement.parse("MutPredPy"), "MutPredPy/resources/memory_usage.npy")))) 
    

def time_estimate_function():
    #"/Users/bergqt01/Research/MutPredPy/resources/sequence_time.npy"))
    return np.poly1d(np.load(os.path.abspath(resource_filename(Requirement.parse("MutPredPy"), "MutPredPy/resources/sequence_time.npy"))))


def read_mutpred_input_fasta(faa_file):

    output = []

    temp = {
            "ID": "",
            "line":"",
            "mutations":"",
            "sequence": ""
        }
    
    cur_job = faa_file.split("/")[-2].split("_")[-1]
    
    with open(faa_file, "r") as fasta:
        cur_line = 1
        for line in fasta:
            
            if line.startswith(">"):
                
                if temp["ID"] != "":
                    output.append(temp.copy())

                    
                header_info = line.replace(">","").replace("\n","").split(" ")

                key = header_info[0]
                mutations = header_info[1:]

                temp["ID"] = key
                temp["mutations"] = mutations
                temp["line"] = cur_line
                cur_line += 1
                
                
            else:
                temp["sequence"] += line.replace("\n","")
                
        else:
            output.append(temp)
                

    output = pd.DataFrame(output)

    output["num_mutations"] = output["mutations"].apply(len)

    output["mutations"] = output["mutations"].apply(lambda x: ",".join(x))
    #print (output)
    output = output.drop_duplicates()

    output["job"] = cur_job

    return output


def collect_fasta(pep_file, primary="Ensembl_proteinid_v", drop_dups=True):

    output = []

    temp = {
            "Ensembl_proteinid_v": "",
            "gene_symbol":"",
            "sequence": ""
        }
    
    if not os.path.exists(pep_file):
        
        raise Exception(f"Reference FASTA file, {pep_file}, not found")
    

    with open(pep_file, "r") as fasta:
        
        for line in fasta:
            
            if line.startswith(">"):
                
                if temp["Ensembl_proteinid_v"] != "":
                    output.append(temp.copy())
                    
                header_info = line.replace(">","").split(" ")
                keys = {k.split(":")[0]:k.split(":")[1] for k in header_info if ":" in k}
                

                protein = header_info[0]
                try:
                    gene_symbol = keys["gene_symbol"]
                except KeyError:
                    gene_symbol = ""
                    
                gene = keys["gene"]
                transcript = keys["transcript"]

                temp["Ensembl_proteinid_v"] = protein
                temp["Ensembl_geneid"] = gene
                temp["Ensembl_transcriptid"] = transcript
                temp["gene_symbol"] = gene_symbol
                temp["sequence"] = ""
                
            else:
                
                temp["sequence"] += line.replace("\n","")
        else:
            output.append(temp)
                

    output = pd.DataFrame(output)
    
    if drop_dups:
        output = output.drop_duplicates(subset=[primary, "sequence"], keep="first")

    output["version"] = output["Ensembl_proteinid_v"].str.split(".").str[1]
    output["Ensembl_proteinid"] = output["Ensembl_proteinid_v"].str.split(".").str[0]
    output = output.sort_values("Ensembl_proteinid")

    output["length"] = output["sequence"].apply(lambda x: len(x))

    memory_function = memory_estimate_function()
    output["Memory Estimate (MB)"] = output["sequence"].apply(lambda x: memory_function(len(x)))

    time_function = time_estimate_function()
    output["Time per Mutation (hrs)"] = output["sequence"].apply(lambda x: time_function(len(x)))

    return output



def collect_mutations(data, file_format):

    if file_format == "dbNSFP":
        data["mutation"] = data.apply(lambda x: collect_dbNSFP_mutations(x), axis=1)

    else:
        #print (data)
        data[["Ensembl_proteinid","mutation"]] = data["hgvsp"].str.split(":",expand=True)
        data["mutation"] = data["mutation"].str.split(".").str[1].apply(lambda x: mutation_mapping(x))
    
    return data



def sequence_quality_check(data):
    

    data["sequence"] = data["sequence"].apply(lambda x: clean_FASTA_sequence(x))
    data[["status","Sequence Errors", "Mutation Errors"]] = data.apply(lambda x: pd.Series(check_sequences(x)), axis=1)
    
    return data