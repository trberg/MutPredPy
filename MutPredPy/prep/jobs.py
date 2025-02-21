import pandas as pd


def clean_splits(row, col_mapping):
    """
    Adds two new columns to the row: 'num_mutations' and 'Time Estimate (hrs)'.
    - 'num_mutations': Counts the number of mutations in the 'mutation_column'.
    - 'Time Estimate (hrs)': Calculates the total time estimate based on the number of mutations.
    
    Args:
        row (pd.Series): A row from the DataFrame.
        col_mapping (dict): A dictionary mapping column names.
    
    Returns:
        pd.Series: The updated row with additional columns.
    """
    row["num_mutations"] = len(row[col_mapping.get("mutation_column")].split(" "))
    row["Time Estimate (hrs)"] = row["Time per Mutation (hrs)"]*row["num_mutations"]
    
    return row



def split_mutations(x, threshold, num_overflow_mutations, col_mapping):
    """
    Splits mutations into chunks based on the time threshold.
    - If there are overflow mutations, they are kept in the first split.
    - Subsequent splits are evenly divided to ensure each split stays below the time threshold.
    
    Args:
        x (pd.Series): A row from the DataFrame.
        threshold (float): Time threshold for splitting.
        num_overflow_mutations (int): Number of mutations that overflow the threshold.
        col_mapping (dict): Dictionary mapping column names.
    
    Returns:
        list: A list of mutation chunks as strings.
    """
    mutation_column = col_mapping.get("mutation_column")
    mutations = x[mutation_column].split(" ")

    # Separate overflow mutations and leftover mutations
    overflow_mutations = mutations[0:num_overflow_mutations]
    leftover_mutations = mutations[num_overflow_mutations:len(mutations)]

    # Calculate the size of each split for leftover mutations
    size_of_splits = int(len(leftover_mutations)/((x["Time Estimate (hrs)"] - (num_overflow_mutations * x["Time per Mutation (hrs)"]))/threshold))
    
    # Return mutation chunks as a list of strings
    if len(overflow_mutations) > 0:
        return [" ".join(overflow_mutations)] + [" ".join(leftover_mutations[i:i+size_of_splits]) for i in range(0,len(mutations[num_overflow_mutations:len(mutations)]),size_of_splits)]
    else:
        return [" ".join(leftover_mutations[i:i+size_of_splits]) for i in range(0,len(mutations[num_overflow_mutations:len(mutations)]),size_of_splits)]



def split_sequence(data, cur_number, threshold, col_mapping):
    """
    Splits a sequence into smaller chunks if the time estimate exceeds the threshold.
    
    Args:
        data (pd.Series): A single row representing the sequence data.
        cur_number (float): Current cumulative time estimate.
        threshold (float): Time threshold for splitting.
        col_mapping (dict): Dictionary mapping column names.
    
    Returns:
        pd.DataFrame: A DataFrame with split sequences.
    """
    mutation_column = col_mapping.get("mutation_column")
    overflow = max(0, cur_number - threshold)
    num_overflow_mutations = int(overflow / data["Time per Mutation (hrs)"])

    # Convert the row to a DataFrame and split the mutations
    data = data.to_frame().T
    data[mutation_column] = data.apply(lambda x: split_mutations(x, threshold, num_overflow_mutations, col_mapping), axis=1)
    
    # Explode the mutation column to create multiple rows for each split
    data = data.explode(mutation_column)
    data = data.apply(lambda row: clean_splits(row, col_mapping), axis=1)
    data.reset_index(drop=True, inplace=True)
    
    return data
        


def split_data(prepare, data, col_mapping, file_number=1):
    """
    Splits large sequences into smaller chunks and writes them to files.
    Tracks job information such as time estimates and memory usage for each split.
    
    Args:
        prepare (object): An instance of the Prepare class.
        data (pd.DataFrame): The input data containing sequences and mutations.
        col_mapping (dict): Dictionary mapping column names.
        file_number (int): Starting file number for output files.
    
    Returns:
        pd.DataFrame: Job information with time estimates and memory usage for each split.
    """
    job_time = 0
    job_information = []
    memory = []
    threshold = prepare.get_time()
    time_est  = "Time Estimate (hrs)"
    
    def record_job_info(file_number, time_estimate, memory_list):
        """Records job information for a given file."""
        job_information.append({
            "File": file_number,
            "Time Estimate": time_estimate,
            "Memory Minimum": max(memory_list)
        })
    
    def write_sequence(row, file_number, time_estimate, location):

        """Writes a sequence to a file with a given file number and time estimate."""
        header = f">{'|'.join(list(row[col_mapping['id_column']]))} {row[col_mapping['mutation_column']]}\n"
        sequence = f"{row['sequence']}\n"
        prepare.write_sequence_to_file(time_estimate, file_number, header, sequence)

    for _,row in data.iterrows():
        if job_time + row[time_est] > threshold + 2:

            # Split the sequence if it exceeds the threshold
            splits = split_sequence(row, job_time, threshold, col_mapping)

            for split_index, split in splits.iterrows():
                write_sequence(split, file_number, split[time_est], "Splits loop")
                memory.append(split["Memory Estimate (MB)"])
                job_time = split[time_est]
                
                # Record job info after each split except the last one
                if split_index < len(splits) - 1:
                    record_job_info(file_number, job_time, memory)
                    memory = []
                    file_number += 1
                    
        else:
            # Write the full sequence if it doesn't exceed the threshold
            write_sequence(row, file_number, job_time, "ELSE")
            memory.append(row["Memory Estimate (MB)"])
            job_time += row[time_est]
        
        # Record job info when the threshold is reached
        if job_time >= threshold:
            record_job_info(file_number, job_time, memory)
            job_time = 0
            memory = []
            file_number += 1
            
    # Final job info recording
    #if memory:
    #    record_job_info(file_number, job_time, memory)
    #    write_sequence(row, file_number, job_time, "FINAL JOB")

    # Convert job information to a DataFrame and categorize memory usage
    job_information = pd.DataFrame(job_information)
    job_information["Normal Memory"] = job_information["Memory Minimum"].apply(lambda x: x <= 5000)
    job_information["Middle Memory"] = job_information["Memory Minimum"].apply(lambda x: x > 5000 and x <= 10000)
    job_information["High Memory"]   = job_information["Memory Minimum"].apply(lambda x: x > 10000)

    # Set Time Estimates and Memory Minimum rounding
    job_information["Time Estimate"] = job_information["Time Estimate"].apply(lambda x: round(x, 2))
    job_information["Memory Minimum"] = job_information["Memory Minimum"].apply(lambda x: round(x, 2))

    return job_information