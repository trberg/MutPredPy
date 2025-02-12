
from string import Template
import os
import pandas as pd

# ==========================
# 🚀 GLOBAL CONFIGURATIONS
# ==========================
MEMORY_CUSHION = 2000  # Additional memory buffer
TIME_CUSHION = 6  # Additional time buffer (hours)
CORES = 4  # Number of CPU cores
TEMPLATES_DIR = "templates"  # Directory for storing templates


# ==========================
# 📝 TEMPLATE LOADING
# ==========================
def load_template(template_name):
    """
    Loads a template file from the 'templates/' directory.
    
    Args:
        template_name (str): The template filename.
    
    Returns:
        Template: A string Template object.

    Raises:
        FileNotFoundError: If the template file is missing.
        IOError: If there is an issue reading the file.
        ValueError: If the template is empty.
    """
    template_path = os.path.join(TEMPLATES_DIR, template_name)

    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template file '{template_path}' not found in '{TEMPLATES_DIR}'.")

    try:
        with open(template_path, "r") as file:
            content = file.read().strip()
    except IOError as e:
        raise IOError(f"Error reading template file '{template_path}': {e}")

    if not content:
        raise ValueError(f"Template file '{template_path}' is empty.")

    return Template(content)


# ==========================
# 🔧 UTILITY FUNCTIONS
# ==========================
def build_job_array(jobs):
    """
    Converts a list of job indices into an LSF job array format.
    
    Args:
        jobs (list): List of job indices.
    
    Returns:
        str: Job array string in LSF format.
    """
    jobs = sorted(jobs)
    sequences = []
    prev_job = -10

    for j in jobs:
        if prev_job < 0 or j != prev_job + 1:
            sequences.append([j])
        else:
            sequences[-1].append(j)
        prev_job = j

    return f'[{",".join(f"{min(s)}-{max(s)}" if len(s) > 1 else str(s[0]) for s in sequences)}]'


# ==========================
# 📊 USAGE REPORT
# ==========================
def usage_report(tech_requirements):
    """
    Generates and prints a usage report based on memory and time estimates.
    
    Args:
        tech_requirements (DataFrame): Dataframe containing memory and time requirements.
    """
    # Categorize jobs
    high_memory = tech_requirements[tech_requirements["High Memory"]]
    mid_memory = tech_requirements[tech_requirements["Middle Memory"]]
    normal_jobs = tech_requirements[tech_requirements["Normal Memory"]]

    # Compute memory and time usage for each category
    def compute_usage(jobs):
        if len(jobs) == 0:
            return 0, 0
        return (max(jobs["Memory Minimum"]) + MEMORY_CUSHION) * len(jobs), max(jobs["Time Estimate"]) + TIME_CUSHION

    high_mem_usage, high_time_usage = compute_usage(high_memory)
    mid_mem_usage, mid_time_usage = compute_usage(mid_memory)
    norm_mem_usage, norm_time_usage = compute_usage(normal_jobs)

    # Generate report using template
    report = load_template("report_template.txt").substitute({
        "nmu": norm_mem_usage // 1000, "nt": norm_time_usage, "nj": len(normal_jobs),
        "mmu": mid_mem_usage // 1000, "mt": mid_time_usage, "mj": len(mid_memory),
        "hmu": high_mem_usage // 1000, "ht": high_time_usage, "hj": len(high_memory),
        "tmu": (high_mem_usage + mid_mem_usage + norm_mem_usage) // 1000,
        "tt": max(norm_time_usage, high_time_usage, mid_time_usage),
        "tj": len(normal_jobs) + len(high_memory) + len(mid_memory)
    })

    print(report)


# ==========================
# 👥 MULTI-USER JOB SPLITTING
# ==========================
def split_for_multiple_users(tech_requirements, users):
    """
    Splits jobs among multiple users for parallel execution.

    Args:
        tech_requirements (DataFrame): Dataframe containing job details.
        users (int): Number of users to split jobs among.

    Returns:
        list: List of DataFrames, each containing jobs for a user.
    """
    categories = ["High Memory", "Middle Memory", "Normal Memory"]
    user_outputs = []

    for user in range(users):
        user_jobs = pd.concat([
            tech_requirements[tech_requirements[cat]].iloc[user::users] for cat in categories
        ])
        user_outputs.append(user_jobs)

    return user_outputs


# ==========================
# 🏗️ LSF CONFIG FILE GENERATION
# ==========================
def build_lsf_config_file(tech_requirements, working_dir, base, user, dry_run):
    """
    Generates LSF configuration files based on job requirements.

    Args:
        tech_requirements (DataFrame): Job requirements.
        working_dir (str): Working directory path.
        base (str): Base name for the job.
        user (int): User index.
        dry_run (bool): If True, only prints the output instead of writing to files.
    """
    user += 1  # Adjust user index

    job_types = [("normal", "Normal Memory"), ("middle", "Middle Memory"), ("high", "High Memory")]

    for job_type, category in job_types:
        jobs = tech_requirements[tech_requirements[category]]
        if len(jobs) == 0:
            continue

        # Load LSF template and substitute values
        config_content = load_template("lsf_configuration_template.txt").substitute({
            "mem": int((max(jobs["Memory Minimum"]) + MEMORY_CUSHION) / CORES),
            "time": f"{int(max(jobs['Time Estimate'])) + TIME_CUSHION}:00",
            "job": f"{base}_variants",
            "job_array": build_job_array(jobs["File"]),
            "working_dir": working_dir,
            "base": base,
            "index": "$LSB_JOBINDEX"
        })

        # Ensure scripts directory exists
        scripts_dir = os.path.join(working_dir, "scripts")
        if not os.path.exists(scripts_dir):
            if dry_run:
                print(f"No 'scripts' folder in {working_dir}.")
            else:
                os.makedirs(scripts_dir)

        # Define output file name
        output_file = os.path.join(scripts_dir, f"{base}_{user}_{job_type}.lsf")

        print(config_content)
        print(f"Written to {output_file}")

        # Write to file if not in dry run mode
        if not dry_run:
            with open(output_file, "w") as file:
                file.write(config_content)
