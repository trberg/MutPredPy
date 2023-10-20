
from string import Template
import os
import pandas as pd


## GLOBAL VARIABLES
memory_cushion = 2000
time_cushion = 6
cores = 4


def build_job_array(jobs):

    jobs = sorted(list(jobs))

    job_sequences = []
    previous_job = -10
    for j in jobs:
        if previous_job < 0:
            cur_index = 0
            job_sequences.append([j])
            
        elif j == previous_job+1:
            job_sequences[cur_index].append(j)
            
        else:
            cur_index += 1
            job_sequences.append([j])
        
        previous_job = j

    job_sequence = f'[{",".join([f"{min(s)}-{max(s)}" if min(s) != max(s) else str(max(s)) for s in job_sequences])}]'
    
    return job_sequence


def report_template():
    return Template(
"""========== Normal Jobs ==========
    Jobs            $nj
    Memory Usage    $nmu Gb
    Time            $nt Hrs

========== Middle Jobs ==========
    Jobs            $mj
    Memory Usage    $mmu Gb
    Time            $mt Hrs
                    
========== High Memory Jobs ==========
    Jobs            $hj
    Memory Usage    $hmu Gb
    Time            $ht Hrs
                    
========== Total Usage ==========
    Total Jobs      $tj
    Memory Usage    $tmu Gb
    Time            $tt Hrs
""")


def usage_report(tech_requirements):
    
    #print (tech_requirements)

    high_memory = tech_requirements[tech_requirements['High Memory']]
    mid_memory  = tech_requirements[tech_requirements['Middle Memory']]
    normal_job  = tech_requirements[tech_requirements['Normal Memory']]

    #print ("============= LOW MEMORY =============")
    #print (normal_job)

    #print ("============= MID MEMORY =============")
    #print (mid_memory)
    
    #print ("============= HIGH MEMORY =============")
    #print (high_memory)


    if len(high_memory) > 0:
        high_mem_usage = int((max(high_memory["Memory Minimum"]) + memory_cushion)) * len(high_memory)
        high_time_usage = int(max(high_memory["Time Estimate"])) + time_cushion
    else:
        high_mem_usage = 0
        high_time_usage = 0

    
    if len(mid_memory) > 0:
        mid_mem_usage = int((max(mid_memory["Memory Minimum"]) + memory_cushion)) * len(mid_memory)
        mid_time_usage = int(max(mid_memory["Time Estimate"])) + time_cushion
    else:
        mid_mem_usage = 0
        mid_time_usage = 0

    
    if len(normal_job) > 0:
        norm_mem_usage = int((max(normal_job["Memory Minimum"]) + memory_cushion)) * len(normal_job)
        norm_time_usage = int(max(normal_job["Time Estimate"])) + time_cushion
    else:
        norm_mem_usage = 0
        norm_time_usage = 0
    

    report = report_template().substitute({
        "nmu":int(norm_mem_usage/1000),
        "nt":norm_time_usage,
        "nj":len(normal_job),

        "mmu":int(mid_mem_usage/1000),
        "mt":mid_time_usage,
        "mj":len(mid_memory),

        "hmu":int(high_mem_usage/1000),
        "ht":high_time_usage,
        "hj":len(high_memory),

        "tmu":int((high_mem_usage+norm_mem_usage+mid_mem_usage)/1000),
        "tt":max([norm_time_usage,high_time_usage,mid_mem_usage]),
        "tj":len(normal_job)+len(high_memory) + len(mid_memory)
    })

    print (report)


def split_for_multiple_users(tech_requirements, users):
    splits = users
    
    high_memory = tech_requirements[tech_requirements['High Memory']]
    mid_memory  = tech_requirements[tech_requirements['Middle Memory']]
    normal_job  = tech_requirements[tech_requirements['Normal Memory']]

    high_mem_splits = int(len(high_memory)/splits)
    mid_mem_splits  = int(len(mid_memory)/splits)
    norm_mem_splits = int(len(normal_job)/splits)

    user_outputs = []
    for user in range(users):
        
        high_mem_cur_start = user*high_mem_splits
        mid_mem_cur_start  = user*mid_mem_splits
        norm_mem_cur_start = user*norm_mem_splits

        high_mem_cur_end = (user+1)*high_mem_splits
        mid_mem_cur_end  = (user+1)*mid_mem_splits
        norm_mem_cur_end = (user+1)*norm_mem_splits

        user_high_memory = high_memory.iloc[high_mem_cur_start:high_mem_cur_end]
        user_mid_memory = mid_memory.iloc[mid_mem_cur_start:mid_mem_cur_end]
        user_normal_memory = normal_job.iloc[norm_mem_cur_start:norm_mem_cur_end]

        user_jobs = pd.concat([user_normal_memory, user_high_memory, user_mid_memory])

        user_outputs.append(user_jobs)

    return user_outputs



def config_template():
    return Template("""
#BSUB -P acc_pejaverlab
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=$mem]
#BSUB -W $time
#BSUB -q premium
#BSUB -J $job$job_array
#BSUB -cwd /sc/arion/projects/pejaverlab/lab_software/MutPredPy/tools/mutpred2.0
#BSUB -oo /sc/arion/projects/pejaverlab/lab_software/MutPredPy/logs/$project/out_$base.%J.faa_file_%I
#BSUB -e /sc/arion/projects/pejaverlab/lab_software/MutPredPy/logs/$project/err_$base.%J.faa_file_%I

module load MCR/R2017b;

/sc/arion/projects/pejaverlab/IGVF/src/mutpred2_dev \
-i /sc/arion/projects/pejaverlab/lab_software/MutPredPy/$intermediate_dir/faa/$project/$base.missense_$index.faa \
-o /sc/arion/projects/pejaverlab/lab_software/MutPredPy/$intermediate_dir/scores/$project/$index/$base.missense_output_$index.txt \
-d /sc/arion/projects/pejaverlab/IGVF/data/mutpred2.0/ \
-p 1 -c 1 -b 0 -t 1 -f 2 \

""")

"""
./run_mutpred2.sh
/sc/arion/projects/pejaverlab/IGVF/src/mutpred2_dev \
-i $intermediate_dir/faa/$project/$base.missense_$index.faa \
-o $intermediate_dir/scores/$base.missense_output_$index.txt \
-d /sc/arion/projects/pejaverlab/IGVF/data/mutpred2.0/ \
-p 1 -c 1 -b 0 -t 1 -f 2 

"""

def build_lsf_config_file(tech_requirements, intermediate_dir, project, base, user, dry_run):

    user += 1

    high_memory = tech_requirements[tech_requirements['High Memory']]
    mid_memory  = tech_requirements[tech_requirements['Middle Memory']]
    normal_job  = tech_requirements[tech_requirements['Normal Memory']]

    jobs = [normal_job, mid_memory, high_memory]
    job_type = ["normal", "middle", "high"]

    for i in range(len(jobs)):
    #for job in jobs:
        if len(jobs[i]) > 0:
            template = config_template().substitute({
                'mem': int((max(jobs[i]["Memory Minimum"]) + memory_cushion)/cores),
                'time': "144:00",#f'{int(max(jobs[i]["Time Estimate"])) + time_cushion}:00',
                'job': f"{project}_variants",
                'job_array': build_job_array(jobs[i]['File']),
                'project': project,
                'intermediate_dir': intermediate_dir,
                'base': base,
                'index': "$LSB_JOBINDEX"
            })
            #print (template)

            if not os.path.exists("scripts"):
                if dry_run:
                    print ("No 'scripts' folder.")
                else:
                    os.mkdir("scripts")
            if job_type[i] == "normal":
                output_config_file_name = f"{project}_{user}.lsf"
            elif job_type[i] == "middle":
                output_config_file_name = f"{project}_{user}_mid_mem.lsf"
            elif job_type[i] == "high":
                output_config_file_name = f"{project}_{user}_high_mem.lsf"

            print (template)
            print (f"Written to scripts/{output_config_file_name}")
            
            if dry_run:
                pass
            else:
                with open(f"scripts/{output_config_file_name}","w") as s:
                    s.write(template)
