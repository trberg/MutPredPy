
from string import Template


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
    normal_job     = tech_requirements[tech_requirements['Normal Memory']]

    #print (high_memory)
    #print (normal_job)

    if len(high_memory) > 0:
        high_mem_usage = int((max(high_memory["Memory Minimum"]) + memory_cushion)) * len(high_memory)
        high_time_usage = int(max(high_memory["Time Estimate"])) + time_cushion
    else:
        high_mem_usage = 0
        high_time_usage = 0

    
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
        "hmu":int(high_mem_usage/1000),
        "ht":high_time_usage,
        "hj":len(high_memory),
        "tmu":int((high_mem_usage+norm_mem_usage)/1000),
        "tt":max([norm_time_usage,high_time_usage]),
        "tj":len(normal_job)+len(high_memory)
    })

    print (report)


def config_template():
    return Template("""
#BSUB -P acc_pejaverlab
#BSUB -n 4
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=$mem]
#BSUB -W $time
#BSUB -q premium
#BSUB -J $job$job_array
#BSUB -oo /hpc/users/bergqt01/pejaverlab/lab_software/MutPredPy/logs/$project/out_$base.%J.faa_file_%I
#BSUB -e /hpc/users/bergqt01/pejaverlab/lab_software/MutPredPy/logs/$project/err_$base.%J.faa_file_%I

/sc/arion/projects/pejaverlab/IGVF/src/mutpred2_dev \
-i /hpc/users/bergqt01/pejaverlab/lab_software/MutPredPy/$intermediate_dir/faa/$project/$base.missense_$index.faa \
-o /hpc/users/bergqt01/pejaverlab/lab_software/MutPredPy/$intermediate_dir/scores/$base.missense_output_$index.txt \
-p 1 -c 1 -b 0 -t 1 -f 2 \
-d /sc/arion/projects/pejaverlab/IGVF/data/mutpred2.0/
""")

"""
./run_mutpred2.sh -i $intermediate_dir/faa/$project/$base.missense_$index.faa -p 1 -c 1 -b 0 -t 0.05 -f 2 -o $intermediate_dir/scores/$base.missense_output_$index.txt
"""

def build_lsf_config_file(tech_requirements, intermediate_dir, project, base):

    

    high_memory = tech_requirements[tech_requirements['High Memory']]
    normal_job  = tech_requirements[tech_requirements['Normal Memory']]

    jobs = [normal_job, high_memory]

    for job in jobs:
        if len(job) > 0:
            template = config_template().substitute({
                'mem': int((max(job["Memory Minimum"]) + memory_cushion)/cores),
                'time': f'{int(max(job["Time Estimate"])) + time_cushion}:00',
                'job': f"{project}_variants",
                'job_array': build_job_array(job['File']),
                'project': project,
                'intermediate_dir': intermediate_dir,
                'base': base,
                'index': "$LSB_JOBINDEX"
            })
            print (template)
    

