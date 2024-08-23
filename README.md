# MutPredPy
A tool to manage running the MutPred2 pathogenicity prediction tool on high performance computing systems.

## Installation
```
pip install --user username /path/to/tool/directory
```

## Commands

### prepare
Takes an input VEP annotated file and outputs parallel MutPred2 jobs with LSF scripts ready to run.

```
options:
  -h, --help            Show this help message and exit
  --time                Target time in hours to run the jobs. MutPredPy uses this to calculate how many     
                        parallel jobs needs to be run at once.
  --input               Path to the input file with VEP or other variant annotations. MutPredPy currently
                        requires some form of columns with HGVSp nomenclature.
  --working_dir         Path to the directory where the outputs are written.
  --dry_run             Calculate the resources that would be used with the current parameters, but write 
                        nothing to file. Estimated resources will be printed out.
  --canonical           Only prepare canonical isoforms for mutpred scoring. Assumes VEP annotations.
  --database [DATABASE]
                        Path to optional config.yaml file for linking to mysql database. Include the name of the configuration in the config.yaml file after an @ symbol (ex. /path/to/file@Remote). If no config name included, program will
                        default to 'Local'.
```

### status

```
options:
  -h, --help            Show this help message and exit
  --job_dir             Path to the directory containing the MutPred2 jobs.
  --summary             (Optional) Path to a summary file that was previously generated by the status
                        command. Will only use the summary file and will not compute the latest status of the jobs. Will ignore the job_dir variable.
  --all                 Show the progress status of all files and jobs
  --show_incomplete     When listing the remaining jobs, show all incomplete jobs and not just jobs 
                        with zero output.
```