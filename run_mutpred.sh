#!/bin/bash
while getopts t:b:f: flag
do
    case "${flag}" in
        t) threads=${OPTARG};;
        b) base=${OPTARG};;
        f) filter=${OPTARG};;
    esac
done
echo "Threads: $threads";
echo "File Basename: $base";
echo "Filter Column: $filter";

##cd /sc/arion/projects/pejaverlab/lab_software/MutPredMerge/tools/mutpred2.0/
