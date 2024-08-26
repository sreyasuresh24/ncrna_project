CHARACTERIZATION OF LONG NON-CODING RNA IN SCHISTOSOMA MANSONI
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$
bash Miniforge3-$(uname)-$(uname -m).sh

mamba create -n ncrna --yes
mamba install --yes --file requirements.txt
mamba activate ncrna

nohup snakemake -d output_results -p -j 100 -C sample_sheet=$PWD/sample_sheet.tsv \
   --rerun-trigger mtime \
   --default-resources "mem='5G'" "time='3:00:00'" "cpus='4'" \
    --cluster 'sbatch -A project0014 --time {resources.time} --cpus-per-task={resources.cpus} --mem={resources.mem} --parsable -o slurm/{rule}.{jobid}.out -e slurm/{rule}.{jobid}.err' \
    --cluster-cancel scancel &

 ena download -p '{experiment_alias}.' ERR5727037 ERR5727030 ERR5727036 ERR5727070 ERR5727071 ERR5727073 ERR5727072 ERR5727209 ERR5727210 ERR5727226 ERR5727227 ERR5727234 ERR5727235 ERR5727236

 ena query PRJEB26892
