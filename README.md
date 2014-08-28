### How this works

Here is the general process. You are assumed to be using the ACI cluster.

0. Compile the three executables in the 3D FT folder. (2 minutes)
..* `cd 3dft; make ift -C src; make ft -C src; make stdev -C src`

1. Create a new directory for each model you are running. Generate a 512 pixel FT of your model in this directory. This allows you to visualize your structure in reciprocal space. Modify the slurm.sh file to submit the 3dft executable. This should be a simple comment/uncomment. (up to a few hours)
..* `sbatch path/slurm.sh modelfile outbase 512`

2. Run the Igor spot analysis [here](https://github.com/paul-voyles/Igor/blob/master/3D%20FFT%20Analysis.ipf). Follow the instructions in the header. This analyzes the high intensity diffraction spots in your FT and creates a parameter file that you will use to extract and quantify important parts of your model. (20-60 minutes)

3. Copy the parameter file to the directory you created in (1). Run this through the IFT program by submitting a job to the cluster again. Modify the slurm.sh file to submit the ift executable. This should be a simple comment/uncomment. (~6-24 hours)
..* `sbatch path/slurm.sh paramfile`

4. Run the batch_convert.py program in an "int" shell. (< 30 minutes)
..* `srun -n16 -N1 -p int --pty bash`
..* `python path/batch_convert.py paramfile jobid 512`

5. When this finishes you will have model files (I call these sub-models) for each spot identified (up to 20 spots). Now you will generate FTs for each sub-model. These are used to confirm you found the correct structure. (3-6 minutes per job)
..* `sbatch path/spot_ft.py paramfile jobid 512`

6. Now you analyze the sub-models using a program in the [model_analysis repo](https://github.com/refreshx2/model_analysis). First you need to create spot_fts/ and submodels/ directories within each of your main folders and up the FTs from (5) in the spot_fts/ directory and the sub-models from (4-5) in the submodels/ directory. Then run the analysis program by submitting it to the cluster. (~1 hour)
..* `mkdir submodels; mv *.xyz *.cif submodels/; mkdir spot_fts/; mv *_512_ft.gfx spot_fts/`
..* `sbatch ~/model_analysis/scripts/submit_ift_cluster_analysis.sh paramfile jobid original_ft main_ft_direc VP_categories_paramfile`
..* For example: `sbatch ~/model_analysis/scripts/submit_ift_cluster_analysis.sh ~/3dft/t1updated/paramfile_t1_update.txt 65421 ~/3dft/t1_512_ft.gfx ~/3dft/t1updated/ ~/model_analysis/scripts/categorize_parameters_updated.txt`

7. Copy the last lines of the output file into e.g. Excel and analyze from there. I wont explain what everything means here. You will either have to read the ift_cluster_analysis.py file or talk to me.
