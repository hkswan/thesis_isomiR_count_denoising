### IMPORTANT FOR UNDERSTANDING HOW FILES ARE ORGANIZED 
<\br>
<\br>
Results of various simulation studies for analyzing the performance of my isomiR count denoising algorithm are stored in the sims directory. In the sims directory is a subdirectory labeled with the date that the simulations were run. For example, the folder contains results of denoising simulated length variant count data that were run on January 8th, 2025. The convention for writing the date for the folder name is MM_DD_YYYY. The date folder then contains two more folders: `partition_objs` and `results`. 
</br>
#### partition_objs folder 
<\br> 
The `partition_objs` folder contains .Rds files. Each .Rds file contains the output of running the `denoise_isomiR_counts` function for a specific miRNA on the simulated data. The name of the miRNA that was passed as an argument function is contained in the name of the file, beginning with "Mmu-". Loading this object will give you things like the `partition_df` object, the maximum number of iterations run, etc. 

<\br>
#### the results folder 
<\br> 
The `results` folder also contains .Rds files. Each .Rds file contains a list of vectors of performance metrics calculated on the output of running the `denoise_isomiR_counts` function for a given miRNA on each of the simulated datasets. The name of the given miRNA that was passed as an argument to the `denoise_isomiR_counts` function is given in the name of the file, beginning with "Mmu-". Loading this object will give you things like the total number of partitions created, the number of false positives from each dataset, the total number of sequences in each dataset, etc. The length of each vector element of the list object tells you the total number of datasets that were simulated e.g. a vector of length 100 means that that we simulated 100 datasets. 
