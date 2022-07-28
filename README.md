# coexist
Effect of using random intraspecific variability as a proxy for unknown environmental variation 

To run the analyses as the authors did, adapt the directories in make.R and adapt script.sh to your cluster, then run sbatch script.sh in a terminal. Make the "outputs" directory a subdirectory: "outputs/outputs_cluster".
Then build the figures locally with figure_conceptual.R and figures_result.R.

To run the analyses locally, adapt the make.R script and execute it:
- directory_writing <- here::here()
- directory_reading <- here::here()
- add a loop "for (s in 1:nrow(Simulations)){}
- then to build the figures with figure_conceptual.R and figures_result.R, use "outputs" and not "outputs/outputs_cluster" to charge the data.