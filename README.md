# Introduction
- Packages used and version of R are detailed in [Section 1](#s1).
- [Section 2] contains a guide to fully reproduce the results of simulations and analysis on real data. This may require <b>several days of running</b> on a desktop PC. Alternatively, results (i.e., the MCMC chains and simulated data) can be downloaded following procedure described [here](#s2-ss3).
- Steps described in [Section 3](#s3) reproduce all graphs included in the article.
- All the  files in online shared folder are described in [Section 4](#s5).

# Section 1:   

The required packages are listed below, all are available on `CRAN` and can be easly installed using `install.packages("name-of-package")`.

### Pakages required

```{r c00, eval=TRUE, message=FALSE, warning=FALSE, echo=TRUE}
library(gridExtra)
library(scales)
library(readr)
library(readxl)
library(stringr)
library(cowplot)
library(reshape2)
library(BayesLogit)
library(meta)
library(mvmeta)
library(gdata)
library(varhandle)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(doSNOW)
library(MASS)
library(HDInterval)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggridges)
library(latex2exp)
library(ggforce)
library(mvtnorm)
library(survival)
library(metafor)
library(mixmeta)

```

```{r c0, echo=TRUE, message=FALSE, warning=FALSE}
sessionInfo()
```


# Section 2: 

This section describes the steps to fully reproduce the results of the study, using the original code.
This can be time-consuming.
To reproduce the plots, the intermedie results can also be downloaded following the steps described [here](#s2-ss3).

The steps described in this Section require the installation of packages listed in [Section 1](#s1).
Computation time depends on the cores available for parallel computation. 
By default the script will use the number of detected logical cores minus 2.

### Reproduce MCMC chains simulations {#s2-ss1}
Make sure the  `Scripts` and `Dati` folders and the other files are in the same folder. 
After that, manually (using source button) or via console (using your own correct path) run `source_file.R`.
This file creates the folder needed to systematically save the chains and will set the working directory as the folder where files are stored. E.g:  
```{r c1_1, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("G:/mvPTGP/folder_to_push/source_file.R", echo=TRUE)
```
After that, you can run the simulation script. 
```{r c1_2, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("A_Polya_Tree_Meta_Analysis_simulation_seed_1_50.R")
```

### Reproduce MCMC chain for real data {#s2-ss2}
Running both the data and simulations scripts is not required.
As in previous section, make sure the `Scripts` and `Dati` folders and the other files are stored in the same folder. 
Then, simply run the script `source_file.R` or set the working directory in the folder using `setwd("my-path")`.
E.g.:
```{r c1_3, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
setwd("G:/mvPTGP/folder_to_push")
```
After that run the script `A_Polya_Tree_PFS.R`.
```{r c1_4, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("A_Polya_Tree_PFS.R")
```

### Download results       {#s2-ss3}
Intermediate results can be downloaded at this [link](https://drive.google.com/drive/folders/11LrUmWVdpX2sUCE354xLjKVNywFZNtXR?usp=sharing){target="_blank"}.
You can download the full folder.
Alternatively, you can download  compressed folder that store all the necessary files (i.e. `PTGP_simulation.zip`) and target `.rds` objects for target checks (e.g., `Model_02112023.rds`).

In any case, make sure to follow the correct folder structure created following the structure of the shared folder.

To reproduce Figures follows step in [Section 3](#s3).

# Section 3: {#s3}

After each one of these steps, the desired graphs will be saved in `.PDF` format in the R working directory.

### Simulations results with chains postprocessing (aprox. 2 hours)
Make sure that:

* Folder `PTGP_simulation` contains all 50 chains.
* The structure of those folders recalls the structure of the [shared folder](https://drive.google.com/drive/folders/11LrUmWVdpX2sUCE354xLjKVNywFZNtXR?usp=sharing){target="_blank"}.
* All packages are properly installed.

Now simply run manually (via source button) or via console (using your own corect path) the `A_Polya_Tree_Meta_Analysis_simulation_post.R` script. E.g.:

```{r c1_52, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("A_Polya_Tree_Meta_Analysis_simulation_post.R")
```

This procedure iterates  the 50 chains and  2500 samples that we stored for each chain.
It engages a standard PC for a few hours.

### Simulations results with stored results
Make sure that:

* Folder `RESULTS_GPPT_SIM_v2.rds` is stored in `R` working directory.
* The structure of those folders recalls the structure of the [shared folder](https://drive.google.com/drive/folders/11LrUmWVdpX2sUCE354xLjKVNywFZNtXR?usp=sharing){target="_blank"}.
* All packages are properly installed.

Now simply run manually (via source button) or via console (using your own corect path) the `B_Polya_Tree_Meta_Analysis_simulation_post.R` script. E.g.:

```{r c1_5, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("B_Polya_Tree_Meta_Analysis_simulation_post.R")
```

### Real data results
Make sure that

* `Star_PFSr_02112023.rds`, `Model_02112023.rds` and `Post_real_data.rds`  are stored in `R` working directory.
* The structure of that folder recalls the structure from the  [shared folder](https://drive.google.com/drive/folders/11LrUmWVdpX2sUCE354xLjKVNywFZNtXR?usp=sharing){target="_blank"}. 
* All packages are properly installed.

Now simply run `B_Polya_Tree_PFS.R` script. E.g.:

```{r c1_6, eval=FALSE, message=FALSE, warning=FALSE, echo=TRUE}
source("B_Polya_Tree_PFS.R", echo=TRUE)
```

# Section 4: {#s5}

Here we we briefly describe files in the shared folder:

- Folder `Dati` Contains data collected reported from studies in various formats.
- Folder `Plots` contains all plots.
- Folder `PTGP_simulation` contains the `.rds` files that store the  results for seach simulation seeds.
- Folder `Scripts` the main scripts used for the analysis.
- Zip Folder `PTGP_simulation.zip`  contains the same files as the unzipped folder, but allows easier downloading.
- File `Model_02112023.rds`, `Post_real_data.rds`, `PTGP_simulation.zip`,  `RESULTS_GPPT_SIM_v2.rds` and  `Star_PFSr_02112023.rds` are the intermediate results files needed to reproduce the plots. 
- `README.RMD` produces this `.html`.
