# Supplementary materials for survival variable importance paper

This repository contains code to reproduce the analyses in ["Assessing variable importance in survival analysis using machine learning"](https://arxiv.org/abs/2311.12726) by Wolock, Gilbert, Simon, and Carone (In press, *Biometrika*, 2024+). 

The methods described in the above paper are implemented in the `survML` package (version 1.2.0 or higher), available on CRAN. 

The code depends on the following `R` packages: 

* `cowplot`: Available on CRAN.
* `dichromat`: Available on CRAN.
* `extrafont`: Available on CRAN.
* `ggpubr`: Available on CRAN.
* `grid`: Available on CRAN.
* `gridExtra`: Available on CRAN.`
* `gtools`: Available on CRAN.
* `mboost`: Available on CRAN.
* `randomForestSRC`: Available on CRAN. 
* `SuperLearner`: Available on CRAN.
* `survex`: Available on CRAN. 
* `survival`: Available on CRAN.
* `survML`: Available on CRAN. 
* `survSuperLearner`: Available on Github at https://github.com/tedwestling/survSuperLearner.
* `tidyverse`: Available on CRAN.

Note that `survSuperLearner`, `randomForestSRC`, and `survex` are only used as comparator methods in the simulation studies. If you have difficulty installing them, you can simply leave those methods out. Similarly, `cowplot`, `dichromat`, `extrafont`, `ggpubr`, `grid`, and `gridExtra` are needed only for making figures. 
