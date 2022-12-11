# PBK models for PFAS 
This repository contains R code for developing PBK model for describing PFAS biodistribution in multiple species and deploying as web services on Jaqpot (https://www.jaqpot.org/). This work is part of the SCENARIOS EU H2020 project (https://scenarios-project.eu/). The list of available PBK models will be updated throughout the lifespan of the project.

## List of PBK models
- Humans:
- Rats:
    * PFOA [male, female] (Loccisano et al., 2012)
    * PFOS [male, female] (Loccisano et al., 2012)
    * PFHxS [male, female] (Kim et al., 2018)
    * PFNA [male, female] (Kim et al., 2019)
    * PFDA [male, female] (Kim et al., 2019)
- Mice:
- Monkeys:
- Fish

## Folder Contents
*  Bernstein et al.2021 --> This folder transfers the PFOA, PFOS, PFHxS, PFNA and PFDA PBK code presented in Bernstein et al. (2021) to R code in an appropriate format, so that it can be readily exposed as a web service on Jaqpot.
## Necessary R packages
]* `deSolve`  --> A library containing functions for solving ODE systems
* `jaqpotr` --> A library containing functions for uploading models on the Jaqpot server. Installation instructions can be found in https://www.jaqpot.org/docs/r

## References
* Bernstein, A. S., Kapraun, D. F., & Schlosser, P. M. (2021). A Model Template Approach for Rapid Evaluation and Application of Physiologically Based Pharmacokinetic Models for Use in Human Health Risk Assessments: A Case Study on Per- and Polyfluoroalkyl Substances. Toxicological sciences : an official journal of the Society of Toxicology, 182(2), 215–228. https://doi.org/10.1093/toxsci/kfab063
* Kim, S. J., Shin, H., Lee, Y. B., & Cho, H. Y. (2018). Sex-specific risk assessment of PFHxS using a physiologically based pharmacokinetic model. Archives of toxicology, 92(3), 1113–1131. https://doi.org/10.1007/s00204-017-2116-5
* Kim, S. J., Choi, E. J., Choi, G. W., Lee, Y. B., & Cho, H. Y. (2019). Exploring sex differences in human health risk assessment for PFNA and PFDA using a PBPK model. Archives of toxicology, 93(2), 311–330. https://doi.org/10.1007/s00204-018-2365-y
* Loccisano, A. E., Campbell, J. L., Jr, Butenhoff, J. L., Andersen, M. E., & Clewell, H. J., 3rd (2012). Comparison and evaluation of pharmacokinetics of PFOA and PFOS in the adult rat using a physiologically based pharmacokinetic model. Reproductive toxicology (Elmsford, N.Y.), 33(4), 452–467. https://doi.org/10.1016/j.reprotox.2011.04.006
