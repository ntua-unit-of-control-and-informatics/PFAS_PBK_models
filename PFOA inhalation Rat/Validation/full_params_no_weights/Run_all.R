path <- "C:/Users/dpjio/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation/full_params_no_weights"
setwd(path)
load("full_params_no_weights.RData")
source("Validate_Cui2008.R")
rm(list = setdiff(ls(), c("path", "AAFE_Cui")))

setwd(path)
source("XY_Plot_tissues.R")


setwd(path)
load("full_params_no_weights.RData")
source("Validate_Kemper2003_serum.R")
rm(list = setdiff(ls(),  c("path", "AAFE_Cui", "AAFE_Kemper_serum")))
setwd(path)
source("XY_Plot_serum.R")

setwd(path)
load("full_params_no_weights.RData")
source("Validate_Kemper2003_excreta_Worley.R")
source("Validate_Kemper2003_excreta_Loccisano.R")
rm(list = setdiff(ls(),  c("path", "AAFE_Cui", "AAFE_Kemper_serum",
                                            "AAFE_Worley", "AAFE_Loccisano")))
setwd(path)
source("XY_Plot_excreta.R")

df_AAFE <- data.frame("AAFE" = c(AAFE_Cui, AAFE_Kemper_serum, AAFE_Worley, AAFE_Loccisano ))
rownames(df_AAFE) <- c("tissues", "serum", "excreta_worley", "excreta_loccisano")
write.csv(df_AAFE,"AAFE_validation.csv")


