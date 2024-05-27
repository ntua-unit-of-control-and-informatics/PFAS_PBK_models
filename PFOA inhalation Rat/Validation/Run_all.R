path <- "C:/Users/ptsir/Documents/GitHub/PFAS_PBK_models/PFOA inhalation Rat/Validation"

setwd(path)
source("Validate_Cui2008.R")
setwd(path)
source("XY_Plot_tissues.R")

setwd(path)
source("Validate_Kemper2003_serum.R")
setwd(path)
source("XY_Plot_serum.R")

setwd(path)
source("Validate_Kemper2003_excreta.R")
setwd(path)
source("XY_Plot_excreta.R")

