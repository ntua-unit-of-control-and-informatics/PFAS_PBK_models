Papp_RYU = 1.46e-6*3600 # cm/h, at pH = 7.4 from Ryu et al. (2024) [https://doi.org/10.1016/j.chemosphere.2024.142390]
D_25 <- (5.46*1e-6)*3600 # at 25o,  [cm2/s]* *3600 --> [cm2/h], from Gauthier & Mabury (2024) [https://doi.org/10.1021/acsestwater.4c00631]
eta_25 <- 	0.89 #dynamic viscosity at 25oC
eta_37 <- 	0.6913 #dynamic viscosity at 37oC
D_37 <- D_25*	eta_25/eta_37 #Corrected diffusion coefficient, via the Enstein-Stokes equation
h_monolayer = 16*e-4 #[um]*1e-4-->[cm] Reference [https://doi.org/10.1002/bit.28362] 
h_uwl <- #thickness of unstirred water monolayer