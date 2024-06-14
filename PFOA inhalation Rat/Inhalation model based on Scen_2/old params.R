if (sex == "M"){
  RAFOatp_k <- 1.000000e+02
  RAFOat1 <- 9.778311e+00
  RAFOat3 <- RAFOat1
  RAFUrat <- RAFOat1
  
  RAFOatp_l <- 6.308075e-01
  RAFOatp2_l <- RAFOatp_l
  RAFNtcp <- RAFOatp_l
  
  KmK_baso <- 9.297172e+03
  VmK_baso <- 7.541167e+03
  
  
}else if(sex == "F"){
  RAFOatp_k <- 1.613936e-03
  RAFOat1 <- 1.000000e+02 
  RAFOat3 <- RAFOat1
  RAFUrat <- RAFOat1
  
  RAFOatp_l <- 1.430072e-01
  RAFOatp2_l <- RAFOatp_l
  RAFNtcp <- RAFOatp_l
  
  KmK_baso <- 1.000000e+03
  VmK_baso <- 1.370572e+03
  
}
RAFOatp2_Int <- 4.104426e+01
#permeabilities correction factor
CF_Peff <- 2.602230e+06 
P_liver_bile <- 2.000878e+00 