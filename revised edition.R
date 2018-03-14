

# R-FUNCTION FOR CALCULATING BIRD FLIGHT PERFORMANCE PARAMETERS##
# SIX FUNCTIONS HAVE BEEN GIVEN IN THIS SCRIPT #

# THE FIRST FUNCTION "BIRDMORPHPARAM" TAKES AS ARGUMENTS, THE MORPHOLOGICAL PARAMETERS OF A USER-DEFINED BIRD AND 
# ITS OUTPUT IS A DATAFRAME OF THESE PARAMETERS#

# THE SECOND FUNCTION "FLIGHTSPEEDCOMPONENTS" CALCULATES THE X,Y,Z-COMPONENTS OF THE FLIGHT SPEED FROM THE DATA GIVEN BY THE USER

# THE THIRD FUNCTION "TAS" EVALUATES THE TRUE AIRSPEED GIVEN THAT THE USER INPUTS THE X,Y,Z-COMPONENTS OF THE WIND SPEED

# THE FOURTH FUNCTION "TAS2" EVALUATES TRUE AIRSPEED GIVEN THAT THE USER INPUTS THE VALUE OF WINDSPEED AND 2 DIRECTIONAL ANGLES

# THE FIFTH FUNCTION EVALUATES THE FLIGHT PERFORMANCE PARAMETERS SUCH AS POWER FROM THE BIRD, MECHANICAL AND CHEMICAL POWER
###############################################################################################################################################

# FUNCTION 1
# Morphological data from bird
BirdMorphParam = function(BMass, WSpan, WArea, C_db=0.2, C_dpro=0.2, ADensity = 1.23,
                          grav=9.8,k=1.2){
  WLoading = BMass / WArea # calculate wing loading
  AR = WSpan^2/WArea  # calculate aspect ratio of wing
  BWeight = BMass*grav # calculate weight of bird
  
  ## flight performance parameters
  ## the following formulations are adapted from C.J. Pennycuick's FLIGHT 1.25 program
  Sb=0.00813*BMass^(0.666) # calculate body frontal area
  flap_freq = BMass^(3/8)*grav^(1/2)*WSpan^(-23/24)*WArea^(-1/3)*ADensity^(-3/8)  #### calculate flapping frequency
  Pbmr = 10^(log10(3.79*BMass^(0.723))) # calculate basal metabolic rate 
  Pmet = 0.23*Pbmr # calculate metabolic Power
  Vmp = (0.807*k^0.25*BMass^0.5*grav^0.5)/(1.23^0.5*WSpan^0.5*Sb^0.25*C_db^0.25)  # calculate speed at minimum Power
  Pam = (1.05*k^0.75*BMass^(3/2)*grav^(3/2)*Sb^(1/4)*C_db^0.25)/(ADensity^0.5*WSpan^(3/2))   # calculate absolute minimum Power
  Mmusc = 0.17*BMass # calculate muscle mass
 
  #print results
  BirdParam = data.frame(BMass = BMass, BWeight, WSpan, WArea, AR, WLoading, C_db, C_dpro, flap_freq, Sb, Pmet, Pbmr, Vmp, Pam)
  return(BirdParam)
}
# END OF FUNCTION 1 

# FUNCTION 2
## To calculate flight speed components from the position of the bird (Vuser)
FlightSpeedComponents = function(t,x,y,z){
  dataday11=data.frame(t, x=x, y=y, z=z)
  change_x=c(dataday11[1,2],diff(dataday11$x))
  change_y=c(dataday11[1,1],diff(dataday11$y))
  change_z=c(dataday11[1,3],diff(dataday11$z))
  step_length = sqrt(change_x^2+change_y^2+change_z^2)
  beta = atan(change_z/sqrt(change_x^2+change_y^2))  #angle of climb
  beta= ifelse(is.nan(beta), 0, beta)
  dataday11=data.frame(dataday11,change_x,change_y, change_z, step_length,beta)
  attach(dataday11)
  FlightSpeed.x = change_x/t # speed of flight in the x-direction
  FlightSpeed.y = change_y/t # speed of flight in the y-direction
  FlightSpeed.z = change_z/t # speed of flight in the z-direction
  FlightSpeed = sqrt(FlightSpeed.x^2 + FlightSpeed.y^2 + FlightSpeed.z^2) #speed of flight due to position vectors of the bird
  Vframe = data.frame(FlightSpeed.x, FlightSpeed.y, FlightSpeed.z, FlightSpeed)
  return(Vframe)
}
# END OF FUNCTION 2

# FUNCTION 3
## To calculate effect of the wind on flight speed of the bird
# Alternative 1 (user inputs x,y,z components of wind speed and the function, TAS, calculates true airspeed) 
TAS = function(FlightSpeedComponents, WindSpeed.x, WindSpeed.y, WindSpeed.z){
   FlightSpeed.x = FlightSpeedComponents$FlightSpeed.x
   FlightSpeed.y = FlightSpeedComponents$FlightSpeed.y
   FlightSpeed.z = FlightSpeedComponents$FlightSpeed.z
   if(missing(WindSpeed.x))   WindSpeed.x = 0;
   if(missing(WindSpeed.y))   WindSpeed.y = 0;
   if(missing(WindSpeed.z))   WindSpeed.z = 0;
  WindSpeed.total = sqrt((WindSpeed.x^2 + WindSpeed.y^2 + WindSpeed.z^2))
  TrueAirSpeed = sqrt((FlightSpeed.x + WindSpeed.x)^2 + (FlightSpeed.y + WindSpeed.y)^2 + (FlightSpeed.z + WindSpeed.z)^2)
  vuframe = data.frame(TrueAirSpeed)
  return(vuframe)
}
# END OF FUNCTION 3

# FUNCTION 4
# Alternative 2(user inputs wind speed(in m/s) and wind directions(in degrees) and then function calculates true airspeed)
TAS2 = function(FlightSpeedComponents, V_WIND, theta, phi){  #V_WIND is wind speed, theta and phi are wind directions
  FlightSpeed.x = FlightSpeedComponents$FlightSpeed.x
  FlightSpeed.y = FlightSpeedComponents$FlightSpeed.y
  FlightSpeed.z = FlightSpeedComponents$FlightSpeed.z
  if(missing(phi)) WindSpeed.z=0;
  WindSpeed.x = V_WIND*cos(theta)
  WindSpeed.y = V_WIND*sin(theta)
  WindSpeed.z = V_WIND*sin(phi)
  WindSpeed.total = sqrt((WindSpeed.x^2 + WindSpeed.y^2 + WindSpeed.z^2))
  TrueAirSpeed = sqrt((FlightSpeed.x + WindSpeed.x)^2 + (FlightSpeed.y + WindSpeed.y)^2 + (FlightSpeed.x + WindSpeed.z)^2)
  vuframe = data.frame(TrueAirSpeed)
  return(vuframe)
}
# END OF FUNCTION 4


# FUNCTION 5
# Calculate flight parameters (mechanical Power, chemical Power, flight type)
FlightPerformance = function(BirdParam,FlightSpeedComponents, TAS, TAS2, t, x, y, z,
                C_l, C_t, grav = 9.8, ADensity = 1.23, k=1.2) {
  if(missing(TAS)) TAS2;
  if(missing(TAS2)) TAS;
  
  BMass = BirdParam$BMass
  WSpan = BirdParam$WSpan
  WArea = BirdParam$WArea
  BWeight = BirdParam$BWeight
  AR = BirdParam$AR
  C_db = BirdParam$C_db
  C_dpro = BirdParam$C_dpro
  Sb = BirdParam$Sb
  Pam = BirdParam$Pam
  Pmet = BirdParam$Pmet
  TrueAirSpeed = TAS$TrueAirSpeed
  #beta = FlightSpeedComponents$beta
  ## Prepare dataframe
   dataday11=data.frame(t, x=x, y=y, z=z)
   change_x=c(dataday11[1,2],diff(dataday11$x))
   change_y=c(dataday11[1,1],diff(dataday11$y))
  change_z=c(dataday11[1,3],diff(dataday11$z))
 # dataday11 = FlightSpeedComponents$dataday11
  dataday11=data.frame(dataday11,change_x,change_y, change_z, TrueAirSpeed,beta)
 #
  ## drag components and total drag
  d_ind = (2*k*(BMass*grav)^2)/(TrueAirSpeed^2*pi*WSpan^2*ADensity)
  d_pro = 0.5*ADensity*WArea*C_dpro*TrueAirSpeed^2 
  d_par = 0.5*ADensity*Sb*C_db*TrueAirSpeed^2 
  
  drag_total = d_ind + d_par + d_pro #total aerodynamic drag
  #################################################################
  transform(dataday11)
  ### flight types
  for(i in 1:length(dataday11$x)){
    if(dataday11$change_z[i]<0){
      dataday11$flighttype[i] = "1"  ##### descent
    } else if(dataday11$change_z[i]>0){
      dataday11$flighttype[i] = "2"   #### climb
    }  else if(dataday11$change_z[i]==0){
      dataday11$flighttype[i] = "3"  #### steady
    }}
  
  ### Descending 
  drag_desc = drag_total+BWeight*sin(beta)
  
  ## Straight
  StrSpeed = sqrt((2*BWeight)/(C_l*ADensity*WArea))  # velocity in straight flight
  d_ind_str = (2*k*(BMass*grav)^2)/(StrSpeed^2*pi*WSpan^2*ADensity)
  d_pro_str = 0.5*ADensity*WArea*C_dpro*StrSpeed^2 
  d_par_str = 0.5*ADensity*Sb*C_db*StrSpeed^2 
  StrDrag = d_ind_str + d_par_str + d_pro_str # total drag in a straight flight
  
  # Climbing 
  thrust_climb = drag_total + BWeight*sin(beta)  ## thrust in a climb
  #Calculate climbing velocity
  attach(dataday11)
  ClimbSpeed = c()
  for(i in 1:length(dataday11$x)){
    if((dataday11$flighttype[i]==2)){
      ClimbSpeed[i] = sqrt((2*BWeight*cos(beta[i]))/(C_l*ADensity*WArea))}  ### climbing velocity
    else{
      ClimbSpeed[i]= NA
    }
  }

  ###Checks for climbing and descending
  #Climb is feasible if TrueAirSpeed is greater or equals climbing velocity otherwise it is not feasible.
  Power=c()
  count = 0
  for(i in 1:length(dataday11$x)){
    if((dataday11$flighttype[i]==2)&(TrueAirSpeed[i]<as.numeric(ClimbSpeed[i]))){
      dataday11$feasib[i] = "not feas climb"
      Power[i] = NA
      count = count+1
    }
    else{
      dataday11$feasib[i] = "feas climb"
      Power[i] = TrueAirSpeed[i]*thrust_climb[i]
    }
    if(dataday11$flighttype[i]==1){
      #drag_desc = drag_total+BWeight*sin(beta)
      Power[i] = TrueAirSpeed[i]*drag_desc[i]
      dataday11$feasib[i] = "desc"
    }
    if((dataday11$flighttype[i]==3)&(TrueAirSpeed[i]>StrSpeed)){
      dataday11$feasib[i] = "feas str"
      Power[i] = TrueAirSpeed[i]*StrDrag
    }
    else if((dataday11$flighttype[i]==3)&(TrueAirSpeed[i]<StrSpeed)){
      dataday11$feasib[i] = "notfeas str"
      Power[i] = NA
      count = count+1
    }
  }

  

  power2 = c()
  for(i in 1:length(dataday11$x)){
    if(dataday11$flighttype[i]==2){
      power2[i] = TrueAirSpeed[i]*thrust_climb[i]
    }
    if(dataday11$flighttype[i]==1){
      #drag_desc = drag_total+weight*sin(beta)
      power2[i] = abs(TrueAirSpeed[i]*drag_desc[i])
      
    }
    if(dataday11$flighttype[i]==3){
      
      power2[i] = TrueAirSpeed[i]*StrDrag
    }
  }
  #dataday11=data.frame(dataday11,power,power2)
  
  ###### end of checks
  
  # ## evaluate mechanical and chemical Power
 Pmech_data = (2*k*(BMass*grav)^2)/(TrueAirSpeed*pi*WSpan^2*ADensity) + (ADensity*TrueAirSpeed^3*Sb*C_db)/(2)# + (8.4/Ra)*Pam  ##mechanical Power
 Pchem_data = 1.1*(Pmech_data+Pmet)/0.23  ##chemical Power
  # 
  # add to dataframe
  mydata = data.frame(dataday11,ClimbSpeed, Power,power2, Pmech_data,Pchem_data,count)

  # filter Inf values from dataset in order to plot
  Newdata = mydata[!is.infinite(Pmech_data), ]

  #calculate minimum Power velocity from data
  Vmp_data= TrueAirSpeed[which(Newdata$Pmech_data==min(Newdata$Pmech_data))]
 # detach(dataday11)#, character.only = TRUE)
  # OUTPUTdetach(pkg, character.only = TRUE)
#SW = suppressWarnings(Newdata)
  return(Newdata) 
  
}
# END OF FUNCTION 5

# FUNCTION 6
# Plot mechanical Powercurve
PowerCurve = function(FlightPerformance, TAS, TAS2){
  if(missing(TAS))  TrueAirSpeed = TAS2$TrueAirSpeed;
  if(missing(TAS2)) TrueAirSpeed = TAS$TrueAirSpeed;
   AA = TrueAirSpeed
   BB = FlightPerformance$Pmech_data
   Vmp = FlightPerformance$Vmp_data
  Powercurve = plot(x= AA, y= BB,type = "p", log="y", xlab="Velocity (m/s)", ylab="Mechanical Power (W)", main="Mechanical Powercurve")
              abline(v = FlightPerformance$Vmp_data, col = "blue", lty = 2)
  return(Powercurve)

}
# END OF FUNCTION 6

