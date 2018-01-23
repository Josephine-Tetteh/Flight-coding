setwd("~/Documents/Project/R-Coding/Flight function")
# R-FUNCTION FOR CALCULATING BIRD FLIGHT PERFORMANCE PARAMETERS##
# TWO FUNCTIONS HAVE BEEN GIVEN IN THIS SCRIPT #
# THE FIRST FUNCTION "BIRDMORPHPARAM" TAKES AS ARGUMENTS, 
# THE MORPHOLOGICAL PARAMETERS OF A USER-DEFINED BIRD AND 
# ITS OUTPUT IS A DATAFRAME OF THESE PARAMETERS#
# THE SECOND FUNCTION "VELO" CALCULATES THE PERFORMANCE PARAMETERS#
# BY TAKING THE MORPHOLOGICAL PARAMETERS DEFINED IN BIRDMORPHPARAM#
# AS INPUTS. ITS OUTPUTS ARE THE PERFORMANCE PARAMETERS IN A DATAFRAME#


# THIS IS THE FIRST FUNCTION
BirdMorphParam = function(BMass, WSpan, WArea, C_d=0.2, ADensity = 1.23,
                          grav=9.8,k=1.2){
  WLoading = BMass / WArea # calculate wing loading
  AR = WSpan^2/WArea  # calculate aspect ratio of wing
  BWeight = BMass*grav # calculate weight of bird
  
  ## flight performance parameters
  ## the following formulations are adapted from C.J. Pennycuick's FLIGHT 1.25 program
  Sb=0.00813*BMass^(0.666) # calculate body frontal area
  flap_freq = BMass^(3/8)*grav^(1/2)*WSpan^(-23/24)*WArea^(-1/3)*ADensity^(-3/8)  #### calculate flapping frequency
  Pbmr = 10^(log10(3.79*BMass^(0.723))) # calculate basal metabolic rate 
  Pmet = 0.23*Pbmr # calculate metabolic power
  Vmp = (0.807*k^0.25*BMass^0.5*grav^0.5)/(1.23^0.5*WSpan^0.5*Sb^0.25*C_d^0.25)  # calculate speed at minimum power
  Pam = (1.05*k^0.75*BMass^(3/2)*grav^(3/2)*Sb^(1/4)*C_d^0.25)/(ADensity^0.5*WSpan^(3/2))   # calculate absolute minimum power
  Pmech = (2*k*(BMass*grav)^2)/(Vmp*pi*WSpan^2*ADensity) + (ADensity*Vmp^3*Sb*C_d)/(2) + (8.4/AR)*Pam  # calculate mechanical power
  Mmusc = 0.17*BMass # calculate muscle mass
  Pm = Pmech/Mmusc  # calculate mass-specific power for the whole muscle
  Pchem = 1.1*(Pmech+Pmet)/0.23  # calculate chemical power
  qlim = (0.30*400000*0.26)/1060 # calculate upper limit to specific power
  Fmito = k*Pm*1060*(10^-6) # calculate  fraction of the muscle mass (or volume) that consists of mitochondria
  Pmax = qlim*Mmusc*(1-Fmito)*flap_freq # calculate maximum power
  Paramdf = data.frame(BMass = BMass, BWeight, WSpan, WArea, AR, WLoading, flap_freq,
                       Vmp, Pmech, Pchem)
  return(Paramdf)
}
# END OF FUNCTION 1 

# NOW THE SECOND FUNCTION

VELO = function(BirdParam, t, x, y, z,
               C_l, C_t, grav = 9.8, ADensity = 1.23, k=1.2) {
  
  BMass = BirdParam$BMass
  WSpan = BirdParam$WSpan
  WArea = BirdParam$WArea
  BWeight = BirdParam$BWeight
  AR = BirdParam$AR
  C_d = BirdParam$C_d
  # t is the time for each step
  # x is the position on the x-axis
  # y is the position on the y-axis
  # z is the position on the z-axis (height)
  # NB: x,y,z should be in metres
  
  # BMass is the body mass in kg
  # WSpan is the wingspan in m
  # WArea is the wing area in m2
  
  # C_d is the drag co-efficient
  # C_l is the lift co-efficient
  # C_t is the thrust co-efficient
  # grav is the acceleration due to gravity; default value is 9.8m/s2
  # ADensity is the Air density; default value is 1.23 kgm/s2
  # k is the induced power factor
  #############################################
  
  ####################################
  if(length(t) != length(x)) stop("t must be the same length as x. t is length: ", length(t), ". x is length: ", length(x))
  object.pos=data.frame(t,x,y,z) # put position values and time into a dataframe
  ####################################
  
 
  #evaluate changes in x,y,z values in order to calculate steplength and climbing angle
  change_x=c(object.pos[1,2],diff(object.pos$x))
  change_y=c(object.pos[1,3],diff(object.pos$y))
  change_z=c(object.pos[1,4],diff(object.pos$z))
  step_length=sqrt(change_x^2+change_y^2+change_z^2) # calculate steplength

  bet = atan(change_z/sqrt(change_x^2+change_y^2))  #angle of climb
  obj.pos = data.frame(object.pos,change_x,change_y, change_z, step_length,bet)
  obj.pos <- transform(obj.pos, beta= ifelse(bet=='NaN', 0, bet))
  
  attach(obj.pos)
  # user components
  # evaluate velocity, lift, drag, thrust, power and energy of user from dataframe (obj.pos)
  vel_user = step_length/t # user-defined velocity
  lift_user = 0.5*ADensity*WArea*C_l*vel_user^2 # user-defined lift
  drag_user = 0.5*ADensity*WArea*C_d*vel_user^2 # user-defined drag
  thrust_user =  0.5*ADensity*WArea*C_t*vel_user^2 + BWeight*sin(beta) # user-defined thrust
  power_user = vel_user*thrust_user # user-defined power
  energy_user = thrust_user*step_length # user-defined energy
  
  # climbing components
  # evaluate climbing velocity, thrust and power
  vel_climb = c()
  thrust_climb =c()
  power_climb = c()
  vel_climb = sqrt((2*BWeight*cos(beta))/(C_l*ADensity*WArea)) # climbing velocity
  thrust_climb = 0.5*ADensity*WArea*C_t*vel_climb^2 + BWeight*sin(beta) # thrust in climb
  power_climb = vel_climb*thrust_climb # power for climb
  
 
  # OUTPUT
  morph_parameters = data.frame(BMass, WSpan, WArea, AR)
  dframe = data.frame(vel_user,thrust_user,power_user,energy_user, vel_climb, thrust_climb, power_climb)
  MyList<- list("a"=morph_parameters, "b"=dframe) 
  newlist = data.frame(list(MyList))
  return(list(MyList)) #, morph_parameters)
}

# END OF FUNCTION 2
##########################################################

# AN EXAMPLE
 df = read.csv("DAY1-converted.csv")  # sample data

 colnames(df) = c("y", "x", "z","change.t") # change column names of data to suit function
 
 # implement function 1
 BirdPa = BirdMorphParam(BMass=2,WSpan=3,WArea=4,C_d=0.2)

 # implement function 2
 q=VELO(BirdPa, t=df$change.t, x=df$x, y=df$y, z=df$z, C_l=0.5, C_t=0.1)
 
 # END OF EXAMPLE.
 

 