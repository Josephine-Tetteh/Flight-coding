# FUNCTION 8
#' Plot the input trajectory according to the types of flight
#' @param x x-position from data
#' @param y y-position from data
#' @param z z-position from data
#' @param FlightPerformance requires flighttype information from \code{FlightPerformance}
#' @return 3D scatterplot of flight types from the input trajectory
#' @import scatterplot3d
#' @examples
#' utils::data(trajectory)
#' trajectory = trajectory
#' bird = BirdMorphParam(BMass = 0.7710, WSpan = 0.98, WArea = 0.119, C_db = 0.1)
#' FSComponents=FlightSpeedComponents(t=trajectory[,4],
#'                    x=trajectory[,1], y=trajectory[,2], z=trajectory[,3])
#' Tairspeed = TrueAirSpeed1(FSComponents)
#' FPerformance = FlightPerformance(bird,FSComponents, Tairspeed,
#'                    C_l = 0.5, C_t = 0.1)
#'
#' PlotFlightType(x=trajectory[,1], y=trajectory[,2], z=trajectory[,3], FPerformance)

#' @export

PlotFlightType = function(x,y,z, FlightPerformance){
  pol = c("blue", "red", "green")
  pol <- pol[as.numeric(FlightPerformance$flighttype)]
  scatterplot3d(x=x[1:length(FlightPerformance$x)], y=y[1:length(FlightPerformance$y)],
                z=z[1:length(FlightPerformance$z)], color=pol,type = 'o',cex.symbols = 0.05,
                               xlab = 'Latitude (m)', ylab = 'Longitude (m)', zlab = 'Height (m)',cex.lab = 1.1, cex.axis = 0.9,
                               main = 'Different flight types',cex.main=1.2)
  legend("right", legend = c("Descent","Climb","Steady"),
         col =  c("blue", "red","green"), pch = 16, inset = -0.06, xpd = TRUE, horiz = FALSE,cex=0.8 )
}


