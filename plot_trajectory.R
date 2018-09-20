# FUNCTION 7
#' Plot input trajectory
#' @param x x-position from data
#' @param y y-position from data
#' @param z z-position from data
#' @return 3D scatterplot of input trajectory
#' @import scatterplot3d
#' @examples
#' utils::data(trajectory)
#' trajectory = trajectory
#' PlotTrajectory(x=trajectory[,1], y=trajectory[,2], z=trajectory[,3])
#' @export

PlotTrajectory  = function(x,y,z){
  trajplot = scatterplot3d(x=x, y=y, z=z, pch=20, type = "o",cex.symbols = 0.05,
                                          cex.lab = 1.1, cex.axis = 0.9, main = 'Flight trajectory',cex.main=1.0)
}
# END OF FUNCTION 7
