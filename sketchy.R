# This file defines the necessary functions to implement the "sketchy lines
# algorithm" by Wood et al 2012 (http://dx.doi.org/10.1109/TVCG.2012.262). No
# external libraries are used. The final product is a function `sketchylines(x,
# y, roughness, strokes, ...)` that can be used as a replacement for the base
# graphics function `lines(x, y, ...)`

#' Helper function used by catmull_rom_spline
tp1 = function(ti, Pi, Pj, alpha) {
  sum((Pi - Pj)^2) ^ (alpha/2) + ti
}

#' Calculate Catmull-Rom spline between points P1 and P2 
#' see https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline
catmull_rom_spline = function(P0, P1, P2, P3, n_points=50, alpha=.5) {

  t0 = 0
  t1 = tp1(t0, P0, P1, alpha)
  t2 = tp1(t1, P1, P2, alpha)
  t3 = tp1(t2, P2, P3, alpha)

  tt = seq(t1, t2, len=n_points)

  A1 = (t1-tt)/(t1-t0)*P0 + (tt-t0)/(t1-t0)*P1
  A2 = (t2-tt)/(t2-t1)*P1 + (tt-t1)/(t2-t1)*P2
  A3 = (t3-tt)/(t3-t2)*P2 + (tt-t2)/(t3-t2)*P3
  B1 = (t2-tt)/(t2-t0)*A1 + (tt-t0)/(t2-t0)*A2
  B2 = (t3-tt)/(t3-t1)*A2 + (tt-t1)/(t3-t1)*A3
  CC = (t2-tt)/(t2-t1)*B1 + (tt-t1)/(t2-t1)*B2

  return(CC)

}


#' Interpolate n points by a catmull rom splines (P is a matrix with 2 columns,
#' x and y coords)
catmull_rom_chain = function(P, n_points=50, alpha=.5) {
    nP = nrow(P)
    C = c()
    for (i in 1:(nP-3)) {
      C = c(C, catmull_rom_spline(P[i, ], P[i+1, ], P[i+2, ], 
                                  P[i+3, ], n_points, alpha))
    }
    C = t(matrix(C, nrow=2))
    return(C)
}

#' Connect points xy1, xy2 with a "sketchy" line as described in Wood et al
#' 2012 (minor modification in calculation of roughness)
sketchyline = function(xy1, xy2, roughness=0.01, strokes=2, ...) {

  ll = sqrt(sum((xy1 - xy2)^2))
  
  # calculate new start and end points
  xy0 = xy1 - (xy2 - xy1)
  xy3 = xy2 + (xy2 - xy1)
  
  # calculate midpoints
  xymid1 = (xy2 + xy1) / 2
  xymid2 = (3 * xy2 + xy1) / 4
  
  # combine points
  xynew = rbind(xy0, xy1, xymid1, xymid2, xy2, xy3)
  
  # create copies of xy with added noise
  noise = replicate(strokes, c(rep(0,2), runif(8, -roughness*ll, roughness*ll), rep(0,2)), simplify=FALSE)
  xynewn = lapply(1:strokes, function(ii) xynew + noise[[ii]])
  
  # calculate splines
  xynewns = lapply(1:strokes, function(ii) catmull_rom_chain(xynewn[[ii]]))
  
  for (ii in 1:strokes) {
    lines(xynewns[[ii]][, 1], xynewns[[ii]][, 2], ...)
  }

}

#' Connect n points by "sketchy" lines, `sketchylines` can be used as a
#' replacement for `lines`
sketchylines = function(x, y, roughness=0.01, strokes=2, ...) {
  n = length(x) 
  for (ii in 1:(n-1)) {
    sketchyline(c(x[ii], y[ii]), c(x[ii+1], y[ii+1]), roughness, strokes, ...)
  }
}

# infinite loop that plots randomly rotated squares with random roughness,
# strokes and colours
while(1) {
  plot(NULL, xlim=c(-1, 1), ylim=c(-1, 1))
  phi = pi * runif(1, -1, 1)
  offset = runif(1, -.3, .3)
  xy = t(matrix(c(0,0,0,1,1,1,1,0,0,0)-.5, nrow=2)) %*% 
       matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), nrow=2) + 
       offset
  xx = xy[,1]
  yy = xy[,2]
  rr = 0.03 * runif(1)
  st = sample(1:10, 1)
  cl = sample(rainbow(10), 1)
  sketchylines(xx, yy, roughness=rr, strokes=st, col=cl)
  Sys.sleep(1)
}


