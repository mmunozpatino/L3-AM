
#################################################################################################
# Plot
#################################################################################################
plot.double.pendulum.init <- function(L=L, H=1.5, W=2, e=0.05, pen.width=0.01, name="Pendulo"){

	#
  x11() #plot in a new window

	plot(x=c(), y=c(),xlim=c(-1,1), ylim=c(0,2.5),main=name,xlab ="X", ylab ="Y")

	#title(main="Péndulo", sub="Red:ANN,Blue:function",xlab="X", ylab="Y")

	#(-1,1.5), (1,1.5), (1,1.45), (-1,1.45)

  init.polygon( H=1.5, W=2, e=0.05, pen.width=0.01)

}
#-----------------------------------------
rotate <- function(center=c(0,0), x, y,T) {
	x <- x - center[1]
	y <- y - center[2]

	xp1 <- x*cos(T) - y*sin(T)
	yp1 <- x*sin(T) + y*cos(T)

	list(x=xp1 + center[1], y=yp1 + center[2])
}
#-----------------------------------------
init.polygon <- function(H=1.5, W=2, e=0.05, pen.width=0.01)
{
  H <<- H
  W <<- W
  e <<- e   # espesor del techo
  d <<- pen.width
  polygon(x=c(-H/2,H/2,H/2,-H/2), y=c(H,H,H-e,H-e), col="black")
}
#-----------------------------------------
plot.double.pendulum <- function(T1, T2, L, col="blue"){

	#PRIMER PENDULO
		#(-d/2,H-e), (d/2,H-e), (d/2,H-e-L), (-d/2, H-e-L)


		x <- c(-d/2,d/2,d/2,-d/2)
		y <- c(H-e,H-e, H-e-L, H-e-L)

		p1 <- rotate(center=c(0, H-e), x=x, y=y, T1)

		polygon(x=p1$x, y=p1$y, col=col, border=NA)

	#SEGUNDO PENDULO (simplemente L mas abajo)
		#(-d/2,H-e), (d/2,H-e), (d/2,H-e-L), (-d/2, H-e-L)
		# Es decir, igual que el pendulo anterior pero L mas abajo
		x2 <- x
		y2 <- y-L

		#Rotamos el punto de rotacion (0,H-e-L)
		x0 <- 0
		y0 <- H-e-L
		p0 <- rotate(center=c(0, H-e), x=x0, y=y0, T1)

		x2 <- x2 + (p0$x-x0)
		y2 <- y2 + (p0$y-y0)
		p2 <- rotate(center=c(p0$x, p0$y), x=x2, y=y2, T2)

		polygon(x=p2$x, y=p2$y, col=col, border=NA)

}

#t[index.t] <- t[index.t] - t.pi[index.t] * (2 * pi)
plot.thetas <- function(it, t1, t2, tn1, tn2, e1, e2, dynamic=FALSE) {
	ymax <- pi
	i.len <- length(it)
	col <- c(1,3,2)
	#xmax <- 0

	# Imprime los 6 vectores en dos ventanas diferentes.
	m1 <- build.matrix(t1[it], tn1[it], e1[it])
	m2 <- build.matrix(t2[it], tn2[it], e2[it])

  X11()
	matplot(it, m1, xlim=c(min(it),max(it)), ylim=c(-ymax,ymax), type="l", main="Theta1", xlab="#iterations", ylab="Theta - Error", col=col )
	legend("topright", legend=c("Theta1", "ThetaN1", "Error"), col=col, cex=.8, lty=c(1,2,3))

	X11()
	matplot(it, m2, xlim=c(min(it),max(it)), ylim=c(-ymax,ymax), type="l", main="Theta2", xlab="#iterations", ylab="Theta - Error", col=col )
	legend("topright", legend=c("Theta2", "ThetaN2", "Error"), col=col, cex=.8, lty=c(1,2,3))

}



#################################################################################################
#
#################################################################################################
build.matrix <- function(t, tn, e) {
		t <- normalize.angles(t)
		tn <- normalize.angles(tn)

		m <- as.matrix(cbind(t,tn,e))

		m
}

#################################################################################################
#
#################################################################################################
normalize.angles <- function(t) {
	t.pi <- round(t/(2*pi))
	index.t <- t >= pi
	t[index.t] <- t[index.t] - t.pi[index.t]*(2*pi)
	index.t <- t < -pi
	t[index.t] <- t[index.t] + t.pi[index.t]*(2*pi)

	t
}
#################################################################################################
#################################################################################################
# Method to split the set data into two subset where one of subsets has perc
# percent of elements of data
#
split.data <- function(data, perc=1) {

  index <- 1:(nrow(data) * perc)

  train <- data[index,]
  test <- data[-index,]

  return(list(train=train, test=test))
}
