uniformPoints.cube = function(N, l) {
	# Uniformly distributes N points in a cube of side l
	# The origin is located in the center of the cube
	# 
	# Args:
	# 	N: number of points
	# 	l: cube side
	# 
	# Returns:
	# 	A 3-columns data.frame with the coordinates of the points
	l <- l/2
	points <- as.data.frame(matrix(runif(N*3, -l, l), ncol=3))
	colnames(points) <- c('x', 'y', 'z')
	return(points)
}
uniformPoints.sphere = function(N, r) {
	# Uniformly distribute N points in a sphere with radius r
	# The origin is located at (0,0,0)
	# 
	# Args:
	# 	N: number of points
	# 	r: sphere radius
	# 	
	# Returns:
	# 	A 3-columns data.frame with the coordinates of the points
	points <- data.frame()
	while(N > nrow(points)) {
		candidates <- uniformPoints.cube(N, r)
		ocd <- unlist(lapply(1:nrow(candidates), FUN=function(i, c) { sqrt(c$x[i]**2 + c$y[i]**2 + c$z[i]**2) }, candidates ))
		inside <- candidates[ocd <= r,]
		if ( 0 != nrow(inside) ) points <- rbind(points, inside[1:(min(nrow(inside),N-nrow(points))),])
	}
	return(points)
}
uniformPoints.cuboid = function(N, a, b, c) {
	# Uniformly distributes N points in a cuboid of dimensions a, b and c
	# The origin is located in the center of the cuboid
	# 
	# Args:
	# 	N: number of points
	# 	a, b, c: half-sizes of the cuboid
	# 
	# Returns:
	# 	A 3-columns data.frame with the coordinates of the points
	points <- as.data.frame(cbind(runif(N, -a, a), runif(N, -b, b), runif(N, -c, c)))
	colnames(points) <- c('x', 'y', 'z')
	return(points)
}
uniformPoints.ellipsoid = function(N, a, b, c) {
	# Uniformly distributes N points in an ellipsoid of dimensions a, b and c
	# The origin is located in the center of the ellipsoid
	# 
	# Args:
	# 	N: number of points
	# 	a, b, c: half-sizes of the ellipsoid
	# 
	# Returns:
	# 	A 3-columns data.frame with the coordinates of the points
	points <- data.frame()
	while(N > nrow(points)) {
		candidates <- uniformPoints.cuboid(N, 2*a, 2*b, 2*c)
		ocd <- unlist(lapply(1:nrow(candidates), FUN=function(i, a, b, c, l) { (l$x[i]/a)**2+(l$y[i]/b)**2+(l$z[i]/c)**2 }, a, b, c, candidates ))
		inside <- candidates[ocd <= 1,]
		if ( 0 != nrow(inside) ) points <- rbind(points, inside[1:(min(nrow(inside),N-nrow(points))),])
	}
	return(points)
}