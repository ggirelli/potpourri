#!/usr/bin/Rscript

library(parallel)
library(tiff)

# PARAMETERS
#--------------------------------------------

N <- 100	# Number of points/dots

a <- 4	# x ellipsoid half-size
b <- 4	# y ellipsoid half-size
c <- 4	# z ellipsoid half-size

px <- .13	# x pixel size
py <- .13	# y pixel size
sh <- .2	# slice height

dr <- .3	# dot radius

uom <- 'um'	# units of measure (used in prompts)
decs <- 3	# number of decimals for intensity values
cores <- 4	# number of cores for parallelized analysis

fname <- 'test.tif'	# output image filename

# CODE
#--------------------------------------------

uniformPoints.parallelepiped = function(N, a, b, c) {
	# Uniformly distributes N points in a parallelepiped of dimensions a, b and c
	# The origin is located in the center of the parallelepiped
	# 
	# Args:
	# 	N: number of points
	# 	a, b, c: dimensions of the parallelepiped
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
	# 	a, b, c: dimensions of the ellipsoid
	# 
	# Returns:
	# 	A 3-columns data.frame with the coordinates of the points
	points <- data.frame()
	while(N > nrow(points)) {
		candidates <- uniformPoints.parallelepiped(N, 2*a, 2*b, 2*c)
		ocd <- unlist(lapply(1:nrow(candidates), FUN=function(i, a, b, c, l) { (l$x[i]/a)**2+(l$y[i]/b)**2+(l$z[i]/c)**2 }, a, b, c, candidates ))
		inside <- candidates[ocd <= 1,]
		if ( 0 != nrow(inside) ) points <- rbind(points, inside[1:(min(nrow(inside),N-nrow(points))),])
	}
	return(points)
}
getIntensityFromDistance = function(tc, pc, mean=0, sd=1) {
	# Calculates the intensity of a pixel based
	# on the distance of the pixel centre from
	# the true dot centre.
	# 
	# Args:
	# 	tc: true centre in cartesian coordinates (2D or 3D)
	# 	pc: pixel centre in cartesian coordinates (2D or 3D)
	# 	
	# Returns:
	# 	The intensity as the density of the normal distribution
	# 	based on the distance between true and pixel centre
	# 	
	return(dnorm(sqrt(sum((tc - pc)**2)), mean=mean, sd=sd))
}
getCamera = function(x, y, span, xlim=NULL, ylim=NULL) {
	# Provides the coordinates of the pixels in a given subset
	# around the provided pixel.
	# 
	# Args:
	# 	x, y: central pixel coordinates (integers)
	# 	span: maximum distance of the pixels to be included in the subset
	# 	[xlim, ylim: minimum and maximum x/y pixel coordinates]
	# 
	# Returns:
	# 	A data.frame with two columns: the coordinates of the
	# 	pixels in the subset.
	# 
	grid <- expand.grid((x-span):(x+span), (y-span):(y+span))
	grid <- as.data.frame(grid, stringsAsFactors=F)
	colnames(grid) <- c('x', 'y')
	if ( !is.null(xlim) ) {
		if ( 2 != length(xlim) ) return(NULL)
		remove <- which(grid[,1] < xlim[1] | grid[,1] > xlim[2])
		if ( 0 != length(remove) ) grid <- grid[-remove,]
	}
	if ( !is.null(ylim) ) {
		if ( 2 != length(ylim) ) return(NULL)
		remove <- which(grid[,2] < ylim[1] | grid[,2] > ylim[2])
		if ( 0 != length(remove) ) grid <- grid[-remove,]
	}
	return(grid)
}
calculateCameraIntensity = function(i, px, py, sh, ix, iy, points, span, mean=0, sd=1) {
	# Calculates the intensity of each pixel in the camera
	# The camera is a square with the pixel containing the point in the middle
	# and comprising pixels up to span pixels away from the middle
	# 
	# Args:
	# 	i: the id of the point
	# 	px, py, sh: voxel sizes
	# 	points: the list of points
	# 	span: half-width of the camera
	# 	mean, sd: mean and sd of fitted gaussian
	# 
	# Returns:
	# 	A data.frame containing the pixel coordinates and intensity (from dnorm)
	# 	i.e.: four columns: x, y, z, intensity
	# 
	point <- unlist(points[i,])
	pixel <- data.frame(ceiling(point[1] / px), ceiling(point[2] / py), ceiling(point[3] / sh), stringsAsFactors=F)
	colnames(pixel) <- c('x', 'y', 'z')
	submask <- getCamera(pixel$x, pixel$y, span, xlim=c(1, ix), ylim=c(1, iy))
	getSubmaskIntensity = function(i, submask) {
		sub <- c(unlist(submask[i,]), pixel$z)
		centre <- (sub - 0.5) * c(px, py, sh)
		intensity <- getIntensityFromDistance(point, centre, mean=mean, sd=sd)
		return(c(sub, intensity))
	}
	return(do.call(rbind, lapply(1:nrow(submask), FUN=getSubmaskIntensity, submask)))
}
sumDuplicatedPixels = function(pixels, cores=1) {
	# Sums the intensity of duplicated pixels in the pixel list
	# 
	# Args:
	# 	pixels: data.frame with a row per pixel and 4 columns (x, y, z, intensity)
	# 	cores: the number of cores for parallelized analysis
	#
	# Returns:
	# 	A data.frame with the same format as pixels
	#
	pid <- paste0(pixels$x, '~', pixels$y, '~', pixels$z)
	rid <- pid[duplicated(pid, fromLast=F) | duplicated(pid, fromLast=T)]
	cleanPixels = function(id, pid, pixels) {
		subset <- pixels[which(pid == id),]
		if ( 1 != length(unique(subset$intensity)) ) {
			#return(c(unlist(subset[1,1:3]), min(subset$intensity)))
			return(c(unlist(subset[1,1:3]), sum(subset$intensity)))
		} else {
			return(subset)
		}
	}
	clean <- do.call(rbind, mclapply(unique(rid),
		FUN=cleanPixels,
		pid, pixels,
		mc.cores=cores
	))
	clean <- as.data.frame(clean, stringsAsFactors=F)
	colnames(clean) <- c('x', 'y', 'z', 'intensity')
	clean <- rbind(clean, pixels[!pid %in% rid,])
	return(clean)
}
mkStack = function(ix, iy, ns, pixels) {
	# Fills a 3D array (stack) with the intensity of the provided pixels
	# 
	# Args:
	# 	ix, iy, ns: image dimensions
	# 	pixels: data.frame with a row per pixel and 5 columns
	# 			(x, y, z, intensity, norm_intensity)
	# 			norm_intensity is expected to have values in the [0;1] interval
	# 
	# Returns:
	# 	A 3D array with ix, iy, ns dimensions
	# 	Intensity goes from 0 to 1
	# 
	stack <- array(0, c(ix, iy, ns))
	for (i in 1:nrow(pixels)) {
		pixel <- unlist(pixels[i,])
		stack[pixel[1], pixel[2], pixel[3]] <- pixel[5]
	}
	return(stack)
}
run = function(N, a, b, c, px, py, sh, dr, decs=3, cores=1, fname="test.tif") {
	# Executes the code: simulates N dots in an ellipsoid volume of
	# given dimensions, with uniform distribution.
	# 
	# Args:
	# 	Description of the parameters are in the first part of this page.
	# 	
	# Returns:
	# 	A list with the parameters from the call, the calculated variables,
	# 	clean and raw pixel list and the stack.
	# 	Also a TIF image based on the stack is written as output.
	#
	cat('> Calculate variables.\n')
	ix <- ceiling(2 * a / px)		# x image size in pixels
	iy <- ceiling(2 * b / py)		# y image size in pixels
	ns <- ceiling(2 * c / sh)		# number of slices in the z-stack
	sd <- dr / 2.35482									# sd of the gaussian profile
	span <- ceiling(max(3 * sd / px, 3 * sd / py))		# half-size of the camera
	
	cat('> Retrieve points (origin in the center of the ellipsoid).\n')
	points <- uniformPoints.ellipsoid(N, a, b, c)

	cat('> Traslate origin.\n')
	points$x <- points$x + a
	points$y <- points$y + b
	points$z <- points$z + c

	# Produces a list of pixels with their intensity (span**2 pixels per point)
	cat('> Calculate pixel intensity.\n')
	pixels <- do.call(rbind, mclapply(1:nrow(points),
		FUN=calculateCameraIntensity,
		px, py, sh, ix, iy, points, span, sd=sd,
		mc.cores=cores
	))
	pixels <- as.data.frame(pixels, stringsAsFactors=F)
	colnames(pixels) <- c('x', 'y', 'z', 'intensity')

	# Clean duplicates
	cat('> Sum overlapping signals.\n')
	pixels2 <- sumDuplicatedPixels(pixels, cores=cores)

	cat('> Normalize intensity.\n')
	pixels2 <- cbind(pixels2, pixels2[,4] - min(pixels2[,4]) + 10**(-decs))
	pixels2[,5] <- pixels2[,5] / max(pixels2[,5])

	# Rename columns
	pixels2 <- as.data.frame(pixels2, stringsAsFactors=F)
	colnames(pixels2) <- c('x', 'y', 'z', 'intensity', 'norm_intensity')

	# Round intensity and remove those with intensity equal to 0
	cat('> Round intensity.\n')
	pixels2$norm_intensity <- round(pixels2$norm_intensity, decs)
	rid <- which(0 == pixels2$norm_intensity)
	if ( 0 != length(rid) ) pixels2 <- pixels2[-rid,]

	cat('> Prepare stack.\n')
	stack <- mkStack(ix, iy, ns, pixels2)

	cat('> Print TIFF image.\n')
	writeTIFF(lapply(1:ns, FUN=function(i) { stack[,,i] } ), fname, bits.per.sample=8)

	cat('~ END ~')
	return(list(
		pixel_list=pixels,
		pixel_list_cleaned=pixels2,
		params=list(
			N=N,
			a=a, b=b, c=c,
			px=px, py=py, sh=sh,
			dr=dr,
			decs=decs,
			cores=cores,
			fname=fname
		),
		stack=stack,
		vars=list(
			sd=sd,
			span=span,
			ix=ix,
			iy=iy,
			sh=sh
		)
	))
}

# RUN
#--------------------------------------------

system.time({ r <- run(N, a, b, c, px, py, sh, dr, decs=decs, cores=cores, fname=fname) })
