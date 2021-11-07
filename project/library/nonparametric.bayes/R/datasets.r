
#' Generates a dataset used to exemplify clustering
#' The cluster centers are set relatively far away to
#' see how well the algorithm performs in simple scenarios
#' @param n Number of datapoints to generate
#' @param sd Standard deviation from the cluster center
#'
#' @return Returns the datapoints and the cluster assignments.
#' The cluster assignments can be used to calculate the performance
#' of the clustering.
generate_split_data <- function(n, sd) {
	# cluster centres
	mu <- matrix(c(3, 3, -3, 3, 3, -3, -3, -3, 0, 0), ncol = 2, byrow = TRUE)
	# how many points to consider in each cluster
	rho <- c(0.3, 0.3, 0.2, 0.1, 0.2)

	# assign each data point to a component
	z <- sample(1:length(rho), n, replace = TRUE, prob = rho)

    # draw each data point according to the cluster-specific
	# likelihood of its component
	x <- cbind(rnorm(rep(NA, n), mu[z, 1], rep(sd, n)),
		rnorm(rep(NA, n), mu[z, 2], rep(sd, n)))

	#return data and cluster to compute the accuracy
	list("x" = x, "z" = z)
}

normalise_iris <- function(dataset) {
    numeric_index <- sapply(dataset, is.numeric)
    dataset[numeric_index] <- lapply(dataset[numeric_index], scale)
    return(dataset)
}

data(iris)
iris <- normalise_iris(iris)
petal <- as.matrix(
	subset(iris, select = -c(Sepal.Length, Sepal.Width, Species))
)
width <- as.matrix(
	subset(iris, select = -c(Sepal.Length, Petal.Length, Species))
)
mixed <- as.matrix(
	subset(iris, select = -c(Sepal.Length, Petal.Width, Species))
)

petal_sigma0 <- matrix(
	c(var(petal[, 1]), 0, 0, var(petal[, 2])), byrow = TRUE, nrow = 2
)
width_sigma0 <- matrix(
	c(var(width[, 1]), 0, 0, var(width[, 2])), byrow = TRUE, nrow = 2
)
mixed_sigma0 <- matrix(
	c(var(mixed[, 1]), 0, 0, var(mixed[, 2])), byrow = TRUE, nrow = 2
)

split_data <- generate_split_data(350, 0.5)