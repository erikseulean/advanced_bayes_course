library(mvtnorm)

ex7_gen_data <- function(Ndata, sd) {
# generate Gaussian mixture model data for inference later
#
# Args:
#  Ndata: number of data points to generate
#  sd: covariance matrix of data points around the
#      cluster-specific mean is [sd^2, 0; 0, sd^2];
#      i.e. this is the standard deviation in either direction
#
# Returns:
#  x: an Ndata x 2 matrix of data points
#  z: an Ndata-long vector of cluster assignments
#  mu: a K x 2 matrix of cluster means,
#      where K is the number of clusters

	# matrix of cluster centers: one in each quadrant
	mu = matrix(c(3,3, -3,3, 3,-3, -3,-3), ncol=2, byrow=TRUE)
	# vector of component frequencies
	rho = c(0.5,0.3,0.2,0.1)

	# assign each data point to a component
	z = sample(1:length(rho), Ndata, replace=TRUE, prob=rho)
	# draw each data point according to the cluster-specific
	# likelihood of its component
	x = cbind(rnorm(rep(NA,Ndata), mu[z,1], rep(sd,Ndata)),
		rnorm(rep(NA,Ndata), mu[z,2], rep(sd,Ndata)))
	
	# return the data.
	# also return the cluster centers and means in case
	# that is useful for comparison
	list("x" = x, "z" = z, "mu" = mu)
}

cluster_datapoints  = function(data, sd, initialisation, sigma0) {

    # these are for color blind people, 
    # plus a few others in case the number of clusters is high at any point
    colors = c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#F79CD4",
        "#B5BBE3","#E6AFB9","#300730","#760021", "#3811c6"
    )
    palette(sample(colors, replace=FALSE))
    min_iterations = 0

    alpha = 0.1

    # number of datapoints
    N = nrow(data)

    # prior covariance matrix for the clusters
    sigma = diag(sd^2, 2)
    
    # precision = sigma^(-1)
    # solve will find the inverse (basically what we need here) 
    precision = solve(sigma)
    precision0 = solve(sigma0)

    # priod mean
    mu0 = matrix(c(5.8, 1.2), ncol=2, byrow=TRUE)

    cluster_assignments = initialisation
    # how many datapoints we have in each cluster
    # initially, 1 cluster will contain all the datapoints
    counts = as.vector(table(cluster_assignments))

    # store the number of clusters that we have at any time
    nr_clusters = length(counts)    
    plot(data, col=cluster_assignments, pch=19)
    line = readline()
    # running the sampler on 1000 iterations
    # will see if this needs to be bigger
    for(iter in 1:1000) {
        print(iter)
        # for each datapoint, remove it from the cluster
        # calculate the probability to be in each of the clusters
        # and the probability of starting a new cluster
        # assign the datapoint to a cluster, and update the
        # cluster counts
        for(datapoint in 1:N) {
            # find the cluster
            cluster = cluster_assignments[datapoint]
            # remove it from the cluster
            counts[cluster] = counts[cluster] - 1

            if(counts[cluster] == 0) {
                counts[cluster] = counts[nr_clusters]
                location = (cluster_assignments==nr_clusters)
                cluster_assignments[location] = cluster
                counts = counts[-nr_clusters]
                nr_clusters = nr_clusters - 1
            }            
            cluster_assignments[datapoint] = -1

            # will contain the log probability for the datapoint 
            # to be part of each cluster OR to start a new cluster
            log_probabilities = rep(NA, nr_clusters + 1)

            for(cluster in 1:nr_clusters) {
                # this is found by using the conjugate prior for the normal distribution                
                cluster_precision = precision0 + counts[cluster]*precision
                cluster_sigma = solve(cluster_precision)
                datapoints_in_cluster = which(cluster_assignments==cluster)
                if(length(datapoints_in_cluster) > 1) {                          
                    points_sum = colSums(data[cluster_assignments == cluster,])
                }
                else {
                    points_sum = data[cluster_assignments == cluster, ]
                    
                }
                # found using the conjugate prior                
                cluster_mean = cluster_sigma %*% (precision %*% points_sum + precision0 %*% t(mu0))
                # note that we don't need the denominator here, since it's common to all the clusters
                # alpha + N - 1 - we can just sum everything and divide by that value
                log_probabilities[cluster] = log(counts[cluster]) + dmvnorm(data[datapoint, ], mean=cluster_mean, sigma=cluster_sigma + sigma, log=TRUE)            
            }

            # add the probability to start a new cluster
            log_probabilities[nr_clusters + 1] = log(alpha) + dmvnorm(data[datapoint, ], mean=mu0, sigma=sigma0 + sigma, log=TRUE)

            # transform and normalise log probabilities
            probabilities = exp(log_probabilities - max(log_probabilities))
            probabilities = probabilities / sum(probabilities)

            #sample which cluster this point belongs to
            new_cluster = sample(1: (nr_clusters+1), 1, replace=TRUE, prob=probabilities)

            #check if the point is starting a new cluster
            if (new_cluster == nr_clusters + 1) {
                counts = c(counts, 0)
                nr_clusters = nr_clusters + 1
            }

            cluster_assignments[datapoint] = new_cluster
            counts[new_cluster] = counts[new_cluster] + 1
        
        plot(data, col=cluster_assignments, pch=19)
        # if(iter >= min_iterations) {
        #     plot(data, col=cluster_assignments, pch=19)
        # }

        # line <- readline()
        # if(line == "x") {
        # 	dev.off()
        # 	return("done")
        # } else if(line == "") {

        # } else {
        # 	min_iterations = iter + as.numeric(line)
        # }
    }
    }

}