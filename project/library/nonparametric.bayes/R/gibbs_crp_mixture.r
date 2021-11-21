#' @importFrom mvtnorm dmvnorm
NULL

library(mvtnorm)
library(progress)

#' Gibbs sampling for the Chinese Restaurant Process
#' Implementation details can be found in the associated paper
#' The algorithm stops at every 1000th iteration and prints
#' the current cluster configuration.
#'
#' @param data A matrix of nx2 containing the datapoints
#' @param sd Prior standard deviation
#' @param initialisation Cluster initialisation for each datapoint.
#' Default initialisation is to set every point in the same cluster.
#' @param sigma0 Covariance matrix for the points. Default initialisation
#' is set to matrix(c(1, 0, 0, 1), mrow=2, byrow=TRUE)
#' @return Returns the cluster assignments after the last iteration.
#' Examples
#' cluster_datapoints(generate_split_data(350, 0.5)$x, sigma0=diag(3^2, 2))
#' cluster_datapoints(petal, sigma0=petal_sigma0)
#' cluster_datapoints(width, sigma0=width_sigma0)
#' cluster_datapoints(mixed, sigma0=mixed_sigma0)
#' @export
cluster_datapoints  <- function(
    data,
    sd=1,
    initialisation=rep(1, nrow(data)),
    sigma0=matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
) {

    # these are for color blind people (or at least I can see them),
    # plus a few others in case the number of clusters is high at any point
    colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
        "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#1f001f",
        "#02137d", "#4d000e", "#300730", "#760021", "#3811c6"
    )
    palette(sample(colors, replace = FALSE))

    alpha <- 0.1

    # number of datapoints
    N <- nrow(data)

    # prior covariance matrix for the clusters
    sigma <- diag(sd^2, 2)

    # solve will find the inverse (basically what we need here)
    precision <- solve(sigma)
    precision0 <- solve(sigma0)

    # priod mean
    mu0 <- matrix(c(0, 0), ncol = 2, byrow = TRUE)

    cluster_assignments <- initialisation
    # how many datapoints we have in each cluster
    # initially, 1 cluster will contain all the datapoints
    counts <- as.vector(table(cluster_assignments))

    # store the number of clusters that we have at any time
    nr_clusters <- length(counts)
    plot(data, col = cluster_assignments, pch = 19)
    print("Initial cluster configuration, press enter to continue")
    line <- readline()

    # create a progress bar as it might be unintuitive
    # that you have to wait
    pb <- progress_bar$new(
        format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull]",
        total = 10000,
        complete = "=",   # Completion bar character
        incomplete = "-", # Incomplete bar character
        current = ">",    # Current bar character
        clear = FALSE,    # If TRUE, clears the bar when finish
        width = 100)      # Width of the progress bar

    # start running the sampler.
    # setting max iterations to 10000, but always can stop early
    for (iter in 1:10000) {
        pb$tick()
        if (iter %% 1000 == 0) {
            print(iter)
        }
        # for each datapoint, remove it from the cluster
        # calculate the probability to be in each of the clusters
        # and the probability of starting a new cluster
        # assign the datapoint to a cluster, and update the
        # cluster counts
        for (datapoint in 1:N) {
            # find the cluster
            cluster <- cluster_assignments[datapoint]
            # remove it from the cluster
            counts[cluster] <- counts[cluster] - 1

            if (counts[cluster] == 0) {
                counts[cluster] <- counts[nr_clusters]
                location <- (cluster_assignments == nr_clusters)
                cluster_assignments[location] <- cluster
                counts <- counts[-nr_clusters]
                nr_clusters <- nr_clusters - 1
            }
            cluster_assignments[datapoint] <- -1

            # will contain the log probability for the datapoint
            # to be part of each cluster OR to start a new cluster
            log_probabilities <- rep(NA, nr_clusters + 1)

            for (cluster in 1:nr_clusters) {
                # this is found by using the conjugate prior for
                # the normal distribution
                cluster_precision <- precision0 + counts[cluster] * precision
                cluster_sigma <- solve(cluster_precision)
                datapoints_in_cluster <- which(cluster_assignments == cluster)
                if (length(datapoints_in_cluster) > 1) {
                    points_sum <- colSums(
                        data[cluster_assignments == cluster, ]
                    )
                }
                else {
                    points_sum <- data[cluster_assignments == cluster, ]
                }
                # found using the conjugate prior
                cluster_mean <- cluster_sigma %*% (
                    precision %*% points_sum + precision0 %*% t(mu0)
                )
                # note that we don't need the denominator here, since
                # it's common to all the clusters
                # alpha + N - 1 - we can just sum everything and
                # divide by that value
                log_probabilities[cluster] <- log(counts[cluster]) + dmvnorm(
                    data[datapoint, ],
                    mean = cluster_mean,
                    sigma = cluster_sigma + sigma,
                    log = TRUE
                )
            }

            # add the probability to start a new cluster
            log_probabilities[nr_clusters + 1] <- log(alpha) + dmvnorm(
                data[datapoint, ],
                mean = mu0,
                sigma = sigma0 + sigma,
                log = TRUE
            )

            # transform and normalise log probabilities
            probabilities <- exp(log_probabilities - max(log_probabilities))
            probabilities <- probabilities / sum(probabilities)

            #sample which cluster this point belongs to
            new_cluster <- sample(
                1: (nr_clusters + 1), 1, replace = TRUE, prob = probabilities
            )

            #check if the point is starting a new cluster
            if (new_cluster == nr_clusters + 1) {
                counts <- c(counts, 0)
                nr_clusters <- nr_clusters + 1
            }

            cluster_assignments[datapoint] <- new_cluster
            counts[new_cluster] <- counts[new_cluster] + 1
        }
        if (iter %% 1000 == 0) {
            plot(data, col = cluster_assignments, pch = 19)
            print("Enter to continue or x to stop")
            line <- readline()
            if (line == "x") {
                dev.off()
                return(cluster_assignments)
            }

        }
    }

}