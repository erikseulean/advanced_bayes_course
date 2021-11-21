#' Sequentially generate draws from a Dirichlet process mixture model,
#' by showing step by step the iterations taken.
#' The plot is centered at 0, with x and y from -5 to 5.
#' The mixture draws the centres for clusters from a Normal distribution
#' with mean mu and standard deviation sigma_0
#' Additional to plotting the points, it also returns
#' the points sampled.
#'
#' Hit enter to keep drawing until max n or type "x" to exit.
#' @param n Number of observations.
#' @param alpha Alpha corresponding to GEM(alpha) used to draw the rho vector.
#' @param mu Mean of the Normal distribution used to draw the clusters.
#' @param sigma_0 Standard deviation of the Normal distribution used
#' to draw the points around the cluster centre.
#' @param sigma Standard deviation for cluster centers
#' @return Returns the n observations sampled from the DPMM distribution.
#' @examples
#' rDPM(n=30, alpha=3, mu=0, sigma_0=1.5, sigma=0.7)
#' @export
rDPM <- function(n, alpha, mu, sigma_0, sigma) {
    rho <- c()
    rhosum <- 0
    rhosum_window <- c(0)
    clusters <- c()
    cluster_assignments <- c()
    cluster_bars <- c()
    samples <- c()

    for (sample in 1: n) {
        u <- runif(1)

        while (rhosum < u) {
            # stick break a portion
            p <- rbeta(1, 1, alpha)
            rho_draw <- (1 - rhosum) * p

            # add the new sampled rho to the list of all rhos
            rho <- rbind(rho, rho_draw)

            # store the total sum
            rhosum <- rhosum + rho_draw

            # store the sequence of sums up to this point
            # to figure out which cluster we need to pick
            rhosum_window <- c(rhosum_window, rhosum)

            # add a new centre for the cluster, associated with this rho
            clusters <- rbind(clusters, rnorm(2, mu, sigma_0))

            cluster_bars <- c(cluster_bars, "grey")
        }

        # get the index of the cluster for this observation
        cluster_index <- max(which(rhosum_window < u))
        cluster_mean <- clusters[cluster_index, ]

        #draw an observation from this cluster
        x <- rnorm(2, mean = cluster_mean, sd = sigma)

        # store which cluster is every point assigned to
        cluster_assignments <- c(cluster_assignments, cluster_index)

        samples <- rbind(samples, x)

        # plot all the rhos that we sampled so far.
        barplot(rbind(rho, 1 - rhosum),
            beside = FALSE,
            horiz = TRUE,
            col = c(cluster_bars, "white"),  # remaining mass in (0,1)
            ylim = c(0, 1),
            width = 0.7,
            main = bquote(rho~"~GEM("~.(alpha)~")")
        )

        # plot which point we sampled on the GEM
        points(u, 1, pch = 25, col = "red", bg = "red")

        plot(x,
            pch = ".",
            xlim = c(-5, 5),
            ylim = c(-5, 5),
            main = paste("N = ", toString(n),
                ", #means = ", toString(length(rho)),
                ", #clusters = ", toString(
                    length(unique(cluster_assignments))), sep = ""
                )
        )

        points(clusters,
            pch = 15,
            col = "black"
        )

        points(samples,
            pch = 19,
            col = cluster_assignments
        )

        line <- readline()
        if (line == "x") {
            return()
        }
    }

    return(samples)
}