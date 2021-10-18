generate_dirichlet_clusters = function(a, K) {
    dev.new(width=8, height=3, unit="in")
    rho = rdirichlet(1,rep(a,K))
    rmatrix = matrix(as.numeric(rho), ncol=1)
    bar_colors = rep("grey",K)

    barplot(
        rmatrix, 
        horiz=TRUE, 
        xlim=c(0, 1),
        col=bar_colors,
        beside=FALSE,
        main=bquote(rho~"~Dirichlet"~"("~.(a)~",...,"~.(a)~")K="~.(K))
    )
}

generate_dirichlet_clusters_with_sampled_points = function(n, a, K) {
    dev.new(width=8, height=3, unit="in")
    rho = rdirichlet(1,rep(a,K))
    rmatrix = matrix(as.numeric(rho), ncol=1)
    bar_colors = rep("grey",K)
    rho_cum_sum = c(0,cumsum(rho))

    print("Newline to sample one point at a time")
    for (i in 1:n) {
        u = runif(n=1)
        index = max(which(rho_cum_sum < u))
        bar_colors[index] = "blue"

        barplot(
            rmatrix, 
            horiz=TRUE, 
            xlim=c(0, 1),
            col=bar_colors,
            beside=FALSE,
            main=bquote(rho~"~Dirichlet"~"("~.(a)~",...,"~.(a)~") K="~.(K)~", N="~.(i))
        )

        points(u, 1.4, pch=25, col="red", bg="red", cex=3.5)
        line <- readline()
        if(line == "x") {
            return ("done")
        }
    }

}