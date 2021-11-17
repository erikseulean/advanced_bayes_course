
# Generation of a Dirichlet distribution
# Copyright (C) 2021, Erik-Cristian Seulean

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



#' Generate a sample from a Dirichlet distirbution
#' Using:
#'https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
#' @param n Number of observations.
#' @param alpha A vector containing the parameters for
#' the Dirichlet distribution.
#' @return A sample of n observations from the Dirichlet distribution.
#' @examples
#' rdirichlet(n=1, alpha=c(2, 2, 2))
#' @export
rdirichlet <- function(n, alpha) {

    if (n < 0) {
        return()
    }

    v <- list()

    for (i in 1:n) {
        gammas <- sapply(alpha, function(k) {
            return(rgamma(n = 1, shape = k))
        })
        v[[i]] <- gammas / (sum(gammas))
    }
    return(if (length(v) == 1) v[[1]] else v)
}