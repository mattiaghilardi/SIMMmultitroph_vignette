# LICENSE -----------------------------------------------------------------
# 
# MIT License
# 
# Copyright (c) 2026 Mattia Ghilardi
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# 
# 
# Function ---------------------------------------------------------------

#' Estimate error in trophic position
#'
#' @param sources A data frame with mean isotopic values of sources. 
#' Must have a column named "d15N_mean" and, if using two baselines,
#' a column named "d13C_mean".
#' @param baselines A data frame with mean isotopic values of one or two baselines. 
#' Must have a column named "d15N_mean" and, if using two baselines,
#' a column named "d13C_mean".
#' @param TP_b A numeric value. Trophic position of baseline(s). 
#' Default to 2 (i.e., primary consumers).
#' @param TDF_N A numeric value. The TDF for \eqn{\delta^{15}\text{N}}.
#' @param TDF_C A numeric value. The TDF for \eqn{\delta^{13}\text{C}}. Default to NULL. 
#' Only necessary for two baselines.
#' @param consumer_p A numeric vector summing to 1. 
#' The source relative contributions to the consumer.
#'
#' @returns A numeric value
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#' 
#' # 3 sources
#' n_sources <- 3
#' sources <- data.frame(source = paste("source", 1:n_sources),
#'                       d13C_mean = rnorm(n_sources, 0, 5), 
#'                       d15N_mean = rnorm(n_sources, 0, 5))
#' 
#' # TDFs (values from Post 2002 for this example)
#' TDF_C <- 0.4
#' TDF_N <- 3.4
#' 
#' # 2 baselines
#' n_baselines <- 2
#' 
#' # Random number generation from the dirichlet distribution
#' # n: number of random vectors to generate
#' # alpha: vector of positive shape parameters
#' rdirichlet <- function (n, alpha) {
#'     x <- matrix(rgamma(length(alpha) * n, alpha), 
#'                 ncol = length(alpha), 
#'                 byrow = TRUE)
#'     return(x / rowSums(x))
#' }
#' 
#' # Generate baselines isotopic values:
#' # 1. randomly generate source proportions
#' # 2. multiply by source values and summarise
#' # 3. add TDF to set baselines at TP = 2
#' baselines <- rdirichlet(n = n_baselines, alpha = c(1, 1, 1)) |> 
#'   t() |> 
#'   as.data.frame() |> 
#'   setNames(paste("baseline", 1:n_baselines)) |>
#'   mutate(source = paste("source", 1:n_sources)) |> 
#'   pivot_longer(cols = starts_with("baseline"),
#'                names_to = "baseline",
#'                values_to = "p") |>
#'   left_join(sources) |>
#'   mutate(across(ends_with("_mean"), ~.x * p)) |>
#'   group_by(baseline) |>
#'   summarise(d13C_mean = sum(d13C_mean) + TDF_C, 
#'             d15N_mean = sum(d15N_mean) + TDF_N)
#'
#' # Source relative contributions for 7 consumers
#' consumers_p <- list(c(0.8, 0.1, 0.1), # relies mostly on source 1
#'                     c(0.1, 0.8, 0.1), # relies mostly on source 2
#'                     c(0.1, 0.1, 0.8), # relies mostly on source 3
#'                     c(0.45, 0.45, 0.1), # relies mostly on source 1 and 2
#'                     c(0.45, 0.1, 0.45), # relies mostly on source 1 and 3
#'                     c(0.1, 0.45, 0.45), # relies mostly on source 2 and 3
#'                     c(1/3, 1/3, 1/3)) # relies equally on sources
#' names(consumers_p) <- paste("consumer", 1:7)
#' 
#' lapply(consumers_p,
#'        function(i) 
#'          estimate_tp_error(sources, baselines, TP_b = 2, TDF_N, TDF_C, i)
#' )
#' 
estimate_tp_error <- function(sources, baselines, TP_b = 2, TDF_N, TDF_C = NULL, consumer_p) {
  
  n_baselines <- nrow(baselines)
  
  if (!n_baselines %in% 1:2) 
    stop(paste("Please provide mean isotope values for only one or two baselines.",
               "\n\"baselines\" has", n_baselines, "rows."))
  
  if (!is.double(TDF_N) | length(TDF_N) != 1) 
    stop("\"TDF_N\" must be a single numeric value.")
  
  if (n_baselines == 2 & (!is.double(TDF_C) | length(TDF_C) != 1)) 
    stop("When using two baselines \"TDF_C\" must be a single numeric value.")
  
  if (n_baselines == 1 & !"d15N_mean" %in% colnames(sources))
    stop(paste("When using one baseline \"sources\" must be a data frame", 
               "including a column named \"d15N_mean\"."))
  
  if (n_baselines == 2 & !all(c("d13C_mean", "d15N_mean") %in% colnames(sources)))
    stop(paste("When using two baselines \"sources\" must be a data frame", 
               "including two columns named \"d13C_mean\" and \"d15N_mean\"."))
  
  if (n_baselines == 1 & !"d15N_mean" %in% colnames(baselines))
    stop(paste("When using one baseline \"baselines\" must be a data frame", 
               "including a column named \"d15N_mean\"."))
  
  if (n_baselines == 2 & !all(c("d13C_mean", "d15N_mean") %in% colnames(baselines)))
    stop(paste("When using two baselines \"baselines\" must be a data frame", 
               "including two columns named \"d13C_mean\" and \"d15N_mean\"."))
  
  if (!is.double(consumer_p) | length(consumer_p) != nrow(sources) | !isTRUE(all.equal(sum(consumer_p), 1)))
    stop(paste("\"consumer_p\" must be a numeric vector of the same length as the",
         "number of sources and should sum to 1."))
  
  # Compute consumer isotopic values at TP = 1 (consumer TP does not affect pb1 and TP error)
  consumer <- data.frame(d15N_mean = sum(sources$d15N_mean * consumer_p))
  if (n_baselines == 2) 
    consumer$d13C_mean = sum(sources$d13C_mean * consumer_p)
  
  # Compute error
  if (n_baselines == 1) {
    error <- (consumer$d15N_mean - (baselines$d15N_mean - TDF_N * (TP_b - 1))) / TDF_N
  } else {
    # Get p of baseline 1
    pb1 <- (consumer$d13C_mean - 
              baselines$d13C_mean[2] - 
              (TDF_C/TDF_N) * consumer$d15N_mean + 
              (TDF_C/TDF_N) * baselines$d15N_mean[2]) /
      (baselines$d13C_mean[1] - 
         baselines$d13C_mean[2] - 
         (TDF_C/TDF_N) * baselines$d15N_mean[1] + 
         (TDF_C/TDF_N) * baselines$d15N_mean[2])
    
    # Constrain pb1 between 0 and 1
    # when consumer mean_d13C is outside the d13C range of baselines
    pb1 <- if (pb1 > 1) 1 else if (pb1 < 0) 0 else pb1
    
    error <- (consumer$d15N_mean - 
                pb1 * (baselines$d15N_mean[1] - TDF_N * (TP_b - 1)) - 
                (1 - pb1) * (baselines$d15N_mean[2] - TDF_N * (TP_b - 1))) /
      TDF_N
  }
  
  return(error)
}
