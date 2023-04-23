#' Simulated data for the illustration of using spatial-temporal negative-binomial 
#' regression models to estimate attributable deaths
#'
#' A dataset containing the simulated death counts and the percent positive across 
#' 51 US states and 93 epidemic weeks. The variables are as follows:
#'
#' @docType data
#' 
#' @usage data(toydata)
#' 
#' @format A data frame with 4743 rows and 6 variables:
#' \describe{
#'   \item{logPop}{log of populations in each state}
#'   \item{week}{index of epidemic weeks}
#'   \item{state}{index of states}
#'   \item{lag1}{simulated percent positive at lag 1}
#'   \item{lag2}{simulated percent positve at lag 2}
#'   \item{y}{observed death counts}
#' }
#' 
#' @examples 
#' data(toydata)
"toydata"