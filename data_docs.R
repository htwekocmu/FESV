#' Developing-Countries Panel (2000-2023)
#'
#' Balanced panel of 25 developing countries observed over 24 years
#' (2000-2023), used in the empirical application of the FESV paper. Variables
#' follow the definitions in the paper; the response is \code{Income} and the
#' six candidate covariates are \code{FDI}, \code{INF}, \code{POP}, \code{TOP},
#' \code{UR}, and \code{EXR}.
#'
#' @format A data frame with 600 rows and 8 variables:
#' \describe{
#'   \item{id}{country index (1-25).}
#'   \item{year}{year (2000-2023).}
#'   \item{Income}{log GDP per capita.}
#'   \item{FDI}{net foreign direct investment inflows, \% of GDP.}
#'   \item{INF}{inflation rate, annual \%.}
#'   \item{POP}{population growth rate, annual \%.}
#'   \item{TOP}{trade openness (exports + imports), \% of GDP.}
#'   \item{UR}{urbanisation rate, \%.}
#'   \item{EXR}{exchange-rate volatility measure.}
#' }
#'
#' @source World Bank World Development Indicators, 2000-2023.
"developing_panel"

#' Developed-Countries Panel (2000-2023)
#'
#' Balanced panel of 20 developed countries observed over 24 years
#' (2000-2023). Same variable definitions as \code{\link{developing_panel}}.
#'
#' @format A data frame with 480 rows and 8 variables.
#' @source World Bank World Development Indicators, 2000-2023.
"developed_panel"

#' Pooled Panel of All Countries (2000-2023)
#'
#' Pooled panel of 45 countries (25 developing + 20 developed) observed over
#' 24 years (2000-2023). Same variable definitions as
#' \code{\link{developing_panel}}.
#'
#' @format A data frame with 1080 rows and 8 variables.
#' @source World Bank World Development Indicators, 2000-2023.
"all_panel"
