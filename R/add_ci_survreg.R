# Copyright (C) 2017 Institute for Defense Analyses
#
# This file is part of ciTools.
#
# ciTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ciTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ciTools. If not, see <http://www.gnu.org/licenses/>.


## ' @param alpha A real number between 0 and 1. Controls the confidence
## '     level of the interval estimates.
## ' @param names \code{NULL} or character vector of length two. If
## '     \code{NULL}, confidence bounds automatically will be named by
## '     \code{add_ci}, otherwise, the lower confidence bound will be
## '     named \code{names[1]} and the upper confidence bound will be
## '     named \code{names[2]}.
## ' @param yhatName A string. Name of the vector of predictions made
## '     for each observation in tb.
## ' @param method A string. One of either "boot" or "parametic".
## ' @param ... Additional arguments.
## '
## ' @return A tibble, \code{tb}, with predicted values, upper and lower
## '     confidence bounds attached.
## '
## ' @seealso \code{\link{add_pi.glm}} for prediction intervals for
## '     \code{glm} objects, \code{\link{add_probs.glm}} for conditional
## '     probabilities of \code{glm} objects, and
## '     \code{\link{add_quantile.glm}} for response quantiles of
## '     \code{glm} objects.
## '
## ' @examples
## '
## ' @export

## add_ci.survreg <- function(tb, fit, p = 0.5, alpha = 0.05,
##                            names = NULL, yhatName = "pred",
##                            method = "boot", nSims = 2000, ...){

##     if (is.null(names)){
##         names[1] <- paste("LCB", alpha/2, sep = "")
##         names[2] <- paste("UCB", 1 - alpha/2, sep = "")
##     }
##     if ((names[1] %in% colnames(tb))) {
##         warning ("These CIs may have already been appended to your dataframe. Overwriting.")
##     }
## }
