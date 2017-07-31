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

add_pi_lm_log <- function(tb, fit, alpha = 0.05, names = NULL, yhatName) {
 
    if (is.null(names)) {
        names[1] <- paste("LPB", alpha/2, sep = "")
        names[2] <- paste("UPB", 1 - alpha/2, sep = "")
    }
    if ((names[1] %in% colnames(tb))) {
        warning ("These PIs may have already been appended to your dataframe. Overwriting.")
    }
 
  out <- predict(fit, tb, interval = "prediction", level = 1 - alpha)
    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- exp(out[, 1])
  tb[[names[1]]] <- exp(out[, 2])
  tb[[names[2]]] <- exp(out[, 3])
  tibble::as_data_frame(tb)
  
}


