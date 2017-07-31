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


add_quantile_lm_log <- function(tb, fit, p, name = NULL, yhatName) {
    if (p <= 0 || p >= 1)
        stop ("p should be in (0,1)")
    if (is.null(name)) {
        name <- paste("quantile", p, sep = "")
    }

    if ((name %in% colnames(tb))) {
        warning ("These quantiles may have already been appended to your dataframe. Overwriting.")
    }
    
    out <- predict(fit, tb, interval = "prediction", se.fit = TRUE)
    fitted <- out$fit[,1]
    residual_df <- out$df
    se_fitted <- out$se.fit
    resid_var <- out$residual.scale^2
    se_pred <- sqrt(resid_var + se_fitted^2)
    t_quantile <- qt(p = p, df = residual_df)
    out_quantiles <- exp(fitted + se_pred * t_quantile)

    if(is.null(tb[[yhatName]]))
        tb[[yhatName]] <- exp(fitted)
    tb[[name]] <- out_quantiles
    tibble::as_data_frame(tb)
}


