# Copyright 2014-2020 Google Inc. All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# The reporting code below was taken from CausalImpact and adapted for use in the tscausal package.

report_statement <- function(actual, prediction, object, alpha, p, sig, pos, ci_coverage, prediction_lower, prediction_upper, absolute_effect, 
                              relative_effect, relative_effect_lower, relative_effect_upper, 
                              absolute_effect_lower, absolute_effect_upper)
{
    stmt <- NULL
    stmt <- paste0(stmt, "\n\nDuring the post-intervention period, the response ",
                   "variable had an average value of approx. ", actual[1],
                   ". ", if (sig) "By contrast, in " else "In ",
                   "the absence of an intervention, ",
                   "we would have expected an average response of ",
                   prediction[1], ". The ", ci_coverage, " interval of this ",
                   "counterfactual prediction is [", prediction_lower[1],
                   ", ", prediction_upper[1], "]. Subtracting this ",
                   "prediction from the observed response yields an estimate ",
                   "of the causal effect the intervention had on the response ",
                   "variable. This effect is ", absolute_effect[1], " with a ",
                   ci_coverage, " interval of [", absolute_effect_lower[1], ", ",
                   absolute_effect_upper[1], "]. For a discussion of ",
                   "the significance of this effect, see below.")
    if (object$include_cumulative) {
        # Summarize sums
        stmt <- paste0(stmt, "\n\nSumming up the individual data points during ",
                       "the post-intervention period (which can only sometimes be ",
                       "meaningfully interpreted), the response variable had an ",
                       "overall value of ", actual[2], ". ",
                       if (sig) "By contrast, had " else "Had ",
                       "the intervention not taken place, we would have expected ",
                       "a sum of ", prediction[2], ". The ", ci_coverage, " interval of ",
                       "this prediction is [", prediction_lower[2], ", ", prediction_upper[2],
                       "].")
    }
    # Summarize relative numbers (in which case row [1] = row [2])
    stmt <- paste0(stmt, "\n\nThe above results are given in terms of ",
                   "absolute numbers. In relative terms, the response variable ",
                   "showed ", if (pos) "an increase of " else "a decrease of",
                   relative_effect[1], ". The ", ci_coverage, " interval of this ",
                   "percentage is [", relative_effect_lower[1], ", ",
                   relative_effect_upper[1], "].")
    # Comment on significance
    if (sig && pos) {
        stmt <- paste0(stmt, "\n\nThis means that the positive effect observed ",
                       "during the intervention period is statistically ",
                       "significant and unlikely to be due to random ",
                       "fluctuations. ",
                       "It should be noted, however, that the question of whether ",
                       "this increase also bears substantive significance can ",
                       "only be answered by comparing the absolute effect (",
                       absolute_effect[1], ") to the original goal of ",
                       "the underlying intervention.")
    } else if (sig && !pos) {
        stmt <- paste0(stmt, "\n\nThis means that the negative effect observed ",
                       "during the intervention period is statistically ",
                       "significant. If the experimenter had expected a positive ",
                       "effect, it is recommended to ",
                       "double-check whether anomalies in the control variables ",
                       "may have caused an overly optimistic expectation of ",
                       "what should have happened in the response variable in the ",
                       "absence of the intervention.")
    } else if (!sig && pos) {
        stmt <- paste0(stmt, "\n\nThis means that, although the intervention ",
                       "appears to have caused a positive effect, this effect ",
                       "is not statistically significant when considering the ",
                       "entire post-intervention period as a whole. Individual ",
                       "days or shorter stretches within the intervention period ",
                       "may of course still have had a significant effect, as ",
                       "indicated whenever the lower limit of the impact ",
                       "time series (lower plot) was above zero.")
    } else if (!sig && !pos) {
        stmt <- paste0(stmt, "\n\nThis means that, although it may look as ",
                       "though the intervention has exerted a negative effect ",
                       "on the response variable when considering the ",
                       "intervention period as a whole, this effect is not ",
                       "statistically significant, and so cannot be ",
                       "meaningfully interpreted.")
    }
    if (!sig) {
        stmt <- paste0(stmt, " The apparent effect could be the result of ",
                       "random fluctuations that are unrelated to the ",
                       "intervention. This is often the case when the ",
                       "intervention period is very long and includes much ",
                       "of the time when the effect has already worn off. ",
                       "It can also be the case when the intervention period ",
                       "is too short to distinguish the signal from the noise. ",
                       "Finally, failing to find a significant effect can ",
                       "happen when there are not enough control variables or ",
                       "when these variables do not correlate well with ",
                       "the response variable during the learning period.")
    }
    if (p < alpha[1]) {
        stmt <- paste0(stmt, "\n\nThe probability of obtaining this effect by ",
                       "chance is very small (Bayesian one-sided tail-area ",
                       "probability p = ", round(p, 3), "). This means the causal ",
                       "effect can be considered statistically significant.")
    } else {
        stmt <- paste0(stmt, "\n\nThe probability of obtaining this ",
                       "effect by chance is p = ", round(p, 3), ". This ",
                       "means the effect may be spurious and would generally ",
                       "not be considered statistically significant.")
    }
    return(stmt)
}