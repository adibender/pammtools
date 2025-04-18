% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi-state-helpers.R
\name{add_counterfactual_transitions}
\alias{add_counterfactual_transitions}
\title{Add counterfactual observations for possible transitions}
\usage{
add_counterfactual_transitions(
  data,
  from_to_pairs = list(),
  from_col = "from",
  to_col = "to",
  transition_col = "transition"
)
}
\arguments{
\item{data}{Data set that only contains rows for transitions that took place.}

\item{from_to_pairs}{A list with one element for each possible initial state.
The values of each list element indicate possible transitions from that state.
Will be calculated from the data if unspecified.}

\item{from_col}{Name of the column that stores initial state.}

\item{to_col}{Name of the column that stores end state.}

\item{transition_col}{Name of the column that contains the transition identifier (factor variable).}
}
\description{
If data only contains one row per transition that took place, this function
adds additional rows for each transition that was possible at that time
(for each subject in the data).
}
