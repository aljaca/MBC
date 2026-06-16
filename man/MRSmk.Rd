\name{MRSmk}
\alias{MRSmk}
\title{Multivariate Rescaling Using the Linear Monge-Kantorovich Map}
\description{
Applies a multivariate dependence adjustment based on a Gaussian optimal
transport, or linear Monge-Kantorovich, map. Variables are first transformed to
van der Waerden normal scores. Correlation changes are represented either in the
log-Euclidean geometry of positive-definite correlation matrices or by additive
changes in Fisher-z transformed correlations. Separate maps are then applied to
current- and projected-period adjusted simulations, and the original marginal
values are restored by rank reordering.
}
\usage{
MRSmk(o.c, m.c, m.p, mhat.c, mhat.p,
      cor_change = c("logdelta", "fisherz"), ties.method = "random")
}
\arguments{
  \item{o.c}{Matrix of observed current-period data. Rows are cases and columns are variables.}
  \item{m.c}{Matrix of modelled current-period data.}
  \item{m.p}{Matrix of modelled projected-period data.}
  \item{mhat.c}{Matrix of already marginally adjusted current-period model data.}
  \item{mhat.p}{Matrix of already marginally adjusted projected-period model data.}
  \item{cor_change}{Character string specifying how the model-projected change in dependence is transferred. The default, \code{"logdelta"}, applies the change in the log-Euclidean geometry of correlation matrices. The alternative, \code{"fisherz"}, applies the change in Fisher-z transformed pairwise correlations.}
  \item{ties.method}{Method passed to \code{rank} when computing van der Waerden normal scores.}
}
\details{
\code{MRSmk} is designed to be applied after marginal bias adjustment, for example
after \code{QDM} or another univariate adjustment has produced \code{mhat.c} and
\code{mhat.p}. The observed current-period dependence is estimated from
\code{o.c}. The model-projected change in dependence between \code{m.c} and
\code{m.p} is used to define a projected target correlation matrix for the
adjusted future simulation. The linear Monge-Kantorovich map then transforms the
normal-score representation of \code{mhat.c} and \code{mhat.p} to the target
correlation matrices. Finally, each transformed normal-score column is used as a
rank template to reorder the original adjusted values, preserving their marginal
distributions.

All inputs must be numeric matrices with the same number of columns and no
missing values. The method requires at least two variables.
}
\value{
A list with components:
\item{mhat.c}{Current-period adjusted matrix after dependence rescaling.}
\item{mhat.p}{Projected-period adjusted matrix after dependence rescaling.}
}

