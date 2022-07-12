#' Code from Wang et al (2021) and Wang et al (2022)
#'
#' Takes the code from Wang et al (2021) and Wang et al (2022)
#' and places it in package form so that it is easier to compare.
#' I would not recommend using this package for real work.
#' You should use the \href{https://cran.r-project.org/package=hwep}{\code{hwep}}
#' package from Gerard (2021).
#'
#' The original repos with the original code are
#' \url{https://github.com/CCBBeijing/hexaploid} and
#' \url{https://github.com/CCBBeijing/OctoploidDeer}.
#'
#' If those repos are ever deleted or made private, you can see my forks at
#' \url{https://github.com/dcgerard/hexaploid} and
#' \url{https://github.com/dcgerard/OctoploidDeer}.
#'
#' @keywords internal
#'
#' @importFrom stats pchisq
#' @importFrom stats D
#' @importFrom hwep rmlike
#'
#' @references
#' \itemize{
#'   \item{D. Gerard. Double reduction estimation and equilibrium tests in natural autopolyploid populations. \emph{bioRxiv}, 2021. \doi{10.1101/2021.09.24.461731}.}
#'   \item{J. Wang, L. Feng, S. Mu, A. Dong, J. Gan, Z. Wen, J. Meng, M. Li, R. Wu, and L. Sun. Asymptotic tests for Hardy-Weinberg equilibrium in hexaploids. \emph{Horticulture Research}, 9, 05 2022. \doi{10.1093/hr/uhac104}.}
#'   \item{J. Wang, X. Lv, L. Feng, A. Dong, D. Liang, and R. Wu. A tracing model for the evolutionary equilibrium of octoploids. \emph{Frontiers in Genetics}, 12, 2021. \doi{10.3389/fgene.2021.794907}.}
#' }
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
