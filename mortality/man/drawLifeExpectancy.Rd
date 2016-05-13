\name{drawLifeExpectancy}
\alias{drawLifeExpectancy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
drawLifeExpectancy(age)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{age}{
%%     ~~Describe \code{age} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (age) 
{
    qxStochasticSet <- readRDS("data/QxStochasticSet.rds")
    df <- data.frame(life_expectancy = eex(qxStochasticSet, age))
    ggplot(df, aes(x = life_expectancy)) + geom_histogram(aes(y = ..density..), 
        binwidth = 0.5, colour = "blue", fill = "white") + geom_density() + 
        geom_vline(aes(xintercept = mean(life_expectancy, na.rm = T)), 
            color = "red", linetype = "dashed", size = 1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
