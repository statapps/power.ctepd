\name{survexp.ca}
\alias{survexp.ca}

\title{
  Canadian Census Data for the Expected Survival and Person Years
}
\description{
  Data from Statistics Canada for the expected survival and person years functions from year 2010 to 2012.
}

\details{
  survexp.ca for total Canadian population, by age and sex, 2010 to 2012. The table contains the annual hazard rate for a matched subject from the population.
}

\source{
  https://www150.statcan.gc.ca/n1/pub/84-537-x/2018001/xls/2010-2012_Tbl-eng.xlsx
}

\examples{
  data(survexp.ca)
  head(survexp.ca)
  is.ratetable(survexp.ca)
}

