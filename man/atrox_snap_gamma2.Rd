\name{atrox_snap_gamma2}
\alias{ExperimentalParameters}
\alias{GppFst_ExampleResults}
\alias{FDR_result}
\alias{power_results}
\docType{data}
\title{included data}
\description{These files contain data associated with the package GppFst. each file is available as a processed R data file or as a raw file

To access the Rdata version:

\code{data(atrox_snap_gamma2)}
\code{data(ExperimentalParameters)}
\code{data(GppFst_ExampleResults)}
\code{data(FDR_result)}
\code{data(power_results)}

To access the raw version use the \code{system.file} function.

\code{file.loc <- system.file("atrox_snap_gamma2.txt", package="GppFst")}
\code{file.loc <- system.file("ExperimentalParameters", package="GppFst"))}
\code{file.loc <- system.file("GppFst_ExampleResults", package="GppFst"))}
\code{file.loc <- system.file("FDR_result", package="GppFst"))}
\code{file.loc <- system.file("power_results", package="GppFst"))}

Then read the table in with the read.table function
\code{read.table(file.loc, header=T)}


}

\author{Rich Adams}
\keyword{datasets}
