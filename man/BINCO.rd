\name{BINCO}
\alias{BINCO}
\title{Variable selection based on their selection frequencies}
\description{A function to calculate the optimal threshold of the selection 
            frequency for variable selection by directly estimating the false
            discovery rate (FDR).}

\usage{
BINCO(count.mix=NULL, freq.m=NULL, nb=NULL, FDR=0.05, vpr=0.8, conservative=F, 
      niv=3, ini.bound=c(10,10,20))
      }

\arguments{
    \item{count.mix}{a vector of non-negative integers. Its length is the number 
                    of models to be aggregated (e.g. number of boostrap resamples)
                    from which the selection frequencies are generated. The 
                    \code{i}th element represents the number of edges/variables 
                    being selected for \code{i} times across all models to be 
                    aggregated. For example, if the second entry is 50 then it 
                    means 50 edges/variables are selected exactly twice across 
                    all the models. The last entry is the number of 
                    edges/variables that are selected consistently by all models.
                    Note, \code{count.mix} and \code{freq.m} cannot both be NULL.}
                    
    \item{freq.m}{This input is specifically for edge selection in network 
                  inference. It is a (symmetric) matrix recording the selection 
                  frequencies of each edge across all the models, e.g., 
                  \code{freq.m[i,j]} is the selection frequency of the edge 
                  connecting nodes i and j. Note, \code{count.mix} and 
                  \code{freq.m} cannot both be NULL.}
                  
    \item{nb}{a positive integer as the number of models from which the selection
             frequencies are generated. It must be provided along with freq.m if
             count.mix is NULL.}
              
    \item{FDR}{a numeric value between 0 and 1. The desired control level for 
                the false discovery rate. Default is 0.05.}
                
    \item{vpr}{a numeric value which is recommended to be between 0.8 and 0.95, 
              as the rule for the "valley point" of the empirical distribution 
              (see Li, et al.,2012). Default value is 0.8. \code{BINCO} will not
              apply if the calculated "valley point" value is greater than 
              \code{vpr} or is greater than 0.95, because a large "valley point" 
              value may result in a failure of FDR control.}
                
    \item{conservative}{a logic value of TRUE or FALSE. Default is FALSE. Set 
                      conservative=TRUE if the FDR control needs to be 
                      conservative. Under the conservative mode, the count values
                      of the fitted null distribution beyond the "valley point" 
                      are set as a constant.} 
                      
    \item{niv}{a positive integer as the number of sets of initial parameter 
                values to be used for density fitting. Default is 3. In our 
                experience, 10 is large enough for most usual situations.}
               
    \item{ini.bound}{a vector of three positive real numbers. It gives the upper 
                      bound of randomly generated initial parameter values 
                      (a, b, r) of the powered-beta probability density function
                      which is used for density fitting. Default value is 
                      (10,10,20).}
           }           
\details{\code{BINCO} conducts model selection by directly controlling the FDRs 
        of the selected edges/variables. It assumes a mixture model for the 
        distribution of edge/variable selection frequencies obtained from model
        aggregation. To estimated the FDRs, it uses a convolution of powered 
        beta and binomial probability density function to fit the selection 
        frequency distribution of null edges/irrelevant variables. Based on the 
        estimated FDRs, an optimal cutoff value for selection frequencies is 
        calculated, such that the set of edges/variables with selection
        frequencies greater than or equal to the cutoff value has the largest 
        estimated power and its estimated FDR is smaller than or equally to the 
        pre-specified level. Details can be found in Li, et al., (2012).}
\value{
    A list of six components
    \item{cut_off}{a numeric value. It is the optimal cutoff value for the 
                  selection frequencies calculated by \code{BINCO}.}
                  
    \item{estimated_FDR}{a numeric value. It is the estimated FDR for the selection
                        of edges/variables with selection frequencies greater than
                        or equal to \code{cut_off}.}
                          
    \item{estimated_power}{a numeric value. It is the estimated number of true 
                          edges/variables among the edges/variables with selection 
                          frequencies greater than or equal to \code{cut_off}.}
                                               
    \item{fitting range}{a vector of length 2. It gives the (selection frequency)
                         range that is used for fitting the null distribution.}
                          
    \item{empirical_mix}{a vector of length \code{nb}. It contains non-negative
                             integers as the counts of edges/variables at different
                             selection frequencies. This is the same as 
                             \code{count.mix}, which will be calculated if 
                             it is not provided by the input.} 
                             
    \item{estimated_null}{a vector of length \code{nb}. It contains non-negative 
                          integers as the estimated counts of null edges/variables
                          at different selection frequencies.}                    
    }
      
\references{S. Li, L. Hsu, J. Peng and P. Wang (2012) Bootstrap Inference for 
            Network Construction. \code{http://arxiv.org/abs/1111.5028}}
            
\author{S. Li, L. Hsu, J. Peng and P. Wang}

\keyword{methods}

\examples{
################################################################################
###### BINCO applies on selection frequency data (generated using space ########
################################################################################

# load the selection frequency data, an empirical edge count distribution (Y) 
# generated by space (see Peng et al., 2009 for details)
data(BINCOSimulation)

#### 1. run BINCO under default settings
out=BINCO(count.mix=BINCOSimulation$Y)
BINCO.plot(out.BINCO=out)

out=BINCO(freq.m=BINCOSimulation$Y.m,nb=BINCOSimulation$nb)
BINCO.plot(out.BINCO=out)

nb=length(out$empirical_mix)
pick=((1:nb)/nb)>=out$cut_off
writeLines(paste("\n  Number of selected edges: ", sum(out$empirical_mix[pick]), "\n",
" Estimated FDR: ", round(out$estimated_FDR, 3), "\n",
" Estimated number of true edges: ", out$estimated_power, "\n"))

#### 2. run BINCO under the conservative mode
out=BINCO(count.mix=BINCOSimulation$Y, conservative=TRUE) 
BINCO.plot(out.BINCO=out)

nb=length(out$empirical_mix)
pick=((1:nb)/nb)>=out$cut_off
writeLines(paste("\n  Number of selected edges: ", sum(out$empirical_mix[pick]), "\n",
" Estimated FDR: ", round(out$estimated_FDR, 3), "\n",
" Estimated number of true edges: ", out$estimated_power, "\n"))
}

