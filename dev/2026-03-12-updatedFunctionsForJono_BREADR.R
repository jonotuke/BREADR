# List of changes:
# 1. Need to make the version for BREADR not so new.
# 2. Removed commented code from callRelatedness.
# 3. Made prior sum to one a warning and less restrictive for callRelatedness.
# 4. Made sure column names were in input data for functions, not exact.
# 5. Added "by population" option to processEigenstrat to make it faster.
# 6. Added a new plot called plotROLL to compare PMR values between groups.

#' callRelatedness
#'
#'
#' A function that takes PMR observations, and (given a prior distribution for
#' degrees of relatedness) returns the posterior probabilities of all pairs of
#' individuals being (a) the same individual/twins, (b) first-degree related,
#' (c) second-degree related or (d) "unrelated" (third-degree or higher).
#' The highest posterior probability degree of relatedness is also returned as
#' a hard classification.
#' Options include setting the background relatedness (or using the sample
#' median), a minimum number of overlapping SNPs if one uses the sample median
#' for background relatedness, and a minimum number of overlapping SNPs
#' for including pairs in the analysis.
#'
#' @param pmr_tibble a tibble that is the output of the processEigenstrat
#' function.
#' @param class_prior the prior probabilities for same/twin, 1st-degree,
#' 2nd-degree, unrelated, respectively.
#' @param average_relatedness a single numeric value, or a vector of numeric
#' values, to use as the average background relatedness. If NULL, the sample
#' median is used.
#' @param median_co if average_relatedness is left NULL, then the minimum cutoff
#' for the number of overlapping snps to be included in the median calculation
#' is 500.
#' @param filter_n the minimum number of overlapping SNPs for which pairs are
#' removed from the entire analysis. If NULL, default is 1.
#'
#' @return   results_tibble: A tibble containing 13 columns:
#' - row: The row number
#' - pair: the pair of individuals that are compared.
#' - relationship: the highest posterior probability estimate of the degree of relatedness.
#' - pmr: the pairwise mismatch rate (mismatch/nsnps).
#' - sd: the estimated standard deviation of the pmr.
#' - mismatch: the number of sites which did not match for each pair.
#' - nsnps: the number of overlapping snps that were compared for each pair.
#' - ave_re;: the value for the background relatedness used for normalisation.
#' - Same_Twins: the posterior probability associated with a same individual/twins classification.
#' - First_Degree: the posterior probability associated with a first-degree classification.
#' - Second_Degree: the posterior probability associated with a second-degree classification.
#' - Unrelated: the posterior probability associated with an unrelated classification.
#' - BF: A strength of confidence in the Bayes Factor associated with the highest posterior probability classification compared to the 2nd highest. (No longer included)
#' @export
#'
#' @examples
#' callRelatedness(counts_example,
#'   class_prior=rep(0.25,4),
#'   average_relatedness=NULL,
#'   median_co=5e2,filter_n=1
#' )
callRelatedness <- function(pmr_tibble,
                            class_prior=rep(0.25,4),
                            average_relatedness=NULL,
                            median_co=5e2,filter_n=1
){
  
  # Check pmr_tibble is of the form required:
  cols_check <- c('pair','nsnps','mismatch','pmr')
  if(!(all(names(pmr_tibble)%in%cols_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(cols_check,collapse=', '),'.'))
  }
  
  # Check that the prior distribution makes sense:
  if(any(class_prior<0)|any(class_prior>1)){
    stop('Posterior probabilities for degrees of relatedness must be values between 0 and 1.')
  }else if(abs(sum(class_prior)-1)>1e-3){
    warning('Posterior probabilities for degrees of relatedness must sum to 1. Normalising.')
    class_prior <- class_prior/sum(class_prior)
  }
  
  # Check that the average relatedness makes sense:
  if(!is.null(average_relatedness)){
    if(any(average_relatedness<=0) | any(average_relatedness>=1)){
      stop('The average relatedness must be a value between 0 and 1.')
    }
  }
  
  # Check median estimator cut off:
  if((median_co<=0)){
    stop('The median cut off number of overlapping snps must be greater than 0.')
  }else if(median_co>=max(pmr_tibble$nsnps,na.rm=T)){
    stop('The median cut off number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }
  
  # Check median estimator cut off:
  if((filter_n<=0)){
    stop('The cut off for the number of overlapping snps must be greater than 0.')
  }else if(filter_n>=max(pmr_tibble$nsnps,na.rm=T)){
    stop('The cut off for the number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }
  
  # local variables
  class_vec <- c('Same_Twins','First_Degree','Second_Degree','Unrelated')
  
  # If no user-defined value for average_relatedness is given, use the median (above nsnps cut off)
  if(is.null(average_relatedness)){
    M <- pmr_tibble %>%
      dplyr::filter(nsnps>median_co) %>%
      dplyr::pull(pmr) %>%
      stats::median(na.rm=T)
    filter_tibble <- pmr_tibble %>%
      dplyr::filter((nsnps>=filter_n)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(ave_rel=M)
  }else{
    M <- average_relatedness
    filter_tibble <- pmr_tibble %>%
      dplyr::mutate(ave_rel=M) %>%
      dplyr::filter((nsnps>=filter_n))
  }
  
  results_tibble <- filter_tibble %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Same_Twins=(stats::dbinom(mismatch,nsnps,0.5*ave_rel,log=T)+log(class_prior[1])),
                  First_Degree=(stats::dbinom(mismatch,nsnps,(1-0.5^2)*ave_rel,log=T)+log(class_prior[2])),
                  Second_Degree=(stats::dbinom(mismatch,nsnps,(1-0.5^3)*ave_rel,log=T)+log(class_prior[3])),
                  Unrelated=weightedBinom(mismatch,nsnps,ave_rel,log=T)+log(class_prior[4])) %>%
    dplyr::mutate(normConst=matrixStats::logSumExp(c(Same_Twins,First_Degree,Second_Degree,Unrelated))%>%exp()) %>%
    dplyr::mutate(relationship=class_vec[which.max(c(Same_Twins,First_Degree,Second_Degree,Unrelated))]) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Same_Twins=exp(Same_Twins)/normConst,
                  First_Degree=exp(First_Degree)/normConst,
                  Second_Degree=exp(Second_Degree)/normConst,
                  Unrelated=exp(Unrelated)/normConst) %>%
    dplyr::mutate(sd=sqrt(pmr*(1-pmr)/nsnps),
                  relationship=factor(relationship,levels=c('Same_Twins','First_Degree','Second_Degree','Unrelated')),
                  row=1:dplyr::n()) %>%
    dplyr::select(row,pair,relationship,pmr,sd,mismatch,nsnps,ave_rel,everything(),-normConst)
  return(results_tibble)
}

#' plotLOAF
#'
#' Plots all (sorted by increasing value) observed PMR values with
#' maximum posterior probability classifications represented by colour
#' and shape.
#' Options include a cut off for the minimum number of overlapping SNPs,
#' the max number of pairs to plot and x-axis font size.
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param nsnps_cutoff the minimum number of overlapping SNPs for which
#' pairs are removed from the plot. If NULL, default is 500.
#' @param N the number of (sorted by increasing PMR) pairs to plot.
#' Avoids plotting all pairs (many of which are unrelated).
#' @param fntsize the fontsize for the x-axis names.
#' @param verbose if TRUE, then information about the plotting process is sent
#' to the console
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' relatedness_example
#' plotLOAF(relatedness_example)
plotLOAF <- function(in_tibble,nsnps_cutoff=NULL,N=NULL,fntsize=7,verbose=TRUE){
  
  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if((!all(names(in_tibble)%in%col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }
  
  # Test that the overlapping snp cut off is sensible.
  if(is.null(nsnps_cutoff)){
    nsnps_cutoff <- 5e2
    if(verbose){
      cat(sprintf('No minimum number of overlapping SNPs given.\nUsing default minimum of 500.\n'))
    }
  }
  if((nsnps_cutoff<=0)){
    stop('The cut off for the number of overlapping snps must be greater than 0.')
  }else if(nsnps_cutoff>=max(in_tibble$nsnps,na.rm=T)){
    stop('The cut off for the number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }
  
  # test that the number of pairs to plot is sensible.
  if(is.null(N)){
    N <- min(nrow(in_tibble),50)
    if(verbose){
      cat(sprintf('No upper limit on number of pairs to plot given.\nPlotting first %i pairs.\n',N))
    }
  }else if(N>sum(in_tibble$nsnps>nsnps_cutoff)){
    N <- sum(in_tibble$nsnps>nsnps_cutoff)
    if(verbose){
      cat(sprintf('Number of pairs to plot was greater than the number of available pairs.\nPlotting first %i pairs.\n',N))
    }
  }
  
  # Define class breaks
  class_breaks <- tibble::tibble(relationship=factor(c('Same_Twins','First_Degree','Second_Degree','Unrelated'),
                                                     levels=c('Same_Twins','First_Degree','Second_Degree','Unrelated')),
                                 meanPMR=((1-0.5^c(1:3,Inf)))*in_tibble$ave_rel[1],
                                 col=ggcolorhue(4))
  # Produce plot
  gg2 <- in_tibble %>%
    dplyr::arrange(pmr) %>%
    dplyr::mutate(relationship=factor(relationship,levels=c('Same_Twins','First_Degree','Second_Degree','Unrelated'))) %>%
    dplyr::mutate(ymin=pmr-2*sd,ymax=pmr+2*sd) %>%
    dplyr::filter((nsnps>=nsnps_cutoff)) %>%
    dplyr::slice(1:N) %>%
    ggplot2::ggplot(ggplot2::aes(x=forcats::fct_reorder(pair,pmr),y=pmr))+
    ggplot2::geom_hline(ggplot2::aes(yintercept=meanPMR,col=relationship),linetype='dashed',data=class_breaks,linewidth=1)+
    ggplot2::theme_bw()+
    ggplot2::geom_errorbar(ggplot2::aes(forcats::fct_reorder(pair,pmr),ymin=ymin,ymax=ymax))+
    ggplot2::geom_point(ggplot2::aes(fill=relationship,shape=relationship),size=2, show.legend = TRUE)+
    ggplot2::ylab("PMR")+
    ggplot2::xlab(" ")+
    # scale_y_continuous(limits=c(0,0.25))+
    ggplot2::scale_fill_discrete(drop=FALSE,name=NULL)+
    ggplot2::scale_colour_discrete(drop=FALSE,name=NULL,guide='none')+
    ggplot2::scale_shape_manual(drop=FALSE,name=NULL,values=c(21:24))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1, size=fntsize),legend.position='bottom')+
    ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(size=5)))
  
  gg2 %>%
    return()
}

#' plotSLICE
#'
#' A function for plotting the diagnostic information when classifying a
#' specific pair (defined by the row number or pair name) of individuals.
#' Output includes the PDFs for each degree of relatedness (given the number of
#' overlapping SNPs) in panel A, and the normalised posterior probabilities
#' for each possible degree of relatedness.
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param row either the row number or pair name for which the
#' posterior distribution is to be plotted.
#' @param title an optional title for the plot. If NULL, the pair
#' from the user-defined row is used.
#' @param class_prior the prior probabilities for same/twin, 1st-degree,
#' 2nd-degree, unrelated, respectively.
#' @param showPlot If TRUE, display plot. If FALSE, just pass plot as a variable.
#' @param which_plot if 1, returns just the plot of the posterior distributions,
#' if 2 returns just the normalised posterior values. Anything else returns both
#' plots.
#' @param labels a length two character vector of labels for plots. Default is
#' no labels.
#'
#' @return a two-panel diagnostic ggplot object
#' @export
#'
#' @examples
#' plotSLICE(relatedness_example, row = 1)
plotSLICE <- function(
    in_tibble,
    row,title=NULL,
    class_prior=rep(1/4,4),
    showPlot=TRUE,
    which_plot = 0,
    labels = NULL){
  
  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if(!all(names(in_tibble%in%col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }
  
  # check row is a numeric or char variable, and if char, find the row.
  if(is.character(row)){
    rowtmp <- row
    row <- which(row==in_tibble$pair)
    if(length(row)==0){
      stop(paste0('There is no pair called: ',rowtmp))
    }
  }else if(is.numeric(row)&(row>nrow(in_tibble))){
    stop('row number requested was higher than the number of rows in input tibble!')
  }else if(!is.numeric(row)){
    stop('row must be a number or a pair name.')
  }
  
  # Check that the prior distribution makes sense:
  if(any(class_prior<0)|any(class_prior>1)){
    stop('Posterior probabilities for degrees of relatedness must be values between 0 and 1.')
  }else if(sum(class_prior)!=1){
    stop('Posterior probabilities for degrees of relatedness must sum to 1.')
  }
  
  # Check that showPlot is a logical variable
  if(!is.logical(showPlot)){
    stop('showPlot must be a logical (TRUE/FALSE) variable.')
  }
  
  # Check that labels is length two is present
  if(!is.null(labels)){
    stopifnot(length(labels) == 2)
  }
  
  # grab PMR info
  x <- in_tibble$mismatch[row]
  N <- in_tibble$nsnps[row]
  M <- in_tibble$ave_rel[row]
  
  # Set title if NULL
  if(is.null(title)){
    title <- in_tibble$pair[row]
  }
  
  # Set max number of points which are plotted at 5000 (uniformly). Set for efficiency.
  if(N>5e3){
    support <- seq(0,N,length.out=5e3)%>%round()
  }else{
    support <-0:N
  }
  
  # Find limit of x-axis.
  XLIM <- (M)+5*sqrt((M)*(1-M)/N)
  
  # Define possible classes.
  kclasses <- c('Same_Twins','First_Degree','Second_Degree','Unrelated')
  
  # find pdf values across "support" for each relatedness category
  distribution.tib <- tibble::tibble(x=rep(support,4),
                                     p.sim=rep((1-0.5^c(1:3,Inf))*M,each=length(support)),
                                     model=rep(kclasses,each=length(support))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(y=(model=='Same_Twins'|model=='Second_Degree')*stats::dbinom(x,N,p.sim)+
                    (model=='First_Degree')*stats::dbinom(x,N,p.sim)+#dbbinom(x,N,alpha,beta)+
                    (model=='Unrelated')*weightedBinom(x,N,p.sim),
                  model=factor(model,levels=kclasses))
  
  # Define plot height for observed PMR below plot
  h <- max(distribution.tib$y)*0.05
  
  # Find highest probability points for dashed lines
  midpoint.tib <- distribution.tib %>%
    dplyr::group_by(model) %>%
    dplyr::filter(y==max(y))
  
  # find critical value
  zstar <- stats::qnorm(0.975)
  
  # make plot of probability curves
  distribution.plot <- distribution.tib %>%
    ggplot2::ggplot(ggplot2::aes(x=x/N,y=y))+
    ggplot2::scale_x_continuous(limits=c(0,XLIM),expand=c(0,0))+
    ggplot2::scale_y_continuous(limits=c(-h,max(distribution.tib$y)),expand=c(0,0))+
    ggplot2::annotate(geom='point',x=x/N,y=-h/2,pch=19,size=2)+
    ggplot2::geom_segment(x=x/N-zstar*sqrt((x/N)*(1-x/N)/N),xend=x/N+zstar*sqrt((x/N)*(1-x/N)/N),y=-h/2,yend=-h/2)+
    ggplot2::geom_segment(x=x/N-zstar*sqrt((x/N)*(1-x/N)/N),xend=x/N-zstar*sqrt((x/N)*(1-x/N)/N),y=-h/4,yend=-3*h/4)+
    ggplot2::geom_segment(x=x/N+zstar*sqrt((x/N)*(1-x/N)/N),xend=x/N+zstar*sqrt((x/N)*(1-x/N)/N),y=-h/4,yend=-3*h/4)+
    ggplot2::theme_bw()+
    ggplot2::xlab(paste0('Pairwise Mismatch Rate (n=',N,')'))+
    ggplot2::ylab(NULL)+
    ggplot2::theme(legend.position='bottom')+
    ggplot2::geom_ribbon(ggplot2::aes(ymin=0,ymax=y,fill=model),col='black',alpha=0.5)+
    ggplot2::geom_segment(data=midpoint.tib,ggplot2::aes(x=p.sim,xend=p.sim,y=0,yend=y*0.99),linetype='dashed',col='black')
  
  # return just distribution plot is asked for
  if(which_plot == 1){
    return(distribution.plot)
  }
  # Calculate and plot posterior probability values
  posterior.tibble <- tibble::tibble(posterior=c(stats::dbinom(x,N,(1-0.5^c(1:3))*M),weightedBinom(x,N,M))*class_prior/
                                       sum(c(stats::dbinom(x,N,(1-0.5^c(1:3))*M),weightedBinom(x,N,M))*class_prior),
                                     model=factor(kclasses,levels=kclasses))
  posterior.tibble <- posterior.tibble %>%
    dplyr::mutate(
      label =label_formatter(posterior)
    )
  # return(posterior.tibble)
  posterior.plot <- posterior.tibble %>%
    ggplot2::ggplot(ggplot2::aes(x=model,y=posterior))+
    ggplot2::theme_bw()+
    ggplot2::geom_bar(stat='identity',ggplot2::aes(fill=model),col='black',alpha=0.8)+
    ggplot2::geom_label(
      ggplot2::aes(
        label=label,
        colour=model
      ), parse = TRUE
    )+
    ggplot2::ylab("Posterior Probability")+
    ggplot2::xlab('Model')+
    ggplot2::scale_x_discrete(guide=ggplot2::guide_axis(n.dodge=2))+
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))
  # Return just posterior plot is asked for
  if(which_plot == 2){
    return(posterior.plot)
  }
  # Combine plots
  # Old code left in as archaeology. Why would you not in a paper for relationship between samples.
  # gg1 <- gridExtra::arrangeGrob(distribution.plot+ggplot2::theme(legend.position='none'),
  #                               posterior.plot+ggplot2::theme(legend.position='none'),
  #                               nrow=1,
  #                               top = grid::textGrob(title,gp=grid::gpar(fontsize=20,font=3)))
  plot_grid <- ggpubr::ggarrange(
    distribution.plot,
    posterior.plot,
    nrow = 1,
    legend = "none",
    labels = labels
  )
  plot_grid <- ggpubr::annotate_figure(
    plot_grid,
    top = ggpubr::text_grob(title,size = 20)
  )
  # Plot if required
  if(showPlot){
    return(plot_grid)
  }
  invisible(plot_grid)
}

utils::globalVariables(
  c("digit", "snp")
)
#' process Eigenstrat data - alternative version
#'
#'A function that takes paths to an eigenstrat trio (ind, snp and geno file) and
#'returns the pairwise mismatch rate for all pairs on a thinned set of SNPs.
#'Options include choosing thinning parameter, subsetting by population names,
#'and filtering out SNPs for which deamination is possible.
#'
#' @param indfile path to eigenstrat ind file
#' @param genofile path to eigenstrat geno file.
#' @param snpfile path to eigenstrat snp file.
#' @param filter_length the minimum distance between sites to be compared (to reduce the effect of LD).
#' @param pop_pattern a character vector of population names to filter the ind file if only some populations are to compared.
#' @param filter_deam a TRUE/FALSE for if C->T and G->A sites should be ignored.
#' @param outfile (OPTIONAL) a path and filename to which we can save the output of the function as a TSV,
#' if NULL, no back up saved. If no outfile, then a tibble is returned.
#' @param chromosomes the chromosome to filter the data on.
#' @param verbose controls printing of messages to console
#'
#' @return out_tibble: A tibble containing four columns:
# - pair [char]: the pair of individuals that are compared.
# - nsnps [numeric]: the number of overlapping snps that were compared for each pair.
# - mismatch [numeric]: the number of sites which did not match for each pair.
# - pmr [numeric]: the pairwise mismatch rate (mismatch/nsnps).
#' @import data.table
#' @export
#'
#' @examples
#' # Use internal files to the package as an example
#' indfile <- system.file("extdata", "example.ind.txt", package = "BREADR")
#' genofile <- system.file("extdata", "example.geno.txt", package = "BREADR")
#' snpfile <- system.file("extdata", "example.snp.txt", package = "BREADR")
#' processEigenstrat(
#' indfile, genofile, snpfile,
#' filter_length=1e5,
#' pop_pattern=NULL,
#' filter_deam=FALSE
#' )
processEigenstrat <- function(indfile, genofile, snpfile,
                              filter_length=NULL, pop_pattern=NULL, filter_deam=FALSE,
                              outfile=NULL, chromosomes=NULL, verbose=TRUE, byPop=FALSE, minMerge=2){
  
  # Check to see if eigenstrat files exist
  if(!file.exists(indfile)){
    stop(paste0('indfile ',indfile,' does not exist!'))
  }
  if(!file.exists(genofile)){
    stop(paste0('genofile ',genofile,' does not exist!'))
  }
  if(!file.exists(snpfile)){
    stop(paste0('snpfile ',snpfile,' does not exist!'))
  }
  
  # minMerge should be a number ≥2
  if((!is.numeric(minMerge))|(minMerge<2)){
    stop('minMerge must be a number ≥ 2.')
  }
  
  # Just check to see if a length filter was included and is reasonable, and if not,
  # let the user know of the default and/or problem.
  if(is.null(filter_length)){
    filter_length <- 1e5
    if(verbose){
      cat(sprintf('No site distance filter.\nUsing default minimum of 1e5.\n'))
    }
  }
  if(!is.numeric(filter_length)|(filter_length<=0)){
    stop('filter_length must be a positive number.')
  }
  
  # Check if all pop_names are in the possible populations.
  # Change to delim so works in different countries.
  ind_raw <- read_ind(indfile)
  if(!all(pop_pattern%in%ind_raw$pop)){
    paste0('Following populations in pop_pattern do not exist in indfile: ',
           paste0(setdiff(pop_pattern,ind_raw$pop),collapse=', ')) %>%
      stop()
  }
  
  # Check that filter deam is a logical variable.
  if(!is.logical(filter_deam)){
    stop('filter_deam must be a logical (boolean) variable, i.e., TRUE or FALSE.')
  }
  
  # Check that the folder in which we plan to save the outfile actually exists.
  if(!is.null(outfile)){
    if(!file.exists(dirname(outfile))){
      stop(paste0('Output folder ',dirname(outfile),' does not exist!'))
    }
  }
  
  # Collect SNP details
  if(verbose){
    cat('Reading in SNP data.\n')
  }
  ## Read in SNP file
  snps <- read_snp(snpfile)
  ## Add relative position
  snps <-
    snps %>%
    dplyr::mutate(
      relative_pos=1:dplyr::n()
    )
  ## Mutate and filter for DEAM
  if(filter_deam){
    snps <-
      snps %>%
      dplyr::mutate(
        deam=dplyr::case_when(
          anc=='C'&der=='T' ~ T,
          anc=='G'&der=='A' ~ T,
          T ~ F)
      ) %>%
      dplyr::filter(!deam) %>%
      dplyr::select(-deam)
  }
  
  # Filter for requested chromosomes
  available_chr <- snps$chr %>%
    unique()
  if(!is.null(chromosomes)){
    if(!any(chromosomes%in%available_chr)){
      stop(paste0('Chromosome names do not match SNP file:\n',
                  'Requested: ',paste0(chromosomes,collapse=', '),'\n',
                  'Available: ',paste0(available_chr,collapse=', '),'\n'))
    }else{
      snps <- snps %>%
        dplyr::filter(chr%in%chromosomes)
    }
  }
  
  # Check that byPop is logical
  if(!is.logical(byPop)){
    stop('byPop must be logical (TRUE/FALSE)')
  }
  
  # Calculate groupings
  
  
  # Print information about chromosomes
  if(verbose){
    cat(paste0('Analysing chromosomes:\n',paste0(unique(snps$chr),collapse=', '),'\n'))
  }
  # Collect (potentially filtered) ind details
  if(is.null(pop_pattern)){
    ind <-  ind_raw %>%
      dplyr::mutate(row_number=1:dplyr::n())
  }else{
    ind <- ind_raw %>%
      dplyr::mutate(row_number=1:dplyr::n()) %>%
      dplyr::filter(pop%in%pop_pattern)
  }
  
  # Throw error if ind has no rows.
  n_ind <- nrow(ind)
  if(n_ind<=1){
    if(is.null(pop_pattern)){
      stop(paste0(n_ind,' individuals in final, filtered ind file. Check ind file.'))
    }else{
      stop(paste0(n_ind,' individuals in final, filtered ind file. Check ind file, and check pop_pattern.'))
    }
  }
  geno_list <- vector(mode='list',length=n_ind)
  
  if(verbose){
    cat('Starting to read in genotype data.\n')
  }
  pb1 <- utils::txtProgressBar(min=1,max=length(geno_list),initial=1,style=3)
  dt <- data.table::fread(genofile, col.names = "snp", colClasses = "character")
  for(i in 1:n_ind){
    dt[, digit := stringr::str_sub(snp, ind$row_number[i], ind$row_number[i]),]
    geno_list[[i]] <- (dt[[2]] %>% dplyr::na_if("9"))[snps$relative_pos]
    utils::setTxtProgressBar(pb1,i)
    
  }
  if(verbose){
    cat('  Complete.\n\n')
    
    cat('Starting to compare genotypes and calculate PMR.\n')
  }
  
  if(!byPop){
    pop_group <- rep('Merged',nrow(ind))
  }else{
    pop_group <- ind %>%
      dplyr::group_by(pop) %>%
      dplyr::mutate(N=n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(pop_group=ifelse(N<minMerge,
                                     'Merged',
                                     pop)) %>%
      dplyr::pull(pop_group)
  }
  
  # Check if things work
  if(min(table(pop_group))==1){
    stop('Even after possibly merging, smallest population group size is 1.')
  }else{
    POPS <- pop_group %>%
      base::unique() %>%
      base::sort()
    if(byPop){
      if(sum(pop_group=='Merged')>0){
        ind %>%
          dplyr::slice(which(pop_group=='Merged')) %>%
          dplyr::pull(pop) %>%
          base::unique() %>%
          base::sort() %>%
          base::paste0(collapse='/') %>%
          base::sprintf('Merged group made from: %s\n',.) %>%
          message()
      }
    }
  }
  
  for(k in 1:length(POPS)){
    if(k==1){
      comp_list <- base::vector(mode='list',
                                length=length(POPS))
    }
    
    currInds <- ind %>%
      dplyr::filter(pop_group==POPS[k])
    n_ind <- nrow(currInds)
    if(byPop){
      sprintf('Comparing within population: %s (%i comparisons)',POPS[k],choose(n_ind,2)) %>%
        message()
    }
    counter <- 1
    out_tibble <- tibble::tibble(pair=rep('x',choose(n_ind,2)),nsnps=0,mismatch=0,pmr=0)
    pb2 <- utils::txtProgressBar(min=0,max=choose(n_ind,2),initial=1,style=3)
    for(i in 1:(n_ind-1)){
      ind_i <- currInds$ind[i]
      vi <- geno_list[[i]]
      for(j in (i+1):n_ind){
        ind_j <- currInds$ind[j]
        vj <- geno_list[[j]]
        
        snps_filtered_prefiltered <- snps  %>%
          dplyr::mutate(row_num=1:dplyr::n()) %>%
          dplyr::filter(!is.na(vi),!is.na(vj))
        
        if(nrow(snps_filtered_prefiltered)>0){
          snps_filtered <- snps_filtered_prefiltered %>%
            dplyr::group_by(chr) %>%
            dplyr::mutate(keep=getfilter(site,filter_length)) %>%
            dplyr::ungroup() %>%
            dplyr::filter(keep) %>%
            dplyr::select(-keep)
          
          out_tibble$pair[counter] <- paste0(ind_i,' - ',ind_j)
          out_tibble$nsnps[counter] <- (!is.na(vi[snps_filtered$row_num])&!is.na(vj[snps_filtered$row_num])) %>% sum()
          out_tibble$mismatch[counter] <- (vi[snps_filtered$row_num]!=vj[snps_filtered$row_num]) %>% sum(na.rm=T)
          out_tibble$pmr[counter] <-  out_tibble$mismatch[counter]/out_tibble$nsnps[counter]
        }else{
          out_tibble$pair[counter] <- paste0(ind_i,' - ',ind_j)
          out_tibble$nsnps[counter] <- 0
          out_tibble$mismatch[counter] <- 0
          out_tibble$pmr[counter] <-  NA
        }
        counter <- counter+1
        utils::setTxtProgressBar(pb2,counter)
      }
    }
    comp_list[[k]] <- out_tibble
  }
  
  OUT_TIBBLE <- comp_list %>%
    dplyr::bind_rows()
  
  if(!is.null(outfile)){
    if(verbose){
      cat(sprintf('\n\nSaving processed data to %s.\n',outfile))
    }
    readr::write_delim(OUT_TIBBLE,file=outfile,delim='\t')
  }
  if(verbose){
    cat('\nComplete.\n\n')
  }
  
  OUT_TIBBLE %>% 
    return()
}

plotROLL <- function(in_tibble,indfile,nsnps_cutoff=NULL,fntsize=7){
  # Check in_tibble is real
  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('pair','nsnps','mismatch','pmr')
  if(!all(names(in_tibble)%in%col_names_check)){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }
  
  # Test that the overlapping snp cut off is sensible.
  if(is.null(nsnps_cutoff)){
    nsnps_cutoff <- 5e2
    if(verbose){
      cat(sprintf('No minimum number of overlapping SNPs given.\nUsing default minimum of 500.\n'))
    }
  }
  if((nsnps_cutoff<=0)){
    stop('The cut off for the number of overlapping snps must be greater than 0.')
  }else if(nsnps_cutoff>=max(in_tibble$nsnps,na.rm=T)){
    stop('The cut off for the number of overlapping snps cannot be greater than the maximum number of overlapping SNPs that were observed!')
  }
  # Check indfile exists
  if(!file.exists(indfile)){
    stop(paste0('indfile ',indfile,' does not exist!'))
  }
  
  IND <- read_tsv(indfile,
                  show_col_types=FALSE,
                  col_names=c('ind','sex','pop')) %>%
    dplyr::select(-sex)
  
  plot_data <- in_tibble %>%
    tidyr::separate(col=pair,
                    into=c('ind1','ind2'),
                    sep=' - ',
                    remove=FALSE) %>%
    dplyr::left_join(IND,
                     by=c('ind1'='ind')) %>%
    dplyr::left_join(IND,
                     by=c('ind2'='ind'),
                     suffix=c('1','2')) %>%
    dplyr::mutate(group=ifelse(pop1==pop2,
                               pop1,
                               'Between Groups'))
  
  # Check that all individuals are accounted for
  fileInds <- c(plot_data$ind1,
                plot_data$ind2) %>%
    base::unique() %>%
    base::sort()
  
  missingInds <- base::setdiff(x=fileInds,
                               y=IND$ind)
  
  if(length(missingInds)>0){
    missing_str <- missingInds %>%
      paste0(collapse=', ') %>%
      sprintf('The following individuals are in the BREADR analysis, but not the ind file: %s',.)
    stop(missing_str)
  }
  
  LEVELS <- plot_data$group %>%
    base::unique() %>%
    base::sort() %>%
    base::setdiff('Between Groups') %>%
    c(.,'Between Groups')
  
  M <- median(plot_data$pmr)
  
  M_pmr <- plot_data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(unrelated=median(pmr)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(second_degree=((1-0.5^c(3)))*unrelated,
                  first_degree=((1-0.5^c(2)))*unrelated,
                  same_twins=((1-0.5^c(1)))*unrelated) %>%
    tidyr::gather(relatonship,
                  pmr,2:5)
  
  outplot <- plot_data %>%
    dplyr::mutate(group=group %>%
                    factor(levels=LEVELS)) %>%
    ggplot2::ggplot(aes(x=group,
               y=pmr,
               fill=group))+
    ggplot2::theme_bw()+
    ggplot2::stat_boxplot(geom='errorbar',
                 width=0.4)+
    ggplot2::geom_boxplot(outlier.shape=21)+
    ggplot2::geom_hline(yintercept=M,
                        linetype='dotted',
                        col='black',
                        linewidth=1)+
    geom_hline(data=M_pmr,
               aes(yintercept=pmr,
                   col=group),
               linetype='dashed',
               linewidth=1)+
    ggplot2::facet_grid(.~group,
                        scales='free_x')+
    ggplot2::ylab("PMR")+
    ggplot2::xlab(" ")+
    harrypotter::scale_fill_hp_d(option='Ravenclaw',drop=FALSE,name=NULL)+
    harrypotter::scale_colour_hp_d(option='Ravenclaw',drop=FALSE,name=NULL,guide='none')+
    ggplot2::scale_shape_manual(drop=FALSE,name=NULL,values=c(21:24))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1, size=fntsize),
                   strip.text.x.top = ggplot2::element_blank(),
                   legend.position='bottom')+
    ggplot2::guides(fill=ggplot2::guide_legend(override.aes=list(size=5)))
  
  return(outplot)
}

#' saveSLICES
#'
#' Plots all pairwise diagnostic plots (in a tibble as output by callRelatedness),
#' as produced by plotSLICE, to a folder.
#' Options include the width and height of the output files, and the units in
#' which these dimensions are measured.
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param outFolder the folder into which all diagnostic plots will be saved
#' @param width the width of the output PDFs.
#' @param height the height of the output PDFs.
#' @param units the units for the height and width of the output PDFs.
#' @param verbose Controls the printing of progress to console.
#'
#' @return nothing
#' @export
#'
#' @examples
#' \donttest{
#' saveSLICES(relatedness_example[1:3, ], outFolder = tempdir())
#' }
saveSLICES <- function(in_tibble,outFolder=NULL,width=297,height=210,units='mm', verbose = TRUE){
  
  
  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if((!all(names(in_tibble)==col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }
  
  # Test to see if outfolder is set and exists
  if(is.null(outFolder)){
    stop('No output folder given.')
  }else{
    if(!dir.exists(outFolder)){
      stop('Output folder does not exist.')
    }
  }
  
  # Make sure output folder string ends in forward slash
  if(stringr::str_sub(outFolder,-1)!='/'){
    outFolder <- paste0(outFolder,'/')
  }
  
  # Test if plotting dimensions are numeric, and if the units are valid.
  if(any(!sapply(c(width,height),is.numeric))){
    stop('Height and width must be numeric values')
  }
  if(!units%in%c(c("in","cm","mm","px"))){
    stop("Plot dimension units must be 'in', 'cm', 'mm' or 'px'")
  }
  
  # Start plotting!
  if(verbose){
    cat('Starting pairwise plot creation.\n')
  }
  
  
  pb = utils::txtProgressBar(min=1,max=nrow(in_tibble),initial=1,style=3)
  for(i in 1:nrow(in_tibble)){
    v.sm <- plotSLICE(in_tibble,i,showPlot=F)
    ggplot2::ggsave(v.sm,
                    filename=paste0(outFolder,gsub(' - ','_',in_tibble$pair[i]),'.pdf'),
                    width=width,height=height,units=units)
    utils::setTxtProgressBar(pb,i)
  }
  if(verbose){
    cat(paste0('\nCompleted\nAll plots in: ',outFolder,'\n'))
  }
  return(TRUE)
}

#' test_degree
#'
#' Test if a degree of relatedness is consistent with an observed PMR
#'
#' @param in_tibble a tibble that is the output of the callRelatedness() function.
#' @param row either the row number or pair name for which the posterior distribution is to be plotted.
#' @param degree the degree of relatedness to be tested.
#' @param verbose a logical (boolean) for whether all test output should be printed to screen.
#'
#' @return the associated p-value for the test
#' @export
#'
#' @examples
#' test_degree(relatedness_example, 1, 1)
test_degree <- function(in_tibble,row,degree,verbose=TRUE){
  # Test that the in_tibble is of the correct form
  if(nrow(in_tibble)==0){
    stop('The input tibble is empty')
  }
  col_names_check <- c('row','pair','relationship','pmr','sd','mismatch','nsnps','ave_rel','Same_Twins','First_Degree','Second_Degree','Unrelated')
  if((!all(names(in_tibble)==col_names_check))){
    stop(paste0('Input tibble/data.frame must have 4 columns named: ',paste0(col_names_check,collapse=', '),'.'))
  }
  
  # check row is a numeric or char variable, and if char, find the row.
  if(is.character(row)){
    rowtmp <- row
    row <- which(row==in_tibble$pair)
    if(length(row)==0){
      stop(paste0('There is no pair called: ',rowtmp))
    }
  }else if(is.numeric(row)&(row>nrow(in_tibble))){
    stop('row number requested was higher than the number of rows in input tibble!')
  }else if(!is.numeric(row)){
    stop('row must be a number or a pair name.')
  }
  
  if(!is.numeric(degree)){
    stop('Degree must be numeric.')
  }
  # if(degree>10|degree<0|(degree%%1!=0)){ # Removed constraint that integer
  if(degree>10|degree<0){
    stop('Degree must be a value k such that 0<k<10.')
  }
  
  M <- in_tibble$ave_rel[row]
  pk <- (1-0.5^(degree+1))*M
  N <- in_tibble$nsnps[row]
  x <- in_tibble$mismatch[row]
  pair <- in_tibble$pair[row]
  # pval <- 2*stats::pbinom(x,N,pk,ifelse(x>pk*N,F,T))
  pval <- min(
    c(
      2*stats::pbinom(x,N,pk,T),
      2*stats::pbinom(x,N,pk,F)
    )
  )
  # CI <- x/N+c(-1,1)*1.96*sqrt(x/N*(1-x/N)/N)
  pvaltext <- ifelse(pval<0.0001,
                     sprintf('p-value          : %s\n',format(pval,scientific=T)),
                     sprintf('p-value          : %.4f\n',pval))
  
  
  degree_string <- dplyr::case_when(degree==1 ~ "1st",
                                    degree==2 ~ "2nd",
                                    degree==3 ~ "3rd",
                                    T ~ paste0(degree,"th"))
  if(verbose){
    cat(sprintf('Testing H0       : "%s" are %s-degree relatives.\n',pair,degree_string))
    cat(sprintf('Expected PMR     : %.4f\n',pk))
    cat(sprintf('Observed PMR     : %.4f\n',x/N))
    cat(sprintf('Estimated degree : %.4f\n',convertP(x/N,M)))
    cat(pvaltext)
    cat(ifelse(pval>=0.05,'Decision         : Retain H0\n','Decision         : Reject H0\n'))
  }
  return(pval)
}

