#main.R

library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(knitr)
library(resample)
library(cluster)
library(ggdendro)
library(ggplot2)

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return dataframe
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data.csv')
#' 
read_expression_table <- function(filename) {
  
  expr_mat <- tibble::as_tibble(t(read.table(
    filename,
    header = TRUE,
    sep = ' '
  )),
  rownames = "sample_ids")
  return(expr_mat)
}

#' First step of filtering noise.
#' Keep only genes expressed in 20% of samples, by saving only genes that have
#' greater than log2(15) gene expression in at least 20% of samples.
#' 
#' @param dataframe unfiltered for noise
#' 
#' @return dataframe filtered for noise
#' 

first_filter <- function(data) {
  #drop all cols that dont meet criteria
  drop <- which(sapply(data, function(x)all(sum(x > log2(15))/35 < 0.2)))
  #combine sample_ids with filtered data
  filtered <- data[,-drop]
  return(filtered)
}

#' Second step of filtering.
#' Have a variance significantly different from the median variance of all 
#' probe sets with p<0.01.
#' 
#' chi_square function calculates t-stat for a vecotr of data then checks if 
#' t-stat falls in the upper and lower range.
#' 
#' sec_filter calculates upper and lower range given degrees of freedom
#' and significance value
#' 
#' Test statistic = (N-1)(s/sigma)^2
#' 
#' @param dataframe
#' 
#' @return dataframe filtered for noise
#' 

chi_square <- function(data, upper, lower, sig, df = 34){
  v = var(data)
  test = df*v/(sig^2)

  if(test > lower || test < upper){
    return_val = FALSE
  } else {
    return_val = TRUE
  }
  return(return_val)
}

sec_filter <- function(data, alpha = 0.01, df = 34){
  
  sample_ids <- data[,1]
  data <- data[,-1]
  
  up <- qchisq(alpha/2,df)
  low <- qchisq(1-alpha/2,df)
  
  sds <- data[-1,] %>%
    lapply(sd)
  sds <- unname(unlist(sds))

  sigma <- median(sds, na.rm = T)
  
  
  drop <- which(sapply(data, function(x)all({chi_square(x, up, low, sigma) == TRUE})))
  data <- data[,-drop]
  
  data <- cbind(sample_ids, data)
  
  return(data)
  
}

#'
#'third filter calculates coefficient of variation for each column in 
#'a datafrmae and checks to see if the coefficient reaches a predetermined 
#'threshold.
#'
#'@param dataframe
#'
#'@return filtered dataframe
#'

third_filter <- function(data, thresh = 0.186){
  sample_ids <- data[,1]
  data <- data[,-1]
  
  drop <- which(sapply(data, function(x){sd(x)/mean(x) < 0.186}))
  
  data <- data[,-drop]
  data <- cbind(sample_ids, data)
  
  return(data)
}

#' clustering
#' 
#' use the cluster package and the associated distance and hclust functions
#' to perform clustering
#' 
#' then generate plot using ggplot2
#' 
#' @param df
#' 
#' @return ggplot object
#' 
best_method <- function(data){
  m <- c( "average", "single", "complete", "ward")
  names(m) <- c( "average", "single", "complete", "ward")
  
  # function to compute coefficient
  ac <- function(x) {
    agnes(data, method = x)$ac
  }
  
  lapply(m, function(x){ac(x)})
  

}

clustering <- function(data){
  
  rownames(data) <- data[,1]
  data <- data[,-1]
  
  model <- agnes(data, method = 'ward')
  
  dhc <- as.dendrogram(model)
  # Rectangular lines
  ddata <- dendro_data(dhc, type = "rectangle")
  plot <- ggplot(segment(ddata)) + 
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(ddata), 
              aes(x = x, y = y, label = label), size = 3, angle = 45) +
    theme_dendro()
  plot

}
