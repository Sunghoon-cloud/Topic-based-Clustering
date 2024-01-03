This repo provides the source codes and the three datasets to help any model users to replicate the results in Kim et al. (2024) "A Topic-based Segmentation Model for Identifying Segment-Level Drivers of Star Ratings from Unstructured Text Reviews", the Journal of Marketing Research  

A. Codes and files related to Datasets.

(A-1) Generating an illustrative simulation study dataset 

 SimGen_Topic-based_Segmentation_Model.R
: This R code generates an illustrative simulation study dataset described in Kim et al. (2024).

(A-2) Processing Data from the original Yelp Data

Process_Yelp_Text.R
: This code helps to extract the dataset used in the Kim et al. (2024). Once any user downloads the Yelp datasets from the provided link, the user can prepare the same datasets used in the Kim et al. (2024) by running the provided R code.

   (A-2-a) This code extracts the dataset (e.g., Arizona-based restaurant brands with more than 10 reviews) from the original Yelp data for restaurant-level segmentation analysis in the first empirical application study section.  

   (A-2-b) This code extract the dataset (e.g., An Italian restaurant with the largest number of reviews) from the original Yelp data for customer-level segmentation study in the second empirical applcation study. If model user wants to see the results relatively quickly, the dataset is recommended given the size is smaller and running time is faster. 

B. Models

(B-1) Preprocessing Text.R
: This is the R code to extract word-frequency matrix (X) from unstructured text reviews. The step should be employed to the outputs of (A-2) step before running the Topic based Segmentation Model to preprocess the unstructured review texts. Any users need to run these codes before running the main model code of (B-2) for all empirical studies (1. restaurant-level segmentation, 2. customer-level segmentation, 3. professor segmentation). 

(B-2) Topic-based Segmentation Model.R
: This is the main R code of the Topic based Segmentation Model estimates the sequential choice model by calling "Est_SCM_202105.R".

(B-3) Identifying_Dimensions.R

: This is the R code of post-hoc analyses to identify dimensions by generating dendrograms (akin to heatmap) to understand co-occurrence structure between words.
