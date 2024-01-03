This repo provides the source codes and the three datasets to help any model users to replicate the results in Kim et al. (2024) "A Topic-based Segmentation Model for Identifying Segment-Level Drivers of Star Ratings from Unstructured Text Reviews", the Journal of Marketing Research  

A. Methods

(A-1) Preprocessing Text.R
: This is the R code to extract word-frequency matrix (X) from unstructured text reviews. The step should be employed to the outputs of (A-2) step before running the Topic based Segmentation Model to preprocess the unstructured review texts. Any users need to run these codes before running the main model code of (B-2) for all empirical studies (1. restaurant-level segmentation, 2. customer-level segmentation, 3. professor segmentation). 

(A-2) Topic-based Segmentation Model.R
This is the main R code of the Topic based Segmentation Model estimates the sequential choice model by calling "Est_SCM_202105.R".

(A-3) Identifying_Dimensions.R
This is the R code of post-hoc analyses to identify dimensions by generating dendrograms (akin to heatmap) to understand co-occurrence structure between words.

B. Simulation study

(B-1) Instruction for replicating simulation study in Kim et al. (2024)

(B-2) Generating an illustrative simulation study dataset 

 SimGen_Topic-based_Segmentation_Model.R
:This R code generates an illustrative simulation study dataset described in Kim et al. (2024).

(B-3) the Generated simulation dataset
Users can download this simulation data and follow the (B-1) instruction to replicate the result in Kim et al. (2024).

C. Empirical application study with YELP
Yelp Dataset Terms of Use, we can't distribute the Yelp data. So, we share the link (https://www.kaggle.com/datasets/yelp-dataset/yelp-dataset/versions/6) where anyone can download the same data. Then, this code enables any users to extract the dataset used in the Kim et al. (2024). Once any user downloads the Yelp datasets from the provided link, the user can prepare the same datasets used in the Kim et al. (2024) by running the provided R code.

   (C-1) Instruction for replicating Yelp Empirical Study in Kim et al. (2024) 

   (C-2) Process_Arizona_Restaurant_Brands.R 
This code extracts the dataset (e.g., Arizona-based restaurant brands with more than 10 reviews) from the original Yelp data for restaurant-level segmentation analysis in the first empirical application study section.  

   (C-1) Process_One_Restaurant_Largest_Reviews.R 
This code extract the dataset (e.g., An Italian restaurant with the largest number of reviews) from the original Yelp data for customer-level segmentation study in the second empirical application study. If model user wants to see the results relatively quickly, the dataset is recommended given the size is smaller and running time is faster. 

D. Empirical application study with ProfessorRatings

   (D-1) Instruction for replicating ProfessorRating Empirical Study in Web Appendix of Kim et al. (2024)

   (D-2) Process_ProfessorRatings.R  
This code helps to prepare ProfessorRating Empirical Study in Web Appendix of Kim et al. (2024)

   (D-3) Professor_Rating.rds
This is the R data for Professor Ratings Review.

