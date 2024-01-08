This repo provides the source codes and the three datasets to help any model users to replicate the results in Kim et al. (2024) "A Topic-based Segmentation Model for Identifying Segment-Level Drivers of Star Ratings from Unstructured Text Reviews", the Journal of Marketing Research  

A. Replicating Illustrative Simulation study in Kim et al. (2024)

  Illustrative_Simulation_Study.R

 From the top to bottom, this R code includes: 
  (1) generating a simulation data described in Kim et al. (2024) 
  (2) estimating the topic-based segmentation model
  (3) identifying topics (frequently co-occurring neighbor words)
  (4) generating hierarchical dendrograms as shown in figure 2 in the paper
 
 At the end, we can check identified dimensions compared with the true dimensions.  


B. Packages 

We developed user-friendly packages with a simple sample simulation data for helping easy application.

(B-1) Preprocessing Text.R
: This is the R code to extract word-frequency matrix (X) from unstructured text reviews. The step should be employed to the outputs of (A-2) step before running the Topic based Segmentation Model to preprocess the unstructured review texts. Any users need to run these codes before running the main model code of (B-2) for all empirical studies (1. restaurant-level segmentation, 2. customer-level segmentation, 3. professor segmentation). 

(B-2) Topic-based Segmentation Model.R
This is the main R code of the Topic based Segmentation Model estimates the sequential choice model by calling "Est_SCM_202105.R".

(B-3) Identifying_Dimensions.R
This is the R code of post-hoc analyses to identify dimensions by generating dendrograms (akin to heatmap) to understand co-occurrence structure between words.


C. Empirical application studies with YELP

Yelp Dataset Terms of Use, we can't distribute the Yelp data. So, we share the link (https://www.kaggle.com/datasets/yelp-dataset/yelp-dataset/versions/6) where anyone can download the same data. Then, this code enables any users to extract the dataset used in the Kim et al. (2024). Once any user downloads the Yelp datasets from the provided link, the user can prepare the same datasets used in the Kim et al. (2024) by running the provided R code.

   (C-1) Instruction for replicating Yelp Empirical Study in Kim et al. (2024) 

   (C-2) Process_Arizona_Restaurant_Brands.R 
This code extracts the dataset (e.g., Arizona-based restaurant brands with more than 10 reviews) from the original Yelp data for restaurant-level segmentation analysis in the first empirical application study section.  

   (C-3) Process_One_Restaurant_Largest_Reviews.R 
This code extract the dataset (e.g., An Italian restaurant with the largest number of reviews) from the original Yelp data for customer-level segmentation study in the second empirical application study. If model user wants to see the results relatively quickly, the dataset is recommended given the size is smaller and running time is faster. 

D. Empirical application study with ProfessorRatings

   (D-1) Instruction for replicating ProfessorRating Empirical Study in Web Appendix of Kim et al. (2024)

   (D-2) Process_ProfessorRatings.R  
This code helps to prepare ProfessorRating Empirical Study in Web Appendix of Kim et al. (2024)

   (D-3) Professor_Rating.rds
This is the R data for Professor Ratings Review.
