This repo provides the R package with sample data and instruction script to help any potential users to use the methods introduced in Kim et al. (2024) "A Topic-based Segmentation Model for Identifying Segment-Level Drivers of Star Ratings from Unstructured Text Reviews", the Journal of Marketing Research. Also, this repo includes source codes and datasets for replicating the results in ther paper. This repo includes 5 parts as below.  

<b> A. Replicating Illustrative Simulation study </b> (this folder includes one file)

 (file)<i> 1. Full codes for replicating Illustrative simulation study.R</i>

: By running this source codes from begining to the bottom in R, you can replicate the illustrative simulation study in the paper. 

This R code file includes (1) codes for generating for a simulation data, (2) codes for estimating the topic-based segmentation model, (3) codes for 	identifying topics (frequently co-occurring neighbor words), and (4) codes for generating hierarchical dendrograms as shown in figure 2 in the paper.
 
<b> B. Model Packages</b>  (this folder includes four files)

We developed user-friendly packages with a simple sample simulation data which can help any model users simply apply the package to their datasets. In this folder we included following four files. 

(file)<i>  2. Instruction script to run Topic-based Segmentation Packages.txt</i>
 	
: This is Instruction for running the two packages of "Package_Topic-based_Segmentation.R" and "Package_Dendrogram_Dimensions.R" with a simple sample data for helping any users to apply the model to their datasets.
	
  (file) <i> Package_Topic-based_Segmentation.R</i>

: This is a function package for topic-based segmentation (i.e., latent class regression with group variable selection).

 (file)  <i> Package_Dendrogram_Dimensions.R</i>

: This is a package for generating dendrogram for identifying topics.

  (file) <i> sample_data.csv</i>

: This is a sample data. True parameter values are presented in the "2. Instruction script to run Topic-based Segmentation Packages.txt".

<b> C. Empirical Studies of Customer-level Segmentation with YELP </b> (this folder includes 2 files)

Note, given the terms of agreement by YELP, we provide the link to download the same Yelp datasets (https://www.kaggle.com/datasets/yelp-dataset/yelp-dataset/versions/6) used in Kim et al. (2024). The provided codes below help replicate the customer-level segmentation results with the raw Yelp data downloaded from the link. 

 (file)  <i> 3-a. Preprocessing raw Yelp Reviews for Customer-level Segmentation.txt</i>

: This is R code for preprocessing the downloaded Yelp data and building DV and IVs matrix for customer-level segmentation study. 
	
   (file)<i> 3-b. Instruction for replicating Customer-level Segmentation analysis.txt</i>

: This is instruction script for replicating customer-level segmentation study with Yelp

<b> D. Empirical Studies of Restaurant-level Segmentation with YELP </b> (this folder includes 2 files)

The provided codes below help replicate the restaurant-level segmentation results with the raw Yelp data downloaded from the link.

 (file)  <i> 4-a. Preprocessing raw Yelp reviews_Restaruant Segmentation.txt</i>

: This is R code for preprocessing the downloaded Yelp data and building DV and IVs matrix for restaurant-level segmentation study.

 (file)  <i> 4-b. Instruction for replicating restaurant-level segmentation analysis.txt</i>

: This is instruction script for replicating restaurant-level segmentation study with Yelp

<b> E. Empirical application study with ProfessorRatings</b>  (this folder includes 2 files)

(file)   <i>5. Instructions for replicating Professor ratings review study.txt</i>

: This is instruction script to replicate the Professor ratings reviews study in Web Appendix.

(file)   <i> Professor_ratings_words_freq.csv</i>

: Data file of Ratings and word frequency matrix extracted from the Professor ratings reviews which were scrapped by author.
