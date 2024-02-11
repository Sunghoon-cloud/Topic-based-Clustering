This repo provides the R package with sample data and instruction script to help any potential users to use the methods introduced in Kim et al. (2024) "A Topic-based Segmentation Model for Identifying Segment-Level Drivers of Star Ratings from Unstructured Text Reviews", the Journal of Marketing Research. Also, this repo includes source codes and datasets for replicating the results in ther paper. This repo includes 5 parts as below.  

<b> A. Replicating Illustrative Simulation study </b> (this folder includes one file)

 (file)<i> 1. Full codes for replicating Illustrative simulation study.R</i>

: By running this source codes from begining to the bottom in R, you can replicate the illustrative simulation study in the paper. [see Table 2 and Figure 2 in main text of the paper]

This R code file includes (1) codes for generating for a simulation data, (2) codes for estimating the topic-based segmentation model, (3) codes for 	identifying topics (frequently co-occurring neighbor words), and (4) codes for generating hierarchical dendrograms as shown in figure 2 in the paper.
 
<b> B. Model Packages</b>  (this folder includes four files)

We developed user-friendly packages with a simple sample simulation data which can help any model users simply apply the package to their datasets. In this folder we included following four files. 

(file)<i>  2. Instruction script to run Topic-based Segmentation Packages.txt</i>
 	
: This is Instruction for running the two packages of "Package_Topic-based_Segmentation.R" and "Package_Dendrogram_Dimensions.R" with a simple sample data for helping any users to apply the model to their datasets.
	
  (file) <i> Package_MixtureRegression_GroupVariableSelection.R</i> 

: This is a function package for topic-based segmentation (i.e., latent class regression with group variable selection). [see Tables 5, 6, 7, and 10 in main text; Tables E-4, E-5, F-1, F-2, F-3, G-1, G-2, G-4 and G-5 in Web Appendices]

 (file)  <i> Dendrogram.R</i>

: This is a function code for generating dendrogram for identifying topics.[see Figure 3 in main text; Figures F-1 and G-1 in Web Appendices]

  (file) <i> sample_data.csv</i>

: This is a sample data. True parameter values are presented in the "2. Instruction script to run Topic-based Segmentation Packages.txt".

(file) <i> web appendix M-H source code for coestimating tuning parameter.txt </i>

: This is a source code for Metropolis Hastings algorithm for coestimating lambda (see Web Appendix B of the referred paper).

<b> C. Empirical Studies of Customer-level Segmentation with YELP </b> (this folder includes 2 files)

Note, given the dataset terms of use by YELP, we provide the link to download the same Yelp datasets (https://www.kaggle.com/datasets/yelp-dataset/yelp-dataset/versions/6) used in Kim et al. (2024). The provided codes below help replicate the customer-level segmentation results with the raw Yelp data downloaded from the link. 

 (file)  <i> 3-a. Preprocessing raw Yelp Reviews for Customer-level Segmentation.txt</i>

: This is R code for preprocessing the downloaded Yelp data and building DV and IVs matrix for customer-level segmentation study. 
	
   (file)<i> 3-b. Instruction for replicating Customer-level Segmentation analysis.txt</i>

: This is instruction script for replicating customer-level segmentation study by using the provided packages. [see Table 10 in main text; Tables F-1, F-2, and F-3 and Figure F-1 in the Web Appendix]


<b> D. Empirical Studies of Restaurant-level Segmentation with YELP </b> (this folder includes 2 files)

The provided codes below help replicate the restaurant-level segmentation results with the raw Yelp data downloaded from the link.

 (file)  <i> 4-a. Preprocessing raw Yelp reviews_Restaruant Segmentation.txt</i>

: This is R code for preprocessing the downloaded Yelp data and building DV and IVs matrix for restaurant-level segmentation study.

 (file)  <i> 4-b. Instruction for replicating restaurant-level segmentation analysis.txt</i>

: This is instruction script for replicating restaurant-level segmentation study by using the provided packages. [see Tables 5, 6 and 7 in main text; Tables E-4 and E-5 and Figure H-1 in Web Appendices]


<b> E. Empirical application study with ProfessorRatings</b>  (this folder includes 2 files)

(file)   <i>5. Instructions for replicating Professor ratings review study.txt</i>

: This is instruction script to replicate the Professor ratings reviews study by using the provided packages in Web Appendix. [see Tables G-1, G-2, G-4 and G-5, and Figures G-1 and H-2 in Web Appendices]


(file)   <i> Professor_ratings_words_freq.csv</i>

: This is data file of Ratings and word frequency matrix extracted from the Professor ratings reviews, which were scrapped by the author from ratemyprofessors.com.


<b>[Disclaimer]</b> 

We developed these R packages and code to be as user-friendly as possible. However, a minimum familiarity with the R language and knowledge of statistics are assumed to properly import the data, run the model, and interpret the results. There is no in-built protection against misuse.

<b>[Technical Help or Problem Report]</b> 

We are expecting the R codes to be easy to use and to be working well. Please send questions or report bugs to Sunghoon Kim (E-mail: Sunghoon.Kim@rutgers.edu).

