## read data
load("geno_data.RData")
# 100 SNP panel
reference_breed<-geno_data[1:50,]
test_breed_1<-geno_data[51:100,]
test_breed_2<-geno_data[101:150,]

# 1) run PCA on reference breed genotype set
geno_pca <- prcomp(reference_breed,center = TRUE, scale. = FALSE)

# geno_pca$rotation is the eigen vector matrix
dim(geno_pca$rotation) # 100  50
# 50 is the principal component number 

# 2) run PCA reconstruction on test breed genotype set
mu = colMeans(reference_breed)  
# center test breed
test_breed_1_center <- scale(test_breed_1, center = mu, scale = FALSE)
# project test breed genotype set to the 50 principal components of reference breed
test_breed_1_pca_score <- as.matrix(test_breed_1_center) %*% geno_pca$rotation
# reproject breed to the orignal centered space
test_breed_1_reproject <- test_breed_1_pca_score %*% t(geno_pca$rotation)
# center back reconstruct test
reconstruct_test_breed_1 <- scale(test_breed_1_reproject, center = -mu, scale = FALSE)


# Finally, we can verify that reference set can rebuild to itself
# geno_pca$x is the reference set PCA score, we can directly reconstruct
reference_breed_reproject <- geno_pca$x %*% t(geno_pca$rotation)
# center back
reconstruct_reference_breed <- scale(reference_breed_reproject, center = -mu, scale = FALSE)
# This is exactly the same set as reference_breed




