#install packages
install.packages("tidyverse")
install.packages("plyr")
install.packages("dplyr")
install.packages("readxl")
install.packages("writexl")
install.packages("factoextra")
install.packages("gplots")
install.packages("ggrepel")
install.packages("pheatmap")

#getwd() - get working directory
#setwd() - set working directory

#load packages
library(tidyverse)
library(plyr)
library(dplyr)
library (readxl)
library(writexl)
library(factoextra)
library(gplots)
library(ggrepel)
library(pheatmap)

#Proteomics data file upload, will not upload if the file is currently open
excel_sheets("mydataset.xlsx")
my_data= excel_sheets("mydataset.xlsx") %>% 
  map(~read_xlsx("mydataset.xlsx",.))
my_data

#Principle Component Analysis (PCA) to visualize the dataset
#For PCA use the transposed data on Sheet 2
#Remove the row with "Protein_Name"
#Add a column "Outcome" to classify the two groups in your dataset
R.PCA= prcomp(my_data[[2]][,3:100], scale=TRUE)
fviz_pca_ind(R.PCA, col.ind=my_data[[2]]$Outcome)
#Save the PCA plot as .png/.jpeg file

#Visualize the proteomics data in Sheet 1
summary(my_data[[1]])
View(my_data[[1]])

#Remove the missing values (Although no missing proteins in this data)
dat_norm=my_data[[1]] %>% drop_na()
#Number of proteins in original data
my_data[[1]]%>% dplyr::summarise(number_of_proteins=n())
#Number of proteins without missing values
dat_norm %>% dplyr::summarise(number_of_proteins=n())

#Select columns and log transform the data
dat_log= dat_norm %>%  
  select(-c(Protein_ID,Protein_Name)) %>% #to save time and space we just run these two commands together 
  log2()

#Bind columns to create log transformed data frame
dat_log_bind= bind_cols(dat_norm[,c(1:2)], dat_log)
view(dat_log_bind)
summary(dat_log_bind)

#we need to check if the data is actually normal in order to apply t-test

#Save log transformed data in excel file
write_xlsx(dat_log_bind, "Log2data.xlsx") #specify your path

# Creating a T-test function for multiple experiments
t_test= function(dt,grp1,grp2){
  # Subset Total Stroke Case group and convert to numeric
  x= dt[grp1] %>% unlist %>% as.numeric()
  # Subset Healthy Control group and convert to numeric
  y= dt[grp2] %>% unlist %>% as.numeric()
  # Perform t-test using the mean of x and y
  result= t.test(x, y)
  # Extract p-values from the results
  p_vals= tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
} 

# Apply t-test function to data using plyr adply
#  .margins = 1, slice by rows, .fun = t_test plus t_test arguments
dat_pvals= plyr::adply(dat_log_bind,.margins = 1, .fun = t_test, grp1 = c(3:12), grp2 = c(13:22)) %>% as_tibble()

#Check the t-test function created above by performing t-test on one protein
t.test(as.numeric(dat_log_bind[1,3:12]), as.numeric(dat_log_bind[1,13:22]))$p.value

# Plot histogram of p-values
dat_pvals %>% 
  ggplot(aes(p_val)) + 
  geom_histogram(binwidth = 0.05, 
                 boundary = 0.5, 
                 fill = "darkblue",
                 colour = "white") +
  xlab("p-value") +
  ylab("Frequency") +
  theme_minimal()

#Bind columns and add p-values to the "dat_log_bind" variable
#Store the data in "dat_combine" variable
dat_combine= bind_cols(dat_log_bind, dat_pvals[,23])
view(dat_combine)
summary(dat_combine)

#Calculating log-fold change
dat_fc= dat_combine %>% 
  group_by(Protein_ID) %>% 
  dplyr::mutate(mean_cases= mean(c(CASE1: CASE10)),
         mean_controls= mean(c(CON1: CON10)),
         log_fc = mean_cases - mean_controls,
         log_pval = -1*log10(p_val))

#Save final data with log-fold change & p-values in excel file
write_xlsx(dat_fc, "finaldataset.xlsx")

#Volcano plot of log-fold change on x-axis and log p-value on y-axis
VP= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val))) + geom_point() + theme_minimal()
VP 
#Add vertical lines for Log2 FC and a horizontal line for p-value threshold
VP2= VP + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP2
#Add a column for non-differential proteins
dat_fc$diffexpressed= "NO"
#Set Log2 FC and p-value cut-offs in the new column
dat_fc$diffexpressed[dat_fc$log_fc>0.58 & dat_fc$p_val<0.05] <- "UP"
dat_fc$diffexpressed[dat_fc$log_fc< -0.58 & dat_fc$p_val<0.05] <- "DOWN"
#Re-plot but this time color the points with "diffexpressed"
VP3= ggplot(data= dat_fc, aes(x=log_fc, y=-log10(p_val), col= diffexpressed)) + geom_point() + theme_minimal()
VP3
#Add lines as before..
VP4= VP3 + geom_vline(xintercept = c(-0.58, 0.58), col= "red") +
  geom_hline(yintercept = -log10(0.05), col="red")  
VP4
#Change point colors
VP5= VP4 + scale_color_manual(values= c("blue", "black", "red"))
mycolors= c("blue", "red", "black")  
names(mycolors) = c("DOWN", "UP", "NO")
VP5= VP4 + scale_color_manual(values=mycolors)
VP5
#Create a new column "proteinlabel" that will contain names of differential protein IDs
dat_fc$proteinlabel= NA
dat_fc$proteinlabel[dat_fc$diffexpressed != "NO"] <- dat_fc$Protein_ID[dat_fc$diffexpressed != "NO"]
ggplot(data=dat_fc, aes(x= log_fc, y= -log10(p_val), col= diffexpressed, label=proteinlabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text()

#Organize the labels using "ggrepel" package and "geom_text_repel" function
#Plot the Volcano plot using all layers used so far
ggplot(data= dat_fc, aes(x=log_fc, y= -log10(p_val), col= diffexpressed, label= proteinlabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline (xintercept = c(-0.58, 0.58), col="red") +
  geom_hline(yintercept = -log10(0.05), col="red")

#Proteins with significant observations
final_data=dat_fc %>%
  # Filter for significant observations
  filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58)) %>% 
  # Ungroup the data
  ungroup() %>% 
  # Select columns of interest
  select(Protein_ID,Protein_Name,mean_cases, mean_controls, log_fc, log_pval, p_val)
view(final_data)

#Save final data with list of differential proteins in excel file
write_xlsx(final_data, "diffproteins.xlsx")

#Filtering the dat_fc to only include differential proteins
dat_filt <- dat_fc %>% filter(log_pval >= 1.3 & (log_fc >= 0.58 | log_fc <= -0.58))

# Convert to matrix data frame
dat_matrix <- as.matrix.data.frame(dat_filt[,3:22]) 
# Name the rows with protein ids
row.names(dat_matrix) <- dat_filt$Protein_ID
# Transpose and scale the data to a mean of zero and sd of one
dat_scaled <- scale(t(dat_matrix)) %>% t()

# Transpose the matrix to calculate distance between experiments, row-wise
d1 <- dat_scaled %>% t() %>%
  dist(.,method = "euclidean", diag = FALSE, upper = FALSE)
# Calculate the distance between proteins row-wise 
d2 <- dat_scaled %>%
  dist(.,method = "euclidean", diag = FALSE, upper = FALSE)

# Set colours for heatmap, 25 increments
my_palette <- colorRampPalette(c("blue","white","red"))(n = 50)

# Plot pheatmap
dat_scaled %>% 
  pheatmap(.,
           fontsize = 7, cluster_cols = FALSE, cluster_rows = TRUE, 
           color= redblue(20), border_color = "black", 
           scale = "none")
