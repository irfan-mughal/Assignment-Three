# Ques 1) Download file "gene_expression.tsv", read in r script and making gene accession number as row number
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv",
              destfile = "gene_expression.tsv")
x<- read.table("gene_expression.tsv")
x<-read.table("gene_expression.tsv", header = TRUE, row.names = 1)
head(x)
str(x)

# Ques 2) Make mean other column and show values of first 6 genes
x$mean<- rowMeans(x)
head(x)

#Ques 3) list of 10 genes with highest mean
order(-x$mean)
ord<- x[order(-x$mean), ]
head(ord, 10)

#Ques 4) mean less than 10
subset(x, mean<10)

#Ques 5) hist of mean
hist(x$mean, main = "Mean values", xlab = "Mean")

#Ques 6) Download file "growth_data.csv" and read into r script and specify column names
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv",
              destfile = "growth_data.csv")
y<- read.csv("growth_data.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(y)
str(y)
colnames(y)

#Ques 7) calculate mean and standart deviation of tree circumference at both sites at start and end
# subsetof northeast(NE) site 
subset(y, Site == "northeast")
NE <- subset(y, Site == "northeast")
head(NE)
str(NE)
#At start
mean(NE$Circumf_2004_cm)
sd(NE$Circumf_2004_cm)
#At end
mean(NE$Circumf_2019_cm)
sd(NE$Circumf_2019_cm)

#subset of southwest(SW) site
subset(y, Site == "southwest")
SW<- subset(y, Site == "southwest")
head(SW)
str(SW)
#At start
mean(SW$Circumf_2004_cm)
sd(SW$Circumf_2004_cm)
#At end
mean(SE$Circumf_2019_cm)
sd(SE$Circumf_2019_cm)

#Ques 8) Boxplot of tree circumference at start and end of both sites
#At northeast site
boxplot(NE$Circumf_2004_cm, NE$Circumf_2019_cm, main = "Boxplot of Tree circumference at northeast", ylab ="circumference of tree", names = c("NE_2004", "NE_2019"))
#At southwest site
boxplot(SW$Circumf_2004_cm, SW$Circumf_2019_cm, main = "Boxplot of Tree circumference at southwest", ylab ="circumference of tree", names = c("SW_2004", "SW_2019"))

#Ques 9) growth of mean over past 10 years at each site
# at NE site
a<- NE$Circumf_2019_cm - NE$Circumf_2009_cm
mean(a)
head(NE)
str(NE)
#at SW site
b<-SW$Circumf_2019_cm - SW$Circumf_2009_cm
mean(b)
head(SW)
str(SW)

#Ques 10) t.test and wilcox.test
t.test(a, b)
wilcox.test(a, b)
