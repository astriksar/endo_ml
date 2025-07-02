# DS project 2024
# Annotation of Methylated Regions, Duplicate removal, merging of counts for the same ensemblid
# Julia Lapucha, Astghik Sarukhanyan, Caroline Forsythe


library(tidyverse)
library(stringr)
library(regioneR)
library(biomaRt)

# Read the count matrix from a text file (raw data, aprox. 3 Mil Regions for 77 Samples)
count_matrix <- read.table("D:/GIT/Master/DataScience/Project/Data/all_Methylation/GSE134052_countdata_mbdseq.txt")

# Check format
head(count_matrix)
row_number <- dim(count_matrix)

#filtering out rows that are 0 only
count_matrix <- count_matrix[rowSums(count_matrix[,1:row_number[2]] > 0) > 0,] 
idRow<- as.data.frame(row.names(count_matrix))

#preparing biomart query - filtering
names(idRow)[1] <- "chromosomes"
idRow[c("chr", "start","end")]<-str_split_fixed(idRow$chromosomes, '_', 3)
idRow$query = paste(gsub("chr",'',idRow$chr),idRow$start,idRow$end, sep = ":")

##################################### Biomart Query ##################################
mart = useMart(
  'ENSEMBL_MART_ENSEMBL',
  host = 'ensembl.org',
  dataset = 'hsapiens_gene_ensembl')

listAttributes(mart, page="feature_page")

# !!! the next code chunk runs approx. 3 hours!
out = getBM(
  attributes = c("chromosome_name",'start_position','end_position','ensembl_gene_id', 'external_gene_name', 'gene_biotype'),
  filters = 'chromosomal_region',
  values = idRow$query,
  mart = mart)
#
out <- out %>% mutate(region= paste("chr",chromosome_name,"_",start_position,"_",end_position, sep=""))


# write.csv(out, "./output/methylationAnnotationwithRegionandChr_allsamples.csv", row.names = F)
# if query saved as csv, it can be loaded at this point:
# out <- read.csv( "./output/Annotation/77Samples/biomart_output.csv")

############################################ Innitial Mapping of Query and regions ###############################

dfcount <- cbind(region = rownames(count_matrix), count_matrix)
rownames(dfcount) <- 1:nrow(dfcount)


# take start and end from idRow and check which out entries start end end there
idRow$start<- as.numeric(idRow$start)
idRow$end<- as.numeric(idRow$end)

out2<-out %>% mutate(Smatch=TRUE,Ematch=TRUE, chr=paste("chr",chromosome_name, sep=""))
chromList <- unique(out2$chr)
interdf<- data.frame()

# map and check overlap of the biomart genes with the chromosmal regions in our dataset 
# might run a bit
for(i in chromList){
  temp_count = idRow %>% filter(chr==i)
  temp_out = out2 %>% filter( chr==i)

  start_join <-temp_count %>%
    left_join(
      temp_out,
      by = join_by(between(y$start_position, x$start, x$end)),
      multiple = "all"
    ) %>%
    mutate(Smatch = !is.na(Smatch))

  end_join <-temp_count %>%
    left_join(
      temp_out,
      by = join_by(between(y$end_position, x$start, x$end)),
      multiple = "all"
    ) %>%
    mutate(Ematch = !is.na(Ematch))
  tempdf<- rbind(start_join, end_join)
  tempdf <- tempdf %>% filter(Smatch==T|Ematch==T)
  interdf <- rbind(interdf,tempdf)
}

head(interdf)
interdf %>% filter(Smatch!=Ematch) # 0, good the genes lie wholly in the regions

#write.csv(interdf, "./output/mappedAnnotation_allsamples.csv", row.names = F)

smallinderdf <- interdf %>% dplyr::select(chromosomes, ensembl_gene_id)

count_matrix$chromosomes <- row.names(count_matrix)

# join the regions
annotatedmatrix <- smallinderdf%>%
  left_join(
    count_matrix,
    by = join_by(chromosomes),
    multiple = "all"
  )
# remove duplicates
annotatedmatrix2 <-  annotatedmatrix[!duplicated(annotatedmatrix), ]

#write.csv(annotatedmatrix3, "./output/annotedMethylationCounts_allsamples.csv")
# cant assign ensemble ids as row names because there are doubles

############################################# Merge Rows with same IDs

# save duplicate Ids
duplicate_ids <-annotatedmatrix2 |>
  add_count(ensembl_gene_id) |>
  filter(n > 1)

# order
orderedDup <- duplicate_ids[order(duplicate_ids$ensembl_gene_id, 
                                  decreasing = TRUE), ]
orderedDup[c("chr", "start","end")]<-str_split_fixed(orderedDup$chromosomes, '_', 3)
# save index as column and reset index!

orderedDup$matrix_rownames<-row.names(orderedDup)
rownames(orderedDup) <- NULL
i=1
# loop to check regional overlap of same ids
# prints only if there is actual overlap
# no output if the regions are overlapping
# !!! This code can run 1-2 hours!!!
while(i < dim(orderedDup)[1]){
  tempid<-orderedDup$ensembl_gene_id[i]
  temp<- orderedDup %>% filter(ensembl_gene_id==tempid)
  if(dim(temp)[1]<=1){
    i <- i+dim(temp)[1]
    # exactly matches cases where rows were exactly duplicate
    next
  }else{
    if(dim(temp)[1]<2){print(paste("more than 2 overlaps found at: ",tempid, sep=""))}
    for (ii in 1:dim(temp)[1]){
      A <- temp[1,] %>% dplyr::select("chr","start","end")
      B <-temp[2,]%>% dplyr::select("chr","start","end")
      # print(length(overlapRegions(A,B)$type)) 
      # can be used to check output  ovf overlapping type
      # if terhe is no overlap there is no output
      # if A overlaps B putput is AinB etc.
      if (length(overlapRegions(A,B)$type)== 0){ # case when they are not overlaping
        # merge them no overlap -> done later since no overlap (faster)
      } else if(overlapRegions(A, B)$type != "equal"){print(paste("overlapping region found at: ",tempid, sep=""))}
    }
  }
  i <- i+dim(temp)[1]
}

# we had no overlap so we merge:
sample_dataset<- annotatedmatrix2  %>%  dplyr::group_by(ensembl_gene_id) %>% 
  dplyr::summarise(across(-chromosomes, sum))


# set ensembleids as rownames
sample_dataset_df <- data.frame(sample_dataset)
rownames(sample_dataset_df)<-sample_dataset_df$ensembl_gene_id
sample_dataset_df<-sample_dataset_df %>% dplyr::select(-(ensembl_gene_id))

write.csv(sample_dataset_df, "D:/GIT/Master/DataScience/Project/Data/
          annotedMethylationCounts_allsamples.csv")


#################################### Boxplot for data exploration
# load sample data for all 77 Patients
sample_data <- read.table("Data/all_Methylation/allMethylationSample_info.txt", header = TRUE, sep=",")
sample_dataset_df$ensembl_gene_id <- rownames(sample_dataset_df)

#convert into long format
count_matrix_melt<- reshape2::melt(sample_dataset_df, 
                                   variaD:/GIT/Master/DataScience/Project/ble.name = "Sample", 
                                   value.name = "Count")
# merghe patient info with counts
count_matrix_merge <- merge(count_matrix_melt, 
                            sample_data, by = "Sample") %>% arrange(Group)

count_matrix_merge <- count_matrix_merge %>%
  mutate(Sample = factor(Sample, levels = unique(Sample[order(Group)])))

# logtransform
count_matrix_merge$Count <- log(count_matrix_merge$Count)

boxplot_rep <- ggplot(count_matrix_merge, 
                      aes(x = reorder(Sample, Group), 
                          y = Count, 
                          color = Group)) +
  geom_boxplot() +
  geom_point(aes(color  = Group))+
  labs(title = "Distribution of Methylation log Counts",
       x = "Sample",
       y = "log Count") + 
  theme_gray()+
  theme(axis.text.x = element_text(face = "bold", 
                                   size = 10, angle = 90, hjust = 1)) +
  theme(axis.text.y = element_text(face = "bold", size =10)) 


# save plot, this will take a few minutes
ggsave(plot=boxplot_rep,
       filename = "D:/GIT/Master/DataScience/Project/boxplot_all.png",
       height = 5, width=9,dpi=320)
