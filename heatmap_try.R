library(tidyverse)
library(scales)
library(ggsci)
library(ggthemr)
library(ComplexHeatmap)
library(circlize)
setwd("C:\\Users\\wtbxs\\Documents\\work\\rice\\search_new\\manuscript\\figure\\Figure5")

data <- read.delim("./blast.txt",header=T,stringsAsFactors = F)

species2tax <- read.delim("./species2taxfull.txt",stringsAsFactors = F)

#species2tax[species2tax$species_c=='Hodeum vulgare',]$species_c <- 'Hordeum Vulgare'
#species2tax[species2tax$species_c=='Oryza indica',]$species_c <- 'Oryza sativa Indica Group'

species2tax_detail <- species2tax %>% separate(taxonomy,as.character(c(seq(1:22))),"; ") 
  

pep_detail <-
  as.data.frame(
    read.delim(
      "./final_table_novel_peptide_index.txt",
      header = T,
      stringsAsFactors = F
    )
  )
pep_ann <-
  pep_detail %>% distinct(Peptide, multiple_sign, Annotation, Associate_Known_Gene) %>%
  filter(Annotation == 'IRGSP-1.0.30') %>% mutate(intergenic_sign = lapply(Associate_Known_Gene, function(x) {
    if (x == 'Intergenic') {
      return('Intergenic')
    } else{
      return("Intragenic")
    }
  }))
multiple_sign <- pep_ann %>% distinct(Peptide, multiple_sign)
intergenic_sign <-
  pep_ann %>% distinct(Peptide, intergenic_sign) %>% group_by(Peptide) %>%
  summarise(sign = paste(sort(unique(
    unlist(intergenic_sign)
  )), collapse = ";")) %>%
  distinct(Peptide, sign)
length_dis <- intergenic_sign %>% select(Peptide) %>%
  mutate(length = str_length(Peptide))
data <- left_join(data, length_dis, by = c("Novel_Peptide"="Peptide"))
data <- left_join(data, multiple_sign, by = c("Novel_Peptide"="Peptide"))
data <- left_join(data, intergenic_sign, by = c("Novel_Peptide"="Peptide"))

clean_data <- data %>% 
  filter(multiple_sign=='unique') %>%
#  select(-contains('abinitio'),-Oryza_sativa) %>%
  arrange(length) %>% 
  select(-Novel_Peptide,-Hit.Count,-length,-multiple_sign,-sign)

row_names <- data %>% 
  filter(multiple_sign=='unique') %>%
#  select(-contains('abinitio'),-Oryza_sativa) %>%
  arrange(length) %>% 
  select(Novel_Peptide)
rownames(clean_data) <- as.character(row_names$Novel_Peptide)

col_names <- colnames(clean_data)
col_names_clean <- gsub("_"," ",col_names)
colnames(clean_data) <- col_names_clean

ha_row = rowAnnotation(
  "Length Distribution" = row_anno_barplot(x = as.numeric((
    data %>%
      filter(multiple_sign ==
               'unique') %>%
      arrange(as.numeric(length)) %>%
      select(length)
  )$length
  ), gp = gpar(col = pal_npg("nrc")(10)[9])),
  width = unit(2, "cm"),
  show_annotation_name = T
)

col_names <- colnames(data)[c(-1,-2,-62,-63,-64)]
col_names_detail <- as.data.frame(col_names) %>% 
  mutate(clean_name=gsub("_"," ",sub("_abinitio","",col_names)))

species2tax_detail_sort <- left_join(col_names_detail, species2tax_detail, by=c("clean_name"="species_c")) %>%
  gather(key="pos",value="genus",`1`:`22`)

Genus <- filter(species2tax_detail_sort,pos==16)$genus
#[hclust(dist(unname(t(clean_data))))$order]
Genus[!Genus=='Poaceae']<-'Others'
Genus[is.na(Genus)]<-'Others'
#Genus_f <- factor(Genus,levels = c("","Poaceae"),order = T)
Genus1 <- filter(species2tax_detail_sort,pos==18)$genus
Genus1[!Genus1=='Oryzoideae']<-'Others'
Genus1[is.na(Genus1)]<-'Others'
ha_col <-
  columnAnnotation(df = data.frame(Oryzoideae = Genus1,Poaceae = Genus),
                   col = list(Oryzoideae = c("Oryzoideae" = pal_npg("nrc")(10)[1], "Others" = "white"),
                              Poaceae = c("Poaceae" = pal_npg("nrc")(10)[4], "Others" = "white")),
                   show_annotation_name = T)

#svg("./test.svg",16,9)
png("./test.png",4096,4096,res=300)
Heatmap(clean_data, 
        #cluster_columns = hclust(dist(unname(t(clean_data)))), 
        cluster_columns = T,
        cluster_rows = F, show_row_names=F,
        split = (data %>% filter(multiple_sign=='unique') %>% arrange(length) %>% select(sign))$sign,
        col = c("white",pal_npg("nrc")(10)[8]), show_heatmap_legend = FALSE,
        top_annotation = ha_col
        #column_order = names(clean_data)[hclust(dist(unname(t(clean_data))))$order]
) + ha_row
dev.off()

