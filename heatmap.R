library(tidyverse)
library(scales)
library(ggsci)
library(ggthemr)
library(ComplexHeatmap)
library(circlize)
library(GenomeInfoDbData)
setwd("C:\\Users\\wtbxs\\Documents\\work\\rice\\search_new\\manuscript\\figure\\Figure5")

data <- read.delim("./blast.txt",header = T)
data(speciesMap)
data(validTaxIds)
data(specData)
find_tax <- function(x){
  tax <- filter(speciesMap, species==x)
  return(tax$taxon)
}
species <- data %>% gather(key="species",value="match_sign",Aegilops_tauschii:Zea_mays) %>%
  mutate(species_c = gsub("[_]"," ", species)) %>% distinct(species_c) %>% 
  filter(!grepl(pattern = 'abinitio', species_c))

species2tax <- species %>%
  mutate(taxon = lapply(species_c, find_tax))
  
species2tax[species2tax$species_c=='Beta Vulgaris',]$taxon <- 161934   
species2tax[species2tax$species_c=='Hodeum vulgare',]$species_c <- 'Hordeum Vulgare'
species2tax[species2tax$species_c=='Hordeum Vulgare',]$taxon <- 4513
species2tax[species2tax$species_c=='Oryza indica',]$species_c <- 'Oryza sativa Indica Group'
species2tax[species2tax$species_c=='Oryza sativa Indica Group',]$taxon <- 39946 

# filter Oryza sativa
species2tax <- species2tax %>% filter(species_c!='Oryza sativa')


find_genus <- function(x){
  tax <- filter(specData, tax_id==x)
  return(paste(unique(as.character(tax$genus)),collapse = ";"))
}
species2genus <- species2tax  %>%
  mutate(genus = lapply(taxon, find_genus)) %>%
  separate(genus,c("first_g","second_g","third_g"),"[;]")

compare_v <- function(x,y){
  oryza <- species2taxtree[species2taxtree$genus=='Oryza',][-1]
  if(y==oryza[as.numeric(x)]){
    return(T)
  }else{
    return(F)
  }
}

species2taxtree <- read.delim("./Lineage.txt") %>% arrange(lineage.abbreviated.) %>%
  separate(lineage.abbreviated.,as.character(c(seq(1:17))),"; ")
species2taxtree[is.na(species2taxtree)] <- ""

species2distance <- species2taxtree %>% gather(key="pos",value="value",-genus) %>%
  group_by(genus) %>%
  mutate(sign=mapply(compare_v,`pos`,`value`)) %>%
  filter(sign==F) %>% 
  group_by(genus) %>%
  summarise(dist=min(as.numeric(`pos`)))

sup <- data.frame(genus=c("Leersia","Oryza"),dist=c(0,0))
species2distance <- rbind(species2distance,sup)

pep_detail <- as.data.frame(read.delim("./final_table_novel_peptide_index.txt",header=T,stringsAsFactors = F))
pep_ann <-
  pep_detail %>% distinct(Peptide, multiple_sign, Annotation, Associate_Known_Gene) %>%
  filter(Annotation == 'IRGSP-1.0.30') %>% mutate(intergenic_sign = lapply(Associate_Known_Gene, function(x) {
    if(x == 'Intergenic'){return('Intergenic')}else{return("Intragenic")}
  }))

multiple_sign <- pep_ann %>% distinct(Peptide,multiple_sign) 
intergenic_sign <-
  pep_ann %>% distinct(Peptide, intergenic_sign) %>% group_by(Peptide) %>%
  summarise(sign = paste(sort(unique(unlist(intergenic_sign))),collapse=";")) %>%
  distinct(Peptide, sign)
length_dis <- intergenic_sign %>% select(Peptide) %>%
  mutate(length=str_length(Peptide))





data <- left_join(data, length_dis, by = c("Novel_Peptide"="Peptide"))
data <- left_join(data, multiple_sign, by = c("Novel_Peptide"="Peptide"))
data <- left_join(data, intergenic_sign, by = c("Novel_Peptide"="Peptide"))

clean_data <- data %>% 
  filter(multiple_sign=='unique') %>%
  select(-contains('abinitio'),-Oryza_sativa) %>%
  arrange(length) %>% 
  select(-Novel_Peptide,-Hit.Count,-length,-multiple_sign,-sign)

row_names <- data %>% 
  filter(multiple_sign=='unique') %>%
  select(-contains('abinitio'),-Oryza_sativa) %>%
  arrange(length) %>% 
  select(Novel_Peptide)
rownames(clean_data) <- as.character(row_names$Novel_Peptide)

ha_row = rowAnnotation(
   Length_Distribution_Of_Peptides = row_anno_barplot(x = as.numeric((data %>% 
                            filter(multiple_sign=='unique') %>%
                            arrange(as.numeric(length)) %>%
                            select(length))$length)),
   width = unit(2, "cm"),show_annotation_name=T
)


species2genus_with_dis <- left_join(species2genus, species2distance, by=c("first_g"="genus")) 
species2genus_with_dis <- left_join(species2genus_with_dis, species2taxtree, by=c("first_g"="genus"))


Genus<-as.character((species2genus_with_dis %>% select(`13`))$`13`)
Genus[Genus==""]<-'Unknown'

ha_col <-
  columnAnnotation(df = data.frame(Genus = Genus[rev(hclust(dist(unname(t(
    clean_data
  ))))$order)]), height = unit(0.2, "grobwidth"))

#svg("./test.svg",16,9)
#png("./test.png",4096,4096,res=300)
Heatmap(clean_data, 
#       cluster_columns = hclust(dist(unname(t(clean_data)))), 
        cluster_columns = T,
       cluster_rows = F, show_row_names=F,
       split = (data %>% filter(multiple_sign=='unique') %>% arrange(length) %>% select(sign))$sign,
       col = c("grey","black"), show_heatmap_legend = FALSE
#       column_order = names(clean_data)[rev(hclust(dist(unname(t(clean_data))))$order)]
) + ha_row
#dev.off()

# 增加长度，位置，chr等信息作为附加标注


