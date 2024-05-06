
# ---- DATA WRANGLING

fasta_to_tibble <- function (fasta_path) {
  seqs <- read_fasta(fasta_path, type = "AA")
  tibble_out <- tibble(label = names(seqs), sequences = seqs, aa_length = seq_nchar(seqs))
  tibble_out
}

make_accession_col <- function (tibble_in) {
  tibble_out <- tibble_in %>%
    mutate(accession = str_extract(label, "\\w+(?=\\|)"))
  tibble_out
}

make_accession_col2 <- function (tibble_in) {
  tibble_out <- tibble_in %>%
    mutate(accession = str_extract(label, "(?<=tr\\|)\\w+"))
  tibble_out
}

make_genus_col <- function(tibble_in){
  tibble_out <- tibble_in %>%
    mutate(GENUS = str_extract(`Taxonomic lineage (SPECIES)`, "\\w+(?=\\s)"))
  tibble_out
}

day_name_string <- function (name_string) {
  today <- Sys.Date()
  today_str <- format(today, format="%b%d")
  #random_number <- as.character(sample(1000:9999, 1))
  string_out <- str_interp("${today_str}_${name_string}")
  string_out
}

connect_data <- function (df, uniprot_file){
  uniprot_df <- read_tsv(uniprot_file)
  df_out <- inner_join(uniprot_df, df, by = c("Entry" = "accession"))
  df_out
}



filter_common_taxons <- function(active_df, ord_min, fam_min){
  
  common_orders <- active_df %>% 
    count(`Taxonomic lineage (ORDER)`) %>% 
    filter(n > ord_min) 
  orders <- common_orders$`Taxonomic lineage (ORDER)` 
  
  common_fams <- active_df %>% 
    count(`Taxonomic lineage (FAMILY)`) %>% 
    filter(n > fam_min) 
  fams <- common_fams$`Taxonomic lineage (FAMILY)` 
  
  popular <- active_df %>% 
    filter(`Taxonomic lineage (ORDER)` %in% orders) %>%
    filter(`Taxonomic lineage (FAMILY)` != "NA") %>%
    filter(`Taxonomic lineage (FAMILY)` %in% fams)
  popular
}

filter_common_taxons2 <- function(active_df, ord_min, fam_min){
  
  common_orders <- active_df %>% 
    count(order) %>% 
    filter(n > ord_min) 
  orders <- common_orders$order
  
  common_fams <- active_df %>% 
    count(family) %>% 
    filter(n > fam_min) 
  fams <- common_fams$family
  
  popular <- active_df %>% 
    filter(order %in% orders) %>%
    filter(family != "NA") %>%
    filter(family %in% fams)
  popular
}


  # ---- FILE EXPORT

export_accessions <- function (tibble_in, name_string, accession_col="accession") {
  unique_name <- day_name_string(name_string)
  tibble_in[,accession_col] %>%
    write_tsv(str_interp("Exported_Accessions_${unique_name}.txt"),
              col_names = FALSE)
}

write_clean_fasta <- function (df, outputName, accession_col="accession"){
  string_out <- str_interp("${outputName}.fasta")
  df %>%
    select(accession_col, sequences) %>%
    deframe() %>%
    write_fasta(string_out)
}

write_clean_fasta2 <- function (df, outputName){
  string_out <- str_interp("${outputName}.fasta")
  df %>%
    select(label, sequences) %>%
    deframe() %>%
    write_fasta(string_out)
}

# ---- GGplot Functions

lengths_graph <- function (tibble_in) {
  ggplot(tibble_in, aes(x=aa_length)) +
    geom_histogram(binwidth = 4, fill = "white", colour = "black") +
    theme_bw()
}


count_graph1 <- function (tibble_in, xVal, fVal) {
  ggplot(tibble_in, 
         aes(x = xVal, fill = fVal)) +
    geom_bar(stat = "count") +
    theme_bw()
}

count_graph2 <- function (df, xv, fv) {
  ggplot(df,aes(x = fct_infreq(xv), fill = fv)) + 
    geom_bar(stat = 'count') +  
    theme_bw()
}
    
count_graph3 <- function (tibble_in, xVal) {
  ggplot(tibble_in, 
         aes(x = fct_infreq(xVal))) +
    geom_bar(stat = "count") +
    theme_bw()
}


  

taxon_graph <- function(dat){
  taxon_graph <- ggplot(dat,
                        aes(x = `Taxonomic lineage (FAMILY)`,
                            fill = `Taxonomic lineage (ORDER)`)) +
    geom_bar() +
    coord_polar() +
    scale_y_continuous(trans='log2') + 
    theme_bw()
  taxon_graph
}

lengths_graph_anno <- function (tibble_in) {
  
  length_statistics <- summary(tibble_in$aa_length)
  min_label <- str_c(names(length_statistics[1]), ":", as.character(length_statistics[[1]]), sep=" ")
  max_label <- str_c(names(length_statistics[6]), ":", as.character(length_statistics[[6]]), sep=" ")
  median_label <- str_c(names(length_statistics[3]), ":", as.character(length_statistics[[3]]), sep=" ")
  mean_label <- str_c(names(length_statistics[4]), ":", round(length_statistics[[4]], digits=0), sep=" ")
  q1_label <- str_c(names(length_statistics[2]), ":", as.character(length_statistics[[2]]), sep=" ")
  q3_label <- str_c(names(length_statistics[5]), ":", as.character(length_statistics[[5]]), sep=" ")
  
  p <- lengths_graph(tibble_in)
  p +
    annotate("segment", x = length_statistics[[1]], xend = length_statistics[[1]], y = 0, yend = 10, colour = "blue") +
    annotate("segment", x = length_statistics[[6]], xend = length_statistics[[6]], y = 0, yend = 10, colour = "blue") +
    annotate("segment", x = length_statistics[[3]], xend = length_statistics[[3]], y = 0, yend = 10, colour = "green") +
    annotate("segment", x = length_statistics[[4]], xend = length_statistics[[4]], y = 0, yend = 10, colour = "orange") +
    #annotate("rect", xmin = length_statistics[[2]], xmax = length_statistics[[5]], ymin = 0, ymax = 100, alpha = 0.1) +
    annotate("text", x = (length_statistics[[1]] + 3), y = 11, label=min_label, colour="blue") +
    annotate("text", x = (length_statistics[[6]] - 3), y = 11, label=max_label, colour="blue") +
    annotate("text", x = (length_statistics[[3]] - 3), y = 11, label=median_label, colour="green") +
    annotate("text", x = (length_statistics[[4]] + 5), y = 9, label=mean_label, colour="orange")
}

# ---- GGtree Functions

