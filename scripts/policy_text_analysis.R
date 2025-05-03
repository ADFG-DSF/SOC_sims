# Find the most common word and word pairs the SSFP and the escapement goal policy.

# Author: Adam
# Version: 2025-04-30

# Packages
packs <- c("tidyverse", "stringr", "tidytext", "readtext", "ggraph") #"textdata", 
lapply(packs, require, character.only = TRUE)
  
# Parameters
# Not Applicable

# ============================================================================


# Import text -------------------------------------------------------------

input <- 
  readtext::readtext(".\\5 AAC 39.222-5 AAC 39.223.docx")

words <- 
  input %>%
  unnest_tokens(output = word, input = text)
  
biwords <- 
  input %>%
  unnest_tokens(bigram, text, token = "ngrams", n = 2)


# Single word analysis -----------------------------------------------------

# * Frequent words -----------------------------------------------------

topwords <-
  words %>%
  mutate(word = ifelse(word == "stocks", "stock", word)) %>% #combine a plural
  filter(word != "means") %>% #delete since it's not informative
  anti_join(stop_words) %>% #stop words is a built in data set of the most common words
  count(word) %>%
  top_n(11) #11 so we can drop one for the stop criteria

topwords %>%
  mutate(stop = min(median(n), min(n))) %>%
  filter(n > stop) %>%
  ggplot(aes(x = n, y = reorder(word, n))) +
    geom_bar(alpha = 0.5, stat = "identity", show.legend = FALSE) +
    xlab("Number of Occurrences") +
    ylab("Most Common Words") +
    scale_y_discrete(labels = function(x) gsub("__.+$", "", x))


# Frequent word pairs -----------------------------------------------------

toppairs <- 
  biwords %>%
  separate(bigram, c("word1", "word2"), sep = " ") %>%
  filter(!word1 %in% stop_words$word,
         !word2 %in% stop_words$word,
         !is.na(word1)) %>%
  count(word1, word2, sort = TRUE)  %>%
  top_n(11)

toppairs %>%
  mutate(word = paste0(word1, " ", word2),
         stop = min(median(n), min(n))) %>%
  mutate(word = ifelse(word == "salmon stocks", "salmon stock", word)) %>% #combine a plural
  filter(word != "means") %>% #delete since it's not informative
  filter(n > stop) %>%
  ggplot(aes(x = n, y = reorder(word, n))) +
    geom_bar(alpha = 0.5, stat = "identity", show.legend = FALSE) +
    xlab("Number of Occurrences") +
    ylab("Most Common Word Pairs") +
    scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) 

# * Network analysis --------------------------------------------------------

GC_network <- 
  biwords %>%
  mutate(bigram = ifelse(bigram == "salmon stocks", "salmon stock", bigram)) %>% #combine a plural
  separate(bigram, c("word1", "word2"), sep = " ") %>%
  filter(!word1 %in% stop_words$word,
         !word2 %in% stop_words$word) %>%
  drop_na(word1, word2) %>%
  count(word1, word2, sort = TRUE)  %>%
  top_n(30) %>%
  mutate(word = paste0(word1, " ", word2),
         stop = min(median(n), min(n))) %>%
  filter(n > stop) %>%
  select(word, n) %>%
  igraph::graph_from_data_frame() #function assumes edges are first 2 columns.

# draw a network graph
set.seed(10) # 76
ggraph(GC_network, layout = "fr") +
  geom_edge_link(show.legend = FALSE, alpha = .5) +
  geom_node_point(color = "#0052A5", size = 3, alpha = .5) +
  geom_node_text(aes(label = name), vjust = 2) +
  ggtitle("Word Network in DSF Gotham Culture survey Responses") +
  theme_void() 
  #theme(plot.title = element_markdown())

