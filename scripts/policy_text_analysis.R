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

input222a <- readLines(".\\5 AAC 39.222(a).txt") %>% as.data.frame() %>% mutate(section = "222(a) -preable")
input222c <- readLines(".\\5 AAC 39.222(c).txt") %>% as.data.frame() %>% mutate(section = "222(c) - principles")
input222d <- readLines(".\\5 AAC 39.222(d).txt") %>% as.data.frame() %>% mutate(section = "222(d) - implementation")
input222f <- readLines(".\\5 AAC 39.222(f).txt") %>% as.data.frame() %>% mutate(section = "222(f) - definitions")
input223 <- readLines(".\\5 AAC 39.223.txt") %>% as.data.frame() %>% mutate(section = "223")
input_proposed <- readLines(".\\SOC_policy_language_April30.txt") %>% as.data.frame() %>% mutate(section = "222(d) - proposed")

input <- 
  do.call(rbind, list(input222a, input222c, input222d, input222f, input223, input_proposed)) %>%
  setNames(c("text", "section")) %>%
  filter(text != "")

words <- 
  input %>%
  unnest_tokens(output = word, input = text)

biwords <- 
  input %>%
  unnest_tokens(bigram, text, token = "ngrams", n = 2)

#dev.new()
# Single word analysis -----------------------------------------------------

# * Frequent words -----------------------------------------------------
topwords <-
  words %>%
  mutate(word = ifelse(word == "stocks", "stock", word)) %>% #combine a plural
  mutate(word = ifelse(word == "habitats", "habitat", word)) %>% #combine a plural
  filter(word != "means") %>% #delete since it's not informative
  anti_join(stop_words) %>% #stop words is a built in data set of the most common words
  group_by(section) %>%
  count(word)

#counts for any word you want
topwords[topwords$word == "phenotypic", ]

topwords %>%
  filter(section %in% c("222(c) - principles", 
                        "222(d) - implementation",
                        "222(d) - proposed", 
                        "222(f) - definitions")) %>% #retain section with many words
  top_n(10) %>% #11 so we can drop one for the stop criteria
  group_by(section) %>%
  mutate(stop = min(median(n), min(n))) %>%
  filter(n > stop) %>%
  ggplot(aes(x = n, y = reorder_within(word, n, section), fill = section)) +
  geom_bar(alpha = 0.5, stat = "identity", show.legend = FALSE) +
  facet_wrap(~ section, ncol = 2, scales = "free_y") +
  xlab("Number of Occurrences") +
  ylab("Most Common Words") +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x))


# * tf_idf ------------------------------------------------------------------
tf_idf <- 
  words %>%
  mutate(word = ifelse(word == "stocks", "stock", word)) %>% #combine a plural
  mutate(word = ifelse(word == "habitats", "habitat", word)) %>% #combine a plural
  mutate(word = ifelse(word == "concerns", "concern", word)) %>% #combine a plural
  filter(word != "means") %>% #delete since it's not informative
  anti_join(stop_words) %>%
  count(section, word, sort = TRUE) %>%
  bind_tf_idf(word, section, n) %>%
  arrange(section, desc(tf_idf))


tf_idf %>%
  # We need to sort the data in descending order so we can create the factors for each term
  arrange(desc(tf_idf)) %>%
  group_by(section) %>%
  top_n(11) %>%
  mutate(stop = min(median(n), min(n))) %>%
  filter(n > stop) %>%
  ggplot(mapping = aes(x = tf_idf, y = reorder_within(word, tf_idf, section), fill = section)) +
  geom_bar(alpha = 0.5, stat = "identity", show.legend = FALSE) +
  facet_wrap(~ section, scales = "free_y") +
  xlab("Weighted Freqency of Occurrence") +
  ylab("Most Common Words") +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x))

# Frequent word pairs -----------------------------------------------------

toppairs <- 
  biwords %>%
  separate(bigram, c("word1", "word2"), sep = " ") %>%
  filter(!word1 %in% stop_words$word,
         !word2 %in% stop_words$word,
         !is.na(word1)) %>%
  mutate(word1 = ifelse(word1 == "stocks", "stock", word1)) %>% #combine a plural
  mutate(word1 = ifelse(word1 == "habitats", "habitat", word1)) %>% #combine a plural
  mutate(word1 = ifelse(word1 == "concerns", "concern", word1)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "stocks", "stock", word2)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "habitats", "habitat", word2)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "concerns", "concern", word2)) %>% #combine a plural
  group_by(section) %>%
  count(word1, word2, sort = TRUE)  %>%
  top_n(6)

toppairs %>%
  group_by(section) %>%
  mutate(word = paste0(word1, " ", word2),
         stop = min(median(n), min(n))) %>%
  filter(n > stop) %>%
  ggplot(aes(x = n, y = reorder_within(word, n, section), fill = section)) +
  geom_bar(alpha = 0.5, stat = "identity", show.legend = FALSE) +
  facet_wrap(~ section, ncol = 2, scales = "free_y") +
  xlab("Number of Occurrences") +
  ylab("Most Common Word Pairs") +
  scale_y_discrete(labels = function(x) gsub("__.+$", "", x)) 

# * Network analysis --------------------------------------------------------

GC_network <- 
  biwords %>% 
  separate(bigram, c("word1", "word2"), sep = " ") %>%
  filter(!word1 %in% stop_words$word,
         !word2 %in% stop_words$word) %>%
  mutate(word1 = ifelse(word1 == "stocks", "stock", word1)) %>% #combine a plural
  mutate(word1 = ifelse(word1 == "habitats", "habitat", word1)) %>% #combine a plural
  mutate(word1 = ifelse(word1 == "concerns", "concern", word1)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "stocks", "stock", word2)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "habitats", "habitat", word2)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "concerns", "concern", word2)) %>% #combine a plural
  mutate(word1 = ifelse(word1 == "means", NA, word1)) %>% #combine a plural
  mutate(word2 = ifelse(word2 == "means", NA, word2)) %>% #combine a plural
  drop_na(word1, word2) %>%
  group_by(section) %>%
  count(word1, word2, sort = TRUE)  %>%
  top_n(6) %>%
  group_by(section) %>%
  mutate(word = paste0(word1, " ", word2),
         stop = min(median(n), min(n))) %>%
  filter(n > stop) %>%
  select(section, word, n) %>%
  igraph::graph_from_data_frame() #function assumes edges are first 2 columns.

# draw a network graph
set.seed(10) # 76
ggraph(GC_network, layout = "fr") +
  geom_edge_link(aes(edge_alpha = n, edge_width = n), show.legend = FALSE, alpha = .5) +
  geom_node_point(color = "#0052A5", size = 3, alpha = .5) +
  geom_node_text(aes(label = name), vjust = 2) +
  ggtitle("Word Network in 5 AAC 39.222-223") +
  theme_void() 
#theme(plot.title = element_markdown())

#dev.off()
