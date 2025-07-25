# Read in Historic SOC listingings as R objects

# Author: Adam Reimer
# Version: 2025-07-7

# Packages
packs <- c("tidyverse", "DT", "flextable")
lapply(packs, require, character.only = TRUE)

# Reads in Andrews table of SOC listings and delistings (current through the 2024 BOF cycle)
# ".\SOC_sims\SOC_report_outline_Mar28working.docx"
# formats for R and creates R tables and figures.
# --------------------------------------------------------------------------------------
# Current SOC listings
SOC_current <- data.frame(
  region = c(rep("Southeast", 7), rep("Central", 8), rep("AYK", 2), rep("Westward", 4)),
  stock = c("King Salmon River", "Stikine River", "Andrew Creek", "Taku River", "McDonald Lake", "Hugh Smith Lake", "Northern Southeast Outside Subregion",
            "McNeil River", "Chuitna River", "Theodore River", "Alexander Creek",  "Eastside Susitna River", "Nushagak River", "Kenai River - late run", "Mikfik Lake", 
            "Yukon River", "Norton Sound Subdistrict 5 & 6", 
            "Karluk River", "Ayakulik River", "Chignik River", "Chignik River - early run"),
  species = c(rep("Chinook", 4), rep("Sockeye", 2), "Chum", 
              "Chum", rep("Chinook", 6), "Sockeye",
              rep("Chinook", 2),
              rep("Chinook", 3), "Sockeye"),
  list_date = as.Date(
    c("2018-1-1", rep("2022-3-1", 3), "2018-1-1", rep("2024-10-1", 2),
      "2016-12-1", rep("2011-2-1", 3), "2020-2-1", "2023-3-1", rep("2023-10-1", 2),
      "2000-9-1", "2004-1-1", 
      "2011-1-1", "2020-1-1", "2023-2-1", "2022-3-1")
  ),
  type = c(rep("Management", 7),
           rep("Management", 8),
           rep("Yield", 2),
           rep("Management", 4)),
  last_review = c(rep(2024, 7),
                  rep(2023, 5), 2022, rep(2023, 2),
                  rep(2022, 2),
                  rep(2023, 2), rep(2022, 2))
)

#Delisted SOC
SOC_delisted <- data.frame(
  region = c(rep("Southeast", 6), rep("Central", 11), rep("AYK", 9), rep("Westward", 1)),
  stock = c("Hugh Smith Lake", "McDonald Lake", "Chilkat River", "Unuk River", "Chickamin River", "Klukshu River",
            "Anchor River", "Lewis River", "Sheep Creek", "Goose Creek", "Kvichak River", "Fish Creek", "Susitna (Yentna) River", "Willow Creek", "Goose Creek", "Kvichak River", "Kvichak River",
            "Norton Sound Subdistrict 1", "Toklat River", "Fishing Branch", "Yukon River - Summer Chum", "Kuskokwim River", "Kuskokwim River", "Yukon River - Fall Chum", "Norton Sound Subdistrict 2 & 3", "Norton Sound Subdistrict 1", 
            "Swanson Lagoon"),
  species = c(rep("Sockeye", 2), rep("Chinook", 3), "Sockeye",
              rep("Chinook", 4), rep("Sockeye", 3), rep("Chinook", 2), rep("Sockeye", 2),
              rep("Chum", 5), "Chinook", rep("Chum", 3),
              rep("Sockeye")),
  list_date = as.Date(
    c("2003-2-1", "2009-2-1", "2017-10-1", "2018-1-1", rep("2022-3-1", 2),
      "2001-11-1", "2011-2-1", rep("2014-1-1", 2), "2003-12-1", "2002-1-1", "2008-2-1", rep("2011-2-1", 2), "2000-9-1", "2009-10-1",
      rep("2000-9-1", 8), "2007-1-1",
      "2013-2-1")
  ),
  type = c(rep("Management", 6),
           rep("Management", 5), rep("Yield", 6),
           rep("Management", 4), rep("Yield", 5),
           "Management"),
  delist_date = as.Date(
    c("2006-1-1", "2012-2-1", rep("2024-10-1", 4),
      "2004-11-1", "2019-10-1", rep("2020-2-1", 2), "2009-10-1", "2005-1-1", rep("2020-2-1", 2), "2014-1-1", "2003-12-1", "2012-12-1",
      "2007-1-1", rep("2004-1-1", 2), rep("2007-1-1", 4), "2019-1-1", "2016-1-1",
      "2019-2-1")
  )
)
SOC_tab <-
  SOC_current %>%
  mutate(delist_date = as.Date(NA)) %>%
  rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
  select(-last_review) %>%
  mutate(years_listed = round(as.numeric(difftime(as.Date(ifelse(is.na(delist_date), today(), delist_date)), list_date, units = "days")) / 365, 1),
         year = as.numeric(format(list_date, format = "%Y")))

multiple_listings <- 
  SOC_tab %>%
  mutate(delisted = !is.na(delist_date)) %>%
  group_by(species, stock) %>%
  summarise(n_listings = n(),
            n_delistings = sum(delisted)) %>%
  filter(n_listings >1 & n_delistings < n_listings) %>%
  select(species, stock)

SOC_current %>%
  mutate(delist_date = as.Date(NA)) %>%
  rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
  select(-last_review) %>%
  mutate(stock = ifelse(stock %in% multiple_listings$stock & 
                          species %in% multiple_listings$species &
                          !is.na(delist_date), paste0(stock, "_temp"), stock)) %>%
  mutate(stock = gsub("(.*)_temp", "\\1", stock))
gsub("(.*)_temp", "\\1", c("dog_temp", "dog_t"))
#Table
# datatable version
SOC_current %>%
  mutate(delist_date = as.Date(NA)) %>%
  rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
  select(-last_review) %>%
  mutate(years_listed = round(as.numeric(difftime(as.Date(ifelse(is.na(delist_date), today(), delist_date)), list_date, units = "days")) / 365, 1),
         list_date = format(list_date, format = "%Y-%b"),
         delist_date = ifelse(is.na(delist_date), "ongoing", format(delist_date, format = "%Y-%b"))) %>%
  arrange(species, list_date) %>%
  mutate(species = as.factor(species),
         region = as.factor(region),
         type = as.factor(type)) %>%
  select(species, region, stock, type, list_date, delist_date, years_listed) %>%
  datatable(class = 'stripe hover compact', 
            filter = "top", 
            rownames = FALSE,
            colnames = c('Species' = 'species',
                         'Region' = 'region',
                         'Stock' = 'stock',
                         'SOC Type' = 'type',
                         'Listing Date' = 'list_date',
                         'Delisting Date' = 'delist_date',
                         'Years as a SOC' = 'years_listed')
            )
  
# flextable version
SOC_tab <-
  SOC_current %>%
  mutate(delist_date = as.Date(NA)) %>%
  rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
  select(-last_review) %>%
  mutate(years_listed = round(as.numeric(difftime(as.Date(ifelse(is.na(delist_date), today(), delist_date)), list_date, units = "days")) / 365, 1))
SOC_tab %>%
  group_by(stock, species) %>%
  mutate(list_date = format(list_date, format = "%Y"),
         delist_date = ifelse(is.na(delist_date), "-", format(delist_date, format = "%Y")),
         first_yr = min(as.numeric(list_date)), 
         sort = paste0(first_yr, stock)) %>%
  arrange(species, sort, list_date) %>% #wanted the table sorted by listing year but also group the same stock by the first year listed.
  select(species, stock, region, type, list_date, delist_date, years_listed) %>%
  flextable() %>%
  set_header_labels(species = "Species", 
                    stock = "Stock", 
                    region = "Region",
                    type = 'SOC \nType',
                    list_date = 'Listing \nYear', 
                    delist_date = 'Delisting \nYear', 
                    years_listed = 'Years \nListed') %>%
  merge_v(~ species + stock) %>%
  valign(j = c("species", "stock"), valign = "top") %>%
  footnote(i = c(7, 12, 13, 18), j = 2, 
           value = as_paragraph(
             c(" Willow, Goose and Sheep Creek SOC listing were consolidated into a single SOC listing (Eastside Susitna River) in 2020 after the escapement goals for these populations were consolidated into a single goal.")),
           ref_symbols = "a",
           part = "body") %>%
  footnote(i = 44, j = 2, 
           value = as_paragraph(
             c(" Chignik River early run sockeye action plan adopted by the BOF in Febuary 2023")),
           ref_symbols = "b",
           part = "body") %>%
  autofit()


table(SOC_tab$type)
table(SOC_tab$region)
table(SOC_tab$species)
#Unique stocks listed
SOC_unique <- 
  SOC_tab %>%
  filter(!duplicated(SOC_tab[, c("species", "stock")])) %>%
  filter(!(stock %in% c("Goose Creek", "Sheep Creek", "Willow Creek"))) %>%
  select(species, stock)
dim(SOC_unique)
SOC_unique %>%
  group_by(species) %>%
  summarise(n = n())

#% delisted
table(is.na(SOC_tab$delist_date))
#% delisted by species
# accounts for stocks which changed type or stock designation but remained withing the SOC process
# Only consider stock listed for 3+ years (since those <3 had no chance to be delisted).
SOC_tab %>%
  filter(!(duplicated(SOC_tab[, c("species", "stock")]) & SOC_tab$stock != "McDonald Lake")) %>%
  filter(!(stock %in% c("Goose Creek", "Sheep Creek", "Willow Creek"))) %>%
  mutate(delisted = ifelse(is.na(delist_date), 0, 1)) %>% 
  filter(as.numeric(format(list_date, format = "%Y")) < 2022) %>%
  group_by(species) %>%
  summarise(pct_delisted = mean(delisted))


# Plot by species
intercept = as.numeric(format(today(), format = "%Y"))
# combine tables
SOC_plot <-
  SOC_current %>%
  mutate(delist_date = as.Date(NA)) %>%
  rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
  select(-last_review) %>%
  #expand East Susitna stocks since none of those stocks have ever been delisted and they have different list dates
  filter(stock != "Eastside Susitna River") %>%
  mutate(delist_date =
           as.Date(ifelse(stock %in% c("Sheep Creek", "Goose Creek", "Willow Creek") & delist_date == SOC_current[SOC_current$stock == "Eastside Susitna River", colnames(SOC_current) == "list_date"],
                          NA,
                          delist_date)),
         year = as.numeric(format(list_date, format = "%Y")),
         years_listed = round(as.numeric(difftime(as.Date(ifelse(is.na(delist_date), today(), delist_date)), list_date, units = "days")) / 365, 1))

SOC_plot %>%
  ggplot(aes(x = year, y = years_listed, color = species, shape = type)) +
  geom_jitter(size = 3, width = 0.1, height = 0.4) +
  geom_abline(intercept = intercept, 
              slope = -1,
              linewidth = 15,
              color = "lightgrey",
              alpha = 0.6) +
  scale_y_continuous(name = "Years as a SOC", 
                     breaks = seq(0, 28, 3),
                     minor_breaks = NULL) +
  scale_shape_manual(values = c("Yield" = 1, "Management" = 16)) +
  labs(x = "Year Listed")

#Plot by type
SOC_plot %>%
  ggplot(aes(x = year, y = years_listed, color = type)) +
  geom_jitter(size = 3) +
  geom_abline(intercept = intercept, 
              slope = -1,
              linewidth = 15,
              color = "lightgrey",
              alpha = 0.6) +
  scale_y_continuous(name = "Years as a SOC", 
                     breaks = seq(0, 28, 3),
                     minor_breaks = NULL) +
  labs(x = "Year Listed")

#Plot by region
SOC_plot %>%
  ggplot(aes(x = year, y = years_listed, color = region)) +
  geom_jitter(size = 3) +
  geom_abline(intercept = intercept, 
              slope = -1,
              linewidth = 15,
              color = "lightgrey",
              alpha = 0.6) +
  scale_y_continuous(name = "Years as a SOC", 
                     breaks = seq(0, 28, 3),
                     minor_breaks = NULL) +
  labs(x = "Year Listed")


  
  
  test <-  
    SOC_current %>%
    mutate(delist_date = as.Date(NA)) %>%
    rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
    select(-last_review) %>%
    mutate(stock = ifelse(stock == "McDonald Lake" & is.na(delist_date), "McDonald Lake2", stock)) %>% #rename since listings were separated in time
    group_by(species, stock) %>%
    arrange(species, stock, list_date) %>%
    summarize(list_date = min(list_date), delist_date = max(delist_date),
              type = paste(type, collapse = " -> ")) %>%
    mutate(stock = ifelse(stock %in% c("Sheep Creek", "Goose Creek", "Willow Creek"), "temp", stock)) %>% #create a temporary group to get the type correct
    group_by(species, stock) %>%
    arrange(species, stock, list_date) %>%
    summarize(list_date = min(list_date), delist_date = max(delist_date),
              type = paste(type, collapse = ", ")) %>%
    mutate(type = ifelse(stock == "temp", paste0("(", type, ")"), type),
           stock = ifelse(stock == "temp", "Eastside Susitna River", stock)) %>%
    group_by(species, stock) %>%
    arrange(species, stock, list_date) %>%
    summarize(list_date = min(list_date), delist_date = max(delist_date),
              type = paste(type, collapse = " -> ")) %>% 
    mutate(stock = ifelse(stock == "McDonald Lake2", "McDonald Lake", stock),
           year = as.numeric(format(list_date, format = "%Y")),
           years_listed = round(as.numeric(difftime(as.Date(ifelse(is.na(delist_date), today(), delist_date)), list_date, units = "days")) / 365, 1))
  
test %>%
  mutate(delist = ifelse(is.na(delist_date), 0, 1)) %>%
  filter(as.numeric(format(list_date, format = "%Y")) < 2022) %>%
  summarise(delist = mean(delist),
            n = n())

test %>%
  ggplot(aes(x = year, y = years_listed, color = species)) +
  geom_jitter(size = 3, widthe = 0.2, height = 0.5) +
  geom_abline(intercept = intercept, 
              slope = -1,
              linewidth = 15,
              color = "lightgrey",
              alpha = 0.6) +
  scale_y_continuous(name = "Years as a SOC", 
                     breaks = seq(0, 28, 3),
                     minor_breaks = NULL) +
  labs(x = "Year Listed")

SOC_plot <-
  SOC_current %>%
  mutate(delist_date = as.Date(NA)) %>%
  rbind(SOC_delisted %>% mutate(last_review = NA)) %>%
  select(-last_review) %>%
  #expand East Susitna stocks since none of those stocks have ever been delisted and they have different list dates
  filter(stock != "Eastside Susitna River") %>%
  mutate(delist_date =
           as.Date(ifelse(stock %in% c("Sheep Creek", "Goose Creek", "Willow Creek") & delist_date == SOC_current[SOC_current$stock == "Eastside Susitna River", colnames(SOC_current) == "list_date"],
                          NA,
                          delist_date)),
         year = as.numeric(format(list_date, format = "%Y")),
         years_listed = round(as.numeric(difftime(as.Date(ifelse(is.na(delist_date), today(), delist_date)), list_date, units = "days")) / 365, 1))

SOC_plot %>%
  mutate(current = ifelse(is.na(delist_date), "Yes", "No")) %>% arrange(species, current, years_listed)
  ggplot(aes(x = years_listed, fill = current)) +
  geom_histogram(binwidth = 1) +
  facet_grid(species ~ .) +
  scale_x_continuous(name = "Years as a SOC", 
                     breaks = seq(0, 28, 3),
                     minor_breaks = NULL) +
  labs(y = "# of stocks")
