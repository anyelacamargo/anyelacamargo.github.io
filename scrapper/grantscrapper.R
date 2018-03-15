library(rvest)
library(stringr)
library(plyr)
library(dplyr)
library(ggvis)
library(knitr)


call <- read_html("bbsrc.html") %>%
  html_nodes("h3 a") %>%
  html_text()

batch_titles <- batches %>%
  html_nodes("href") %>%
  html_text()


url <- read_html("bbsrc.html") %>%
  html_nodes("h3 a") %>%
  html_attr("href")
  
url1 <- paste('https://www.bbsrc.ac.uk', url, sep='')
 
date <- read_html("bbsrc.html") %>%
  html_nodes("time") %>%
  html_text()

calltable <- data.frame(callname = call, deadline = date, call_url = paste('https://www.bbsrc.ac.uk', url, sep=''))
