## -----------------------------------------------------------------------------
library(roadoi)
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                             "10.1103/physreve.88.012814"), 
                    email = "najko.jahn@gmail.com")

## -----------------------------------------------------------------------------
library(dplyr)
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                             "10.1103/physreve.88.012814",
                             "10.1093/reseval/rvaa038",
                             "10.1101/2020.05.22.111294",
                             "10.1093/bioinformatics/btw541"), 
                    email = "najko.jahn@gmail.com", .flatten = TRUE) %>%
  dplyr::count(is_oa, evidence, is_best) 

## -----------------------------------------------------------------------------
roadoi::oadoi_fetch(dois = c("10.1186/s12864-016-2566-9",
                             "10.1103/physreve.88.012814"), 
                    email = "najko.jahn@gmail.com", 
                    .progress = "text")

