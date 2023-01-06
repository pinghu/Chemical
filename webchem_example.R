library(webchem)
library(dplyr)
head(lc50)

lc50_sub <- lc50[1:15, ]

lc50_sub$inchikey <- cts_convert(lc50_sub$cas, from = "CAS", to = "InChIKey", match= "first", verbose = FALSE)
head(lc50_sub)

any(is.na(lc50_sub$inchikey))

x <- get_cid(lc50_sub$inchikey, from = "inchikey", match = "first", verbose = FALSE)
library(dplyr)
x1=data.frame(x)
lc50_sub$newinchikey=unlist(lc50_sub$inchikey)
lc50_sub2 <- full_join(lc50_sub, x1, by = c("newinchikey" = "query"))
head(lc50_sub2)

y <- pc_prop(lc50_sub2$cid, properties = c("IUPACName", "XLogP"))
#> https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/property/IUPACName,XLogP/JSON
y$CID <- as.character(y$CID)
lc50_sub3 <- full_join(lc50_sub2, y, by = c("cid" = "CID"))
head(lc50_sub3)

#####Not working, so 
out <- pan_query(lc50_sub3$cas, verbose = FALSE)
lc50_sub3$who_tox <- sapply(out, function(y) y$`WHO Acute Toxicity`)
lc50_sub3$common_name <- sapply(out, function(y) y$`Chemical name`)
# #equivalent with purrr package:
# lc50_sub3$who_tox <- map_chr(out, pluck, "WHO Acute Toxicity")
# lc50_sub3$common_name <- map_chr(out, pluck, "Chemical name")

#tidy up columns
lc50_done <- dplyr::select(lc50_sub3,  cas, inchikey, XLogP)
head(lc50_done)
lc50_done$CID <- cts_convert(lc50_done$cas, from = "CAS", to = "PubChem CID", match= "first", verbose = FALSE)

