# 1.1 https://developers.google.com/earth-engine/datasets/catalog/MODIS_006_MCD12Q1
IGBPnames_006 <- c(
    "UNC", "ENF", "EBF", "DNF", "DBF", "MF", # 5 types
    "CSH", "OSH", "WSA", "SAV", "GRA", "WET",
    "CRO", "URB", "CNV", "SNO", "BSV", "water"
) # 1:17
IGBPcodes_006 <- setNames(0:17, IGBPnames_006)

# 1.2 https://developers.google.com/earth-engine/datasets/catalog/MODIS_051_MCD12Q1
IGBPnames_005 <- c(
    "water", # 0
    "ENF", "EBF", "DNF", "DBF", "MF",
    "CSH", "OSH", "WSA", "SAV", "GRA", "WET",
    "CRO", "URB", "CNV", "SNO", "BSV", "UNC"
) # 0:16, 254
IGBPcodes_005 <- setNames(c(0:16, 254), IGBPnames_005)

# IGBPnames5 <- c(
#     "water", "ENF", "EBF", "DNF", "DBF", "MF",
#     "CSH", "OSH", "WSA", "SAV", "GRA", "WET",
#     "CRO", "UB", "CNV", "SNO", "BAR", "UNC"
# )
