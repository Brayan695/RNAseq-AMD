library(dplyr)
meta = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/MetaSheet_1_4.csv")

feature_selected = meta %>% select(sample_id, mgs_level, A69S_GG, A69S_TT,
                                   Y402H_CC, Y402H_TT,
                                   mh_Dementia, mh_asthma,
                                   mh_dementia, mh_fibromyalgia,
                                   mh_high.chol, mh_hypothroidism,
                                   mh_smoker..1.pack.day..50.yr., mh_smoker..half.pk.day..50..yrs.,
                                   oc_cataract, oc_cataracts..OU., oc_confimred.pseudophakic..OU.,
                                   sex_bin)

write.csv(feature_selected, "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/MetaSheet_1_4_17.csv", row.names = F)
