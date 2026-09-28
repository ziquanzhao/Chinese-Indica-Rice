library(HaploBlocker)
source("HaploBlocker.skip_homo_check.R")
blocklist <- block_calculation(dhm ="./data/XianCultivar1.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.beagle.OnlyPhaseKeepMissing.Chr01.vcf.gz",adaptive_mode = TRUE,min_majorblock=1000,overlap_remove = FALSE)
saveRDS(blocklist, "XianCultivar1.HaploBlocker.Chr01.overlap.rds")
