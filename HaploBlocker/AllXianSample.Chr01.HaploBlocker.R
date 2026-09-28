library(HaploBlocker)
source("HaploBlocker.skip_homo_check.R")
blocklist <- block_calculation(dhm ="./data/667XianSample.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.beagle.OnlyPhaseKeepMissing.Chr01.vcf.gz",adaptive_mode = TRUE,min_majorblock=1000,overlap_remove = FALSE)
saveRDS(blocklist, "667XianSample.HaploBlocker.Chr01.overlap.rds")
