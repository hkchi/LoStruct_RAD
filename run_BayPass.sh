## Genome-environment association analysis using BayPass
## Kaichi Huang 2019 Sept

# Generate input file from VCF
zcat GSD_RAD.HA412.good_snps2.vcf.gz | perl vcf2baypass.pl pop_info.tsv GSD_RAD.HA412.good_snps2.baypass

# Run BayPass under the core model mode to generate covariance matrix
i_baypass -npop 20 -gfile GSD_RAD.HA412.good_snps2.baypass.txt -outprefix all -nthreads 20

# Run BayPass under the standard covariate model using importance sampling (IS) estimator
i_baypass -npop 20 -gfile GSD_RAD.HA412.good_snps2.baypass.txt -efile soil.baypass.txt -scalecov -omegafile all_mat_omega.out -outprefix all_soilGEA -nthreads 20

# Produce POD samples with 1,000 SNPs and run BayPass
Rscript pod.R
i_baypass -npop 20 -gfile G.GSD_RAD.HA412.good_snps2.baypass.pod -efile soil.baypass.txt -scalecov -omegafile all_mat_omega.out -outprefix pod_soilGEA -nthreads 20
