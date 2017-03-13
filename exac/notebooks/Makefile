# Run all notebooks for ExAC analysis

all: ../reports/exac_data_tidying.html ../reports/exac_conservation_analysis.html \
	../reports/exac_exon_expression.html ../reports/exac_imperfects_analysis.html

../reports/exac_data_tidying.html: exac_data_tidying.Rmd ../produced_data/all_alignments.csv
	@ # produces report and RDS file ../produced_data/exac_data_clean.rds
	@ Rscript -e "rmarkdown::render('exac_data_tidying.Rmd', output_format = 'html_document')"
	@ mv exac_data_tidying.html ../reports/


../reports/exac_conservation_analysis.html: exac_conservation_analysis.Rmd ../produced_data/exac_data_clean.rds
	@ # produces report and RDS file ../produced_data/exac_data_cons_scores.rds
	@ Rscript -e "rmarkdown::render('exac_conservation_analysis.Rmd', output_format = 'html_document')"
	@ mv exac_conservation_analysis.html ../reports/


../reports/exac_exon_expression.html: exac_exon_expression.Rmd ../produced_data/exac_data_cons_scores.rds
	@ # produces reports and RDS file ../produced_data/exac_data_cons_exp.rds
	@ Rscript -e "rmarkdown::render('exac_exon_expression.Rmd', output_format = 'html_document')"
	@ mv exac_exon_expression.html ../reports/

../reports/exac_imperfects_analysis.html: exac_imperfects_analysis.Rmd ../produced_data/all_alignments.csv ../produced_data/exac_data_clean.rds
	@ Rscript -e "rmarkdown::render('exac_imperfects_analysis.Rmd', output_format = 'html_document')"
	@ mv exac_imperfects_analysis.html ../reports/