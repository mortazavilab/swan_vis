set -e
# cd ~/mortazavi_lab/bin/TALON/

rm ~/mortazavi_lab/bin/refactor_splice_graph/testing/combine.db

talon_initialize_database \
	--f ~/mortazavi_lab/bin/refactor_splice_graph/testing/test_combine_2.gtf \
	--g hg38 \
	--a combine \
	--l 0 \
	--o ~/mortazavi_lab/bin/refactor_splice_graph/testing/combine

cd ~/mortazavi_lab/bin/refactor_splice_graph/

python talon_db_to_graph.py \
	-db testing/combine.db \
	--o testing/combine_db 

python plot_graph.py \
	-g testing/test_combine_splicing_graph.p \
	--o ~/mortazavi_lab/bin/refactor_splice_graph/testing/figures/test_combine_2 \
	--color_introns \
	--color_alt_TSS \
	--color_alt_TES

open testing/figures/test_combine_2*
