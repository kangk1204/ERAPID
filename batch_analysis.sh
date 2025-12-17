python erapid.py --gse GSE95587
cp GSE95587_coldata_for_edit.tsv GSE95587
python erapid.py --gse GSE95587 --phase analysis --group_ref control --evidence_keywords alzheimer

python erapid.py --gse GSE80655
cp GSE80655_coldata_for_edit.tsv GSE80655
python erapid.py --gse GSE80655 --phase analysis --group_ref control,bp,md --evidence_keywords schizophrenia

python erapid.py --gse GSE104704
cp GSE104704_coldata_for_edit.tsv GSE104704
python erapid.py --gse GSE104704 --phase analysis --group_ref young,old --evidence_keywords alzheimer

python erapid.py --gse GSE125583
cp GSE125583_coldata_for_edit.tsv GSE125583
python erapid.py --gse GSE125583 --phase analysis --group_ref control --evidence_keywords alzheimer

python erapid.py --gse GSE153873
cp GSE153873_coldata_for_edit.tsv GSE153873
python erapid.py --gse GSE153873 --phase analysis --group_ref young,old --evidence_keywords alzheimer

python erapid.py --gse GSE190185
cp GSE190185_coldata_for_edit.tsv GSE190185
python erapid.py --gse GSE190185 --phase analysis --group_ref 3_3,4_4 --evidence_keywords alzheimer,apoe4

python erapid.py --phase meta --gse GSE104704,GSE125583,GSE153873,GSE95587 --method both --out meta_results/ad_control_only --evidence_keywords alzheimer --evidence_top_n 300 --meta_include GSE153873:ad_vs_old_common_deg
