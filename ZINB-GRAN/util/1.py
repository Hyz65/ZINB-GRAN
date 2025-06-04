import pandas as pd

# Define the four sets of genes
genes_group1 = ["CEBPD", "NFE2", "MAF", "FOSL2", "JAZF1", "KLF10", "MAFB", "ID3", "GATA3", "RORA",
                "ATF3", "HES4", "EGR2", "EOMES", "MYBL1", "FOXP3", "EBF1", "JUN", "KLF4", "CEBPB",
                "BACH1", "HES1", "MEF2C", "JUND", "PAX5", "BCL11A", "MAX", "MYC", "BATF", "BACH2",
                "GATA2", "MYBL2", "ASCL2", "MYB", "EGR1", "FOSL1", "LEF1"]
genes_group2 = ["MYBL1", "ID3", "MAF", "GATA3", "FOXP3", "EBF1", "MYBL2", "BATF", "MAX", "JUN",
                "RORA", "BCL11A", "MYC", "JUND", "FOSL1", "EOMES", "LEF1", "MEF2C", "ATF3", "KLF10",
                "BACH2", "NFE2", "MYB", "PAX5", "HES4", "CEBPD", "FOSL2", "MAFB", "CEBPB", "ASCL2",
                "BACH1", "JAZF1", "HES1", "EGR2", "GATA2", "EGR1", "KLF4"]
genes_group3 = ["JUND", "BCL11A", "EGR1", "LEF1", "FOSL2", "BACH2", "JAZF1", "BATF", "FOSL1", "ASCL2",
                "FOXP3", "MYBL2", "JUN", "MEF2C", "GATA2", "MAX", "KLF10", "MAFB", "PAX5", "ID3",
                "CEBPD", "HES1", "MYC", "MYBL1", "RORA", "HES4", "EGR2", "BACH1", "NFE2", "EBF1",
                "ATF3", "KLF4", "MYB", "GATA3", "MAF", "EOMES", "CEBPB"]
genes_group4 = ["FOSL2", "BATF", "LEF1", "MYBL2", "BACH2", "FOSL1", "MAFB", "GATA2", "HES4", "GATA3",
                "ATF3", "EOMES", "NFE2", "FOXP3", "MYC", "HES1", "EGR2", "JUN", "JAZF1", "MYB",
                "MAF", "KLF4", "EGR1", "EBF1", "KLF10", "BCL11A", "MEF2C", "RORA", "BACH1", "MAX",
                "CEBPB", "ASCL2", "PAX5", "JUND", "MYBL1", "CEBPD", "ID3"]

# Convert to sets for easier comparison
set1 = set(genes_group1)
set2 = set(genes_group2)
set3 = set(genes_group3)
set4 = set(genes_group4)

# Find common genes across all four groups
common_genes = set1.intersection(set2).intersection(set3).intersection(set4)

print(common_genes)
