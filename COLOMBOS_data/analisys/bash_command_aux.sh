cut -f -1,4-440 ../colombos_mtube_exprdata_20151029.txt > colombos_GSM27855.txt
sed 's/NaN/0/g' colombos_GSM27855.txt > colombos_GSM27855_noNaN.txt
