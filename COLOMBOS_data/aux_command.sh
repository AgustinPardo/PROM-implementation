cut -f -1,4- colombos_mtube_exprdata_20151029.txt > colombos_all.txt
cat colombos_mtube_exprdata_20151029.txt| awk 'NR==6' > colombos_mtube_GPL.txt
