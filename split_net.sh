# Separo la red de sanz en archivos de TFs, y targets
rm TF.txt targets.txt;
while read p; 
        do IFS=, read -r a b <<< $p;
                if [[ $a == $b ]]
                then
                        :
                else
                        echo $a >> TF.txt;
                        echo $b >> targets.txt;
                        echo 0  >> litevidence.txt;
                        echo 0  >> prob_prior.txt;
                fi
done < sanz_tb_net.txt;



