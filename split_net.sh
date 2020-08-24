# Separo la red de sanz en archivos de regulator, y targets
rm regulator_filter.txt targets_filter.txt litevidence_filter.txt prob_prior_filter.txt;
filter1='Rv1638A Rv2737A Rv3022A Rv0590A Rv0470A Rv3724A Rv3224A Rv3724B Rv1000c Rv3324A Rv3224B Rv3197A Rv3312A';
filter2=$(<filter.dat);
filter="$filter1 $filter2";

while read p; 
        do IFS=, read -r a b <<< $p;
                if [[ $a == $b ]]
                then
                        :
                else
                        if [[ $filter == *$a* ]] | [[ $filter == *$b* ]]
                        then
                               :
                        else
                                echo $a >> regulator_filter.txt;
                                echo $b >> targets_filter.txt;
                                echo 0  >> litevidence_filter.txt;
                                echo 0  >> prob_prior_filter.txt;
                        fi
                fi
done < sanz_tb_net.txt;





