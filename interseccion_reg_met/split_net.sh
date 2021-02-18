# Separo la red de sanz en archivos de regulator, y targets
rm regulator_filter.txt targets_filter.txt litevidence_filter.txt prob_prior_filter.txt;
# Estos sets son para comprar la red Sanz con el modelo iEK1011
filter13='Rv1638A Rv2737A Rv3022A Rv0590A Rv0470A Rv3724A Rv3224A Rv3724B Rv1000c Rv3324A Rv3224B Rv3197A Rv3312A';
filter41='Rv1909c Rv1990c Rv2358 Rv3164c Rv2669 Rv0353 Rv0001 Rv0182c Rv2710 Rv1657 Rv0485 Rv3286c Rv1931c Rv0117 Rv0823c Rv0348 Rv2034 Rv1404 Rv3414c Rv0981 Rv2720 Rv2745c Rv0212c Rv3849 Rv3223c Rv2069 Rv1343c Rv1267c Rv1379 Rv1956 Rv0967 Rv3133c Rv1221 Rv1846c Rv0465c Rv3080c Rv3648c Rv2711 Rv0586 Rv3291c Rv2359';
filter1030=$(<filter1030.dat);
filter="$filter13 $filter41 $filter1030";

while read p; 
        do IFS=, read -r a b <<< $p;
                if [[ $a == $b ]]
                then
                        :
                else
                        #if [[ $filter == *$a* ]] || [[ $filter == *$b* ]]
                        if [[ $filter == *$b* ]]
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





