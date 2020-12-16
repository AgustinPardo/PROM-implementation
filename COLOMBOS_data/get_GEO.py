from bs4 import BeautifulSoup # BeautifulSoup is in bs4 package
import requests

import re


entrada=open('colombos_mtube_exprdata_20151029.txt','r')
entrada =entrada.readlines()
experimentos=list(set(entrada[3].rstrip().split("\t")))

for experimento in experimentos[1:]:
        print("#####################################")
        print(experimento)
        URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+experimento
        content = requests.get(URL)

        soup = BeautifulSoup(content.text, 'html.parser')

        rows = soup.find_all('td')
        i=0
        for row in rows:          # Print all occurrences
            if row.get_text() == "Experiment type":
                print(row.get_text())
                next=i+1
                print(rows[next].get_text())

            if row.get_text() == "Summary":
                print(row.get_text())
                next=i+1
                print(rows[next].get_text())

            if row.get_text()[:9] == "Platforms":
                print("\t"+"Plataforma")
                print("\t"+"\t"+row.get_text())
                next=i+1
                print("\t"+"\t"+rows[next].get_text())

                plataformas= re.findall(r'GPL\w+', rows[next].get_text())
                for plataforma in plataformas:
                        URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+plataforma
                        content = requests.get(URL)
                        soup = BeautifulSoup(content.text, 'html.parser')
                        rows = soup.find_all('td')
                        x=0
                        for row in rows:
                                if row.get_text() == "Technology type":
                                        print("\t"+"\t"+row.get_text())
                                        sig=x+1
                                        print("\t"+"\t"+rows[sig].get_text())
                                x=x+1

            i=i+1

