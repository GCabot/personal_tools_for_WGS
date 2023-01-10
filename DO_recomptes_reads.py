#!/usr/bin/python
# Writen by Gabriel Cabot, Institut d'Investigació Sanitària Illes Balears (IdISBa)
# email: gabriel.cabot@ssib.es
# Antibiotic Resistance and Pathogenicity of Bacterial Infections Group (https://arpbigidisba.com/)
# GitHub (https://github.com/GCabot)

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pylab as plt
import os
import time
from timeit import default_timer as timer
from os import system
import csv
import math
import seaborn as sns
plt.style.use("seaborn")

###############################################################################
#   1. DEFINING VARIABLES
###############################################################################

numerical="0123456789"
alphabetical="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
alphanumerical=alphabetical+numerical
indelidentifier="+-"
equalasref=",."
trashyones="^$!"
anychar=alphabetical+numerical+indelidentifier+equalasref+trashyones

bases="AaTtCcGg"

###############################################################################
#   2. OPENING DATA_FILE
###############################################################################

system("clear")

print ("Writen by Gabriel Cabot, Institut d'Investigació Sanitària Illes Balears (IdISBa)")
print ()

print ("\x1B[3m" + "This script analyze the relative frequency of each event for a given totalpileup file" + "\x1B[0m")

print ("\x1B[3m" + "Please, use it UNDER YOUR OWN responsability" + "\x1B[0m")

print("{:>60}".format("G. Cabot"))

time.sleep(5)

system("clear")

print ()


datafile=input("Select a totalpileup file to be analysed: ")

outfolder=input("Select an output folder to store analysis: ")

outfile=(outfolder + "my_processed_table.csv")

print (f"File to analyse is: {datafile}")
print ()
print (f"Folder to store analysis is: {outfolder}")

confirmation=input("Is this correct?")

if confirmation == 'yes' :

    pass

if confirmation == 'no':

    quit()

    # Opening & parsing datafile:
dataset = pd.read_csv(datafile, sep='\t', header=None)
                      
dataset.columns=("NODE", "POSITION", "REFERENCE", "CHANGE", "Q1", "Q2", "Q3",
                "READ_N", "READS", "QUAL")

#print("This is the DATASET")
#print(dataset)
#print()

positions=dataset.POSITION.sort_values().unique()
events=("REFERENCE", "POSITION", "NO_CHANGE", "A/a", "T/t", "C/c", "G/g", "Ins", "Del", "READ_N", "READS") # Between brakets to avoid an order change
eventsplot=("NO_CHANGE", "A/a", "T/t", "C/c", "G/g", "Ins", "Del")
counttable=pd.DataFrame(index=positions,columns=events)
plottable=pd.DataFrame(index=positions,columns=eventsplot)

#print("This is your empty table")
#print(counttable)   # Empty table to fill with character count.

#print ()
#print (positions)


###############################################################################
#   3. 
###############################################################################

for i in range (len(dataset)):
    position=dataset.iloc[i].POSITION
    reference=dataset.iloc[i].REFERENCE
    readcount=dataset.iloc[i].READ_N
    sequence=dataset.iloc[i].READS
    print ("Position:", position, reference, readcount, sequence, len(sequence))

    z=0
    NO_CHANGE=0
    A=0
    T=0
    G=0
    C=0
    Ins=0
    Del=0

    while z in range(len(sequence)):
        
        if z < (len(sequence)-1):

            charact=sequence[z]
            checindel=sequence[z+1]

        if z == (len(sequence)-1):

            charact=sequence[z]
            checindel=sequence[z]

        if charact == 'A':
            if checindel not in indelidentifier:

                A+=1
        
        if charact == 'a':
            if checindel not in indelidentifier:
                
                A+=1

        if charact == 'T':
            if checindel not in indelidentifier:

                T+=1
        
        if charact == 't':
            if checindel not in indelidentifier:

                T+=1

        if charact == 'G':
            if checindel not in indelidentifier:

                G+=1
        
        if charact == 'g':
            if checindel not in indelidentifier:

                G+=1

        if charact == 'C':
            if checindel not in indelidentifier:

                C+=1
        
        if charact == 'c':
            if checindel not in indelidentifier:
                
                C+=1

        if charact in equalasref:
            if checindel not in indelidentifier:

                NO_CHANGE+=1

        if charact in trashyones:
            if checindel not in indelidentifier:
                
                pass

        if charact == "+":
            print("Into the Insertion loop")

            inschar1=sequence[z+1]
            
            if inschar1 not in numerical:

                pass

            if inschar1 in numerical:
                
                inschar2=sequence[z+2]

                print (f"The first insertion character is {inschar1}")

                if inschar2 not in numerical:

                    lenght=int(inschar1)
                    Ins+=1
                    z+=int(lenght+1)

                if inschar2 in numerical:
                    print (f"The second insertion character is {inschar2}")

                    lenght=(int(inschar1)*10+int(inschar2))
                    Ins+=1
                    z+=int(lenght+2)

            print (f"Insertion lenght is {lenght}")

        if charact == "-":
            print("Into the Deletion loop")

            delchar1=sequence[z+1]

            if delchar1 not in numerical:

                pass

            if delchar1 in numerical:
                print (f"The first deletion character is {delchar1}")
                
                delchar2=sequence[z+2]

                if delchar2 not in numerical:

                    lenght=int(delchar1)
                    Del+=1
                    z+=int(lenght+1)

                if delchar2 in numerical:
                    print (f"The second deletion character is {delchar2}")
                    lenght=(int(delchar1)*10+int(delchar2))
                    Del+=1
                    z+=int(lenght+2)
                                   
            print (f"Deletion lenght is {lenght}")

        z+=1
    
    #print (f"El número de lecturas coincidentes con la referencia es {NO_CHANGE}")
    Equalperc=round(NO_CHANGE/readcount, 3)
    print (f"Percentage of reads equal to reference is {Equalperc}")
    print ()
    #print (f"El número de lecturas correspondientes a INSerciones es {Ins}")
    Insperc=round(Ins/readcount, 3)
    print (f"Percentage of reads with an Insertion is {Insperc}")
    #print (f"El número de lecturas correspondientes a DELeciones es {Del}")
    Delperc=round(Del/readcount, 3)
    print (f"Percentage of reads with a Deletion is {Delperc}")
    #print (f"El número de lecturas correspondientes a ADENINA es {A}")
    Aperc=round(A/readcount, 3)
    print (f"Percentage of reads of Adenine is {Aperc}")
    #print (f"El número de lecturas correspondientes a TIMINA es {T}")
    Tperc=round(T/readcount, 3)
    print (f"Percentage of reads of Timine is {Tperc}")
    #print (f"El número de lecturas correspondientes a CITOSINA es {C}")
    Cperc=round(C/readcount, 3)
    print (f"Percentage of reads of Citosine is {Cperc}")
    #print (f"El número de lecturas correspondientes a GUANINA es {G}")
    Gperc=round(G/readcount, 3)
    print (f"Percentage of reads of Guanine is {Gperc}")
    print ()

    counttable.loc[[int(position)],['POSITION']]=position
    counttable.loc[[int(position)],['REFERENCE']]=reference
    counttable.loc[[int(position)],['NO_CHANGE']]=(Equalperc*100)
    counttable.loc[[int(position)],['A/a']]=(Aperc*100)
    counttable.loc[[int(position)],['T/t']]=(Tperc*100)
    counttable.loc[[int(position)],['G/g']]=(Gperc*100)
    counttable.loc[[int(position)],['C/c']]=(Cperc*100)
    counttable.loc[[int(position)],['Ins']]=(Insperc*100)
    counttable.loc[[int(position)],['Del']]=(Delperc*100)
    counttable.loc[[int(position)],['READ_N']]=readcount
    counttable.loc[[int(position)],['READS']]=sequence

print(counttable)



counttable.to_csv('/home/micro128g2/Documentos/COUNTREAD_SCRIPT/my_processed_table.csv', sep='\t', header=True, index=True)

counttable.to_excel('/home/micro128g2/Documentos/COUNTREAD_SCRIPT/my_processed_table.xlsx', sheet_name='Hoja1')