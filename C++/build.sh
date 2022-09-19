
FILE=$PWD
IND=$FILE/include
OBJ=$FILE/obj

g++ -Wall -O3 -I$IND $OBJ/urn_distribution.o $OBJ/survival.o $OBJ/vector.o $OBJ/kaplan_meier.o $OBJ/censoring_functions.o $OBJ/likelihood_functions.o $OBJ/brup_em.o BRUP_ER.cpp -o BRUP_ER
