import csv
import struct

from collections import defaultdict

columns = defaultdict(list) # each value in each column is appended to a list

with open('sigma_xx.csv') as f:
    reader = csv.DictReader(f) # read rows into a dictionary format
    for row in reader: # read a row as {column1: value1, column2: value2,...}
        for (k,v) in row.items(): # go over each column name and value 
            columns[k].append(v) # append the value into the appropriate list
                                 # based on column name k

SigmaXX = columns["SigmaXX"]
X = columns["Points:0"]

file = open("sigma_xx.txt","w")
file.write("data {")
file.write("\n x, y\n")
for i in range(0,1001):
	file.write("%s, %s\n" % (X[i] ,SigmaXX[i]))
file.write("\n };\n")	
file.close() 
