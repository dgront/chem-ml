import sys

data = open(sys.argv[1])
header = data.readline()
print(header)
uniq = set()
for line in data:
    tokens = line.strip().split()
    if tokens[2] not in uniq:
        print(line.strip())
        uniq.add(tokens[2])
        
    
    
