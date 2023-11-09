import sys

filename = sys.argv[1] # file ".txt" storing the dataset (not-preprocessed)
storeIdMap = sys.argv[2] # 0 to not store IDs; 1: store map of IDs
fout = filename.strip().split(".txt")[0]
fout = fout +"-remapped.txt"

# Load edges in custom format for node IDs
fileread = open(filename, 'r')
edgesNotSorted = []
for line in fileread:
    chunks = line.strip().split(" ")
    if len(chunks) != 3:
        print("Format is not src dst tim!!")
        exit(0)
    else:
        edgesNotSorted.append((chunks[0], chunks[1], int(chunks[2])))
fileread.close()

if len(edgesNotSorted) < 1:
    print("No edges!!")
    exit(0)

#print(edgesNotSorted)
# sort and remap node IDs also make timestamps start from 1
edgesSorted = sorted(edgesNotSorted, key=lambda edge: edge[2], reverse=False)
#print(edgesSorted)
MINTIM = 1
shift = 0
if edgesSorted[0][2] < MINTIM:
    shift = MINTIM - edgesSorted[0][2]
print("Timestamps will be shifted (in time) by", shift, "unit steps")

mapID = {}
ID = 0
remappedEdges = []

for edge in edgesSorted:
    src = edge[0]
    dst = edge[1]
    tim = edge[2]
    if src == dst: # Self-loops are not admitted!
        continue
    if src not in mapID.keys():
        mapID[src] = ID
        ID = ID + 1
    if dst not in mapID.keys():
        mapID[dst] = ID
        ID = ID + 1
    remappedEdges.append([mapID[src], mapID[dst], tim+shift])
#print(remappedEdges)
outdata = open(fout, 'w')
for edge in remappedEdges:
    outdata.write(str(edge[0])+ " " + str(edge[1]) + " " + str(edge[2]) + "\n")
outdata.close()

if storeIdMap == "1":
    fout = filename.strip().split(".txt")[0]
    fout = fout +"-IDMAP.txt"
    outdata = open(fout, 'w')
    outdata.write("Original-ID New-ID\n")
    for key in sorted(mapID.keys()):
        outdata.write(str(key) + " " + str(mapID[key]) + "\n")
    outdata.close()
