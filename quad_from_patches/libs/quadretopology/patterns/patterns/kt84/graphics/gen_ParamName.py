import sys

if len(sys.argv) != 2:
    print "Usage: gen_ParamName.py ParamNameX.Y.txt"
    sys.exit()

fin = open(sys.argv[1])

lines = fin.readlines()

for i in range(len(lines) / 3):
    name = lines[3 * i].rstrip("\n")
    desc = lines[3 * i + 1].rstrip("\n")
    print "// " + desc
    print name[3:] + " = " + name + ","
