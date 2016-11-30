from struct import *

f = open("test_9.22","rb")
try:
    count = 0
    count1 = 0
    byte = f.read(4)
    while byte != "":
        count1 = count1 +2
        if count < 2:
            print unpack('>i',byte)
            print unpack('>i',byte)
        else:
            print unpack('>d',byte)
            print unpack('>d',byte)
        if count < 1:
            byte = f.read(4)
        else:
            byte = f.read(8)
        count = count+1
    
finally:
    f.close()



print count1
print count1

