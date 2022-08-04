
import numpy as np
my_file = open("./frame_change.txt", "r")
content = my_file.read()
content_list = content.split("\n")
for i in range(1,len(content_list)-1):
       val = int(content_list[i]) - int(content_list[i-1])
       if val > 200:
            print(val)

print(''''''''''\n''''''''')
import numpy as np
my_file = open("./id_hydronium.dat", "r")
content = my_file.read()
content_list = content.split("\n")
value = []
for i in content_list[0:-1]:
    append_me = i.split()[1]
    value.append(append_me)
    #value.append(i.split()[0])

old_count = 0
for i in range(1,len(value)):
    if(value[i] != value[i-1]):
        count = i
        val = count - old_count 
        if val > 200:
            print(val)
        old_count = count
