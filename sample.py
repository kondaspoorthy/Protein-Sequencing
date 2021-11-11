l=[["blue",0.2,0.6],["yellow",0.4,0.3],["orange",0.6,0.5]]
k=["blule","yellow","pink","red"]
lst=[]
for each in k:
   if each in l[0:len(l)][0]:
      lst.append(True)
   else:
      lst.append(False)
print(lst)
