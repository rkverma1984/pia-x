filtered_keyword = 'Only H5'
filtered_keyword2 = 'Mini'
filtered_keyword3 = 'mini'

f = open('extracted_information_from_excel.txt','r')
lines=f.read().split('\n')
f.close()
lines = [line for line in lines if line != '']

#H5 fragment note is in column 15 (1-based)
f = open('filtered_information.txt','w')
for line in lines:
    if not filtered_keyword in line and not filtered_keyword2 in line and not filtered_keyword3 in line:
        f.write(line+'\n')
f.close()

