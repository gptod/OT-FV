import pandas
import sys
import os

file_location =sys.argv[1]

df = pandas.read_csv(file_location,skipinitialspace=True)
df.loc['total']= df.sum()

# iterating the columns
for col in df.columns:
    print(col)

print(df)
#df.iloc['total', df.columns.get_loc('np')] = ' '
df['error'] = df['error'].astype(str)
df.at['total','error']= ' '
print(df)


file_name = os.path.basename(file_location )  #eds_report.csv
location = os.path.abspath(os.path.dirname(file_location ))

outfile=location+'/man_'+file_name

df.to_csv(outfile, index=False)


#df1 = df.set_index('np')
#df1.loc['Total'] = df1.sum(numeric_only=True)
#print(df1)
#df.append(df.sum(numeric_only=True), ignore_index=True)



#print(df)
