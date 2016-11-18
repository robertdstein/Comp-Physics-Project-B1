import csv

source = "lifetime.txt"

datalist = []

with open(source, 'rb') as f:
	reader = csv.reader(f, delimiter=' ', quotechar='|')
	for row in reader:
		vals = [float(row[0]), float(row[1])]
		datalist.append(vals)

data = np.array(datalist)

def run():
	return data
