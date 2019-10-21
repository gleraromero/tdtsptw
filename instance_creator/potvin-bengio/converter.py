import os, json, random, math

speed_profiles = [
	[1.5, 1.0, 1.67, 1.17, 1.33],  # fast
	[1.17, 0.67, 1.33, 0.83, 1.0], # normal
	[1.0, 0.33, 0.67, 0.5, 0.83]   # slow
]
zone_count = 5
cluster_count = 3

instance_files = ["rc_201.1.txt", "rc_202.4.txt", "rc_204.3.txt", "rc_206.2.txt", "rc_208.1.txt", "rc_201.2.txt", "rc_203.1.txt", "rc_204.4.txt", "rc_206.3.txt", "rc_208.2.txt", "rc_201.3.txt", "rc_203.2.txt", "rc_205.1.txt", "rc_206.4.txt", "rc_208.3.txt", "rc_201.4.txt", "rc_203.3.txt", "rc_205.2.txt", "rc_207.1.txt", "rc_202.1.txt", "rc_203.4.txt", "rc_205.3.txt", "rc_207.2.txt", "rc_202.2.txt", "rc_204.1.txt", "rc_205.4.txt", "rc_207.3.txt", "rc_202.3.txt", "rc_204.2.txt", "rc_206.1.txt", "rc_207.4.txt"]

def dist(x1, x2, y1, y2):
	return math.sqrt((x1-x2)**2 + (y1-y2)**2)

# Read instances.
index = []

for file_name in instance_files:
	print(file_name)
	instance = {}
	with open(F"instances/{file_name}", "r") as ff:
		instance = ff.read()

	lines = instance.split("\n")
	a = []
	b = []
	x = []
	y = []
	p = []
	n = int(lines[0])
	for i in range(1, n+1):
		line = [text for text in lines[i].replace("\t", " ").split(" ") if text != '']
		x.append(float(line[1]))
		y.append(float(line[2]))
		a.append(float(line[4]))
		b.append(float(line[5]))
		p.append(float(line[6]))
	x.append(x[0])
	y.append(y[0])
	a.append(a[0])
	b.append(b[0])
	p.append(p[0])

	T = b[0] * 100.0
	speed_zones = [[0.0, 0.2*T], [0.2*T,0.3*T], [0.3*T, 0.7*T], [0.7*T,0.8*T], [0.8*T, T]]
	T = b[0]
	n = len(x)
	I = {
		"digraph": {
			"vertex_count": n,
			"arc_count": (n-1)*(n-2),
			"arcs": [[1 if (i != j and j != 0 and i != n-1 and not(i == 0 and j == n)) else 0 for j in range(0, n)] for i in range(0, n)],
			"coordinates": [[x[i], y[i]] for i in range(0, n)]
		},
		"cluster_count": cluster_count,
		"cluster_speeds": speed_profiles,
		"clusters": [[random.randrange(cluster_count) for j in range(0, n)] for i in range(0, n)],
		"speed_zone_count": zone_count,
		"speed_zones": speed_zones,
		"distances": [[dist(x[i], x[j], y[i], y[j]) for j in range(0, n)] for i in range(0, n)],
		"time_windows": [[a[i], b[i]] for i in range(0, n)],
		"start_depot": 0,
		"end_depot": n-1,
		"horizon": [0.0, T]
	}
	# Add process times to distances.
	for i in range(0, n):
		for j in range(0, n):
			I["distances"][i][j] += p[i]

	# Multiply all values by 100.00
	for i in range(0, n):
		I["time_windows"][i][0] = I["time_windows"][i][0] * 100.0
		I["time_windows"][i][1] = I["time_windows"][i][1] * 100.0
		I["horizon"][0] = 0.0
		I["horizon"][1] = b[0] * 100.0
		for j in range(0, n):
			I["distances"][i][j] = I["distances"][i][j] * 100.0

	# Truncate distances and modify them to satisfy triangle inequality.
	for i in range(0, n):
		for j in range(0, n):
			I["distances"][i][j] = math.floor(I["distances"][i][j])
	# for i in range(0, n):
	# 	for j in range(0, n):
	# 		for k in range(0, n):
	# 			if I["distances"][i][j] > I["distances"][i][k] + I["distances"][k][j]:
	# 				I["distances"][i][j] = I["distances"][i][k] + I["distances"][k][j]

	instance_name = file_name.replace(".txt", "")
	with open(F"output/{instance_name}.json", "w") as ff:
		ff.write(json.dumps(I))

	index.append({"instance_name": instance_name, "file_name": F"{instance_name}.json", "tags": [instance_name, F"N{n}"], "n": n})

def sort_key(item):
	return item["instance_name"]

index.sort(key = sort_key)

for item in index:
	del item["n"]

# Create index.json
with open(F"output/index.json", "w") as ff:
		ff.write(json.dumps(index))

with open(F"output/solutions.json", "w") as ff:
		ff.write(json.dumps([]))