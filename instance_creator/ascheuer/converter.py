import os, json, random, math

speed_profiles = [
	[1.5, 1.0, 1.67, 1.17, 1.33],  # fast
	[1.17, 0.67, 1.33, 0.83, 1.0], # normal
	[1.0, 0.33, 0.67, 0.5, 0.83]   # slow
]
zone_count = 5
cluster_count = 3

instance_files = ["rbg010a.tw", "rbg016a.tw", "rbg016b.tw", "rbg017.2.tw", "rbg017a.tw", "rbg017.tw", "rbg019a.tw", "rbg019b.tw", "rbg019c.tw", "rbg019d.tw", "rbg020a.tw", "rbg021.2.tw", "rbg021.3.tw", "rbg021.4.tw", "rbg021.5.tw", "rbg021.6.tw", "rbg021.7.tw", "rbg021.8.tw", "rbg021.9.tw", "rbg021.tw", "rbg027a.tw", "rbg031a.tw", "rbg033a.tw", "rbg034a.tw", "rbg035a.2.tw", "rbg035a.tw", "rbg038a.tw", "rbg040a.tw", "rbg041a.tw", "rbg042a.tw", "rbg048a.tw", "rbg049a.tw", "rbg050a.tw", "rbg050b.tw", "rbg050c.tw", "rbg055a.tw", "rbg067a.tw", "rbg086a.tw", "rbg092a.tw", "rbg125a.tw", "rbg132.2.tw", "rbg132.tw", "rbg152.3.tw", "rbg152.tw", "rbg172a.tw", "rbg193.2.tw", "rbg193.tw", "rbg201a.tw", "rbg233.2.tw", "rbg233.tw"]

# Read instances.
index = []

for file_name in instance_files:
	print(file_name)
	instance = {}
	with open(F"instances/{file_name}", "r") as ff:
		instance = ff.read()

	lines = instance.split("\n")
	n = int(lines[3].split(" ")[1])
	a = []
	b = []
	p = []
	d = [[0.0 for j in range(0, n)] for i in range(0, n)]
	for i in range(7, 7+n):
		p_i = float(lines[i].split(":")[1].split("[")[0].replace(" ", ""))
		a_i = float(lines[i].split("[")[1].split(",")[0].replace(" ", ""))
		b_i = float(lines[i].split("[")[1].split(",")[1].replace("]", ""))
		p.append(p_i)
		a.append(a_i)
		b.append(b_i)

	base = 7+n+1
	for i in range(base, base+n):
		line = [text for text in lines[i].replace("\t", " ").split(" ") if text != '']
		for j in range(0, n):
			d[i-base][j] = float(line[j]) + p[i-base]
	for i in range(0, n):
		d[i][i] = 0.0
		d[i][0] = 0.0
		d[n-1][i] = 0.0
	T = b[-1]
	a[0] = 0.0 # Tilk et al uses this approach.
	b[0] = 1000.0 # Tilk et al uses this approach.
	a[-1] = 0.0 # Tilk et al uses this approach.
	speed_zones = [[0.0, 0.2*T], [0.2*T,0.3*T], [0.3*T, 0.7*T], [0.7*T,0.8*T], [0.8*T, T]]

	I = {
		"digraph": {
			"vertex_count": n,
			"arc_count": (n-1)*(n-2),
			"arcs": [[1 if (i != j and j != 0 and i != n-1 and not(i == 0 and j == n)) else 0 for j in range(0, n)] for i in range(0, n)]
		},
		"cluster_count": cluster_count,
		"cluster_speeds": speed_profiles,
		"clusters": [[random.randrange(cluster_count) for j in range(0, n)] for i in range(0, n)],
		"speed_zone_count": zone_count,
		"speed_zones": speed_zones,
		"distances": d,
		"time_windows": [[a[i], b[i]] for i in range(0, n)],
		"start_depot": 0,
		"end_depot": n-1,
		"horizon": [0.0, T]
	}

	instance_name = file_name.replace(".tw", "")
	with open(F"output/{instance_name}.json", "w") as ff:
		ff.write(json.dumps(I))

	index.append({"instance_name": instance_name, "file_name": F"{instance_name}.json", "tags": [instance_name, F"N{n}"]})

# Create index.json
with open(F"output/index.json", "w") as ff:
		ff.write(json.dumps(index))

with open(F"output/solutions.json", "w") as ff:
		ff.write(json.dumps([]))