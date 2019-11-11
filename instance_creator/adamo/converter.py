import os, json, random, math

speed_zone_count = 73
cluster_count = 3

# Load Jams.
jam_files = ["Jams70A", "Jams70B", "Jams80A", "Jams80B", "Jams90A", "Jams90B"]
jam_names = ["70A", "70B", "80A", "80B", "90A", "90B"]
jam_types = ["A", "B", "A", "B", "A", "B"]
jam_numbers = ["70", "70", "80", "80", "90", "90"]
jams = []
for jam_file in jam_files:
	with open(F"instances/Pattern/{jam_file}", "r") as ff:
		jams.append([[float(x) for x in line.split("\t")] for line in ff.read().split("\n") if line != ""])

# Read instances.
instance_files = ["15A_10_new.txt", "15A_1_new.txt", "15A_2_new.txt", "15A_3_new.txt", "15A_4_new.txt", "15A_5_new.txt", "15A_6_new.txt", "15A_7_new.txt", "15A_8_new.txt", "15A_9_new.txt", "15B_10_new.txt", "15B_1_new.txt", "15B_2_new.txt", "15B_3_new.txt", "15B_4_new.txt", "15B_5_new.txt", "15B_6_new.txt", "15B_7_new.txt", "15B_8_new.txt", "15B_9_new.txt", "15C_10_new.txt", "15C_1_new.txt", "15C_2_new.txt", "15C_3_new.txt", "15C_4_new.txt", "15C_5_new.txt", "15C_6_new.txt", "15C_7_new.txt", "15C_8_new.txt", "15C_9_new.txt", "20A_10_new.txt", "20A_1_new.txt", "20A_2_new.txt", "20A_3_new.txt", "20A_4_new.txt", "20A_5_new.txt", "20A_6_new.txt", "20A_7_new.txt", "20A_8_new.txt", "20A_9_new.txt", "20B_10_new.txt", "20B_1_new.txt", "20B_2_new.txt", "20B_3_new.txt", "20B_4_new.txt", "20B_5_new.txt", "20B_6_new.txt", "20B_7_new.txt", "20B_8_new.txt", "20B_9_new.txt", "20C_10_new.txt", "20C_1_new.txt", "20C_2_new.txt", "20C_3_new.txt", "20C_4_new.txt", "20C_5_new.txt", "20C_6_new.txt", "20C_7_new.txt", "20C_8_new.txt", "20C_9_new.txt", "25A_10_new.txt", "25A_1_new.txt", "25A_2_new.txt", "25A_3_new.txt", "25A_4_new.txt", "25A_5_new.txt", "25A_6_new.txt", "25A_7_new.txt", "25A_8_new.txt", "25A_9_new.txt", "25B_10_new.txt", "25B_1_new.txt", "25B_2_new.txt", "25B_3_new.txt", "25B_4_new.txt", "25B_5_new.txt", "25B_6_new.txt", "25B_7_new.txt", "25B_8_new.txt", "25B_9_new.txt", "25C_10_new.txt", "25C_1_new.txt", "25C_2_new.txt", "25C_3_new.txt", "25C_4_new.txt", "25C_5_new.txt", "25C_6_new.txt", "25C_7_new.txt", "25C_8_new.txt", "25C_9_new.txt", "30A_10_new.txt", "30A_1_new.txt", "30A_2_new.txt", "30A_3_new.txt", "30A_4_new.txt", "30A_5_new.txt", "30A_6_new.txt", "30A_7_new.txt", "30A_8_new.txt", "30A_9_new.txt", "30B_10_new.txt", "30B_1_new.txt", "30B_2_new.txt", "30B_3_new.txt", "30B_4_new.txt", "30B_5_new.txt", "30B_6_new.txt", "30B_7_new.txt", "30B_8_new.txt", "30B_9_new.txt", "30C_10_new.txt", "30C_1_new.txt", "30C_2_new.txt", "30C_3_new.txt", "30C_4_new.txt", "30C_5_new.txt", "30C_6_new.txt", "30C_7_new.txt", "30C_8_new.txt", "30C_9_new.txt", "35A_10_new.txt", "35A_1_new.txt", "35A_2_new.txt", "35A_3_new.txt", "35A_4_new.txt", "35A_5_new.txt", "35A_6_new.txt", "35A_7_new.txt", "35A_8_new.txt", "35A_9_new.txt", "35B_10_new.txt", "35B_1_new.txt", "35B_2_new.txt", "35B_3_new.txt", "35B_4_new.txt", "35B_5_new.txt", "35B_6_new.txt", "35B_7_new.txt", "35B_8_new.txt", "35B_9_new.txt", "35C_10_new.txt", "35C_1_new.txt", "35C_2_new.txt", "35C_3_new.txt", "35C_4_new.txt", "35C_5_new.txt", "35C_6_new.txt", "35C_7_new.txt", "35C_8_new.txt", "35C_9_new.txt", "40A_10_new.txt", "40A_1_new.txt", "40A_2_new.txt", "40A_3_new.txt", "40A_4_new.txt", "40A_5_new.txt", "40A_6_new.txt", "40A_7_new.txt", "40A_8_new.txt", "40A_9_new.txt", "40B_10_new.txt", "40B_1_new.txt", "40B_2_new.txt", "40B_3_new.txt", "40B_4_new.txt", "40B_5_new.txt", "40B_6_new.txt", "40B_7_new.txt", "40B_8_new.txt", "40B_9_new.txt", "40C_10_new.txt", "40C_1_new.txt", "40C_2_new.txt", "40C_3_new.txt", "40C_4_new.txt", "40C_5_new.txt", "40C_6_new.txt", "40C_7_new.txt", "40C_8_new.txt", "40C_9_new.txt", "45A_10_new.txt", "45A_1_new.txt", "45A_2_new.txt", "45A_3_new.txt", "45A_4_new.txt", "45A_5_new.txt", "45A_6_new.txt", "45A_7_new.txt", "45A_8_new.txt", "45A_9_new.txt", "45B_10_new.txt", "45B_1_new.txt", "45B_2_new.txt", "45B_3_new.txt", "45B_4_new.txt", "45B_5_new.txt", "45B_6_new.txt", "45B_7_new.txt", "45B_8_new.txt", "45B_9_new.txt", "45C_10_new.txt", "45C_1_new.txt", "45C_2_new.txt", "45C_3_new.txt", "45C_4_new.txt", "45C_5_new.txt", "45C_6_new.txt", "45C_7_new.txt", "45C_8_new.txt", "45C_9_new.txt", "50A_10_new.txt", "50A_1_new.txt", "50A_2_new.txt", "50A_3_new.txt", "50A_4_new.txt", "50A_5_new.txt", "50A_6_new.txt", "50A_7_new.txt", "50A_8_new.txt", "50A_9_new.txt", "50B_1_new.txt", "50B_2_new.txt", "50B_3_new.txt", "50B_4_new.txt", "50B_5_new.txt", "50B_6_new.txt", "50B_7_new.txt", "50B_8_new.txt", "50B_9_new.txt", "50B__new.txt", "50C_1_new.txt", "50C_2_new.txt", "50C_3_new.txt", "50C_4_new.txt", "50C_5_new.txt", "50C_6_new.txt", "50C_7_new.txt", "50C_8_new.txt", "50C_9_new.txt", "50C__new.txt", "55A_10_new.txt", "55A_1_new.txt", "55A_2_new.txt", "55A_3_new.txt", "55A_4_new.txt", "55A_5_new.txt", "55A_6_new.txt", "55A_7_new.txt", "55A_8_new.txt", "55A_9_new.txt", "55B_10_new.txt", "55B_1_new.txt", "55B_2_new.txt", "55B_3_new.txt", "55B_4_new.txt", "55B_5_new.txt", "55B_6_new.txt", "55B_7_new.txt", "55B_8_new.txt", "55B_9_new.txt", "55C_10_new.txt", "55C_1_new.txt", "55C_2_new.txt", "55C_3_new.txt", "55C_4_new.txt", "55C_5_new.txt", "55C_6_new.txt", "55C_7_new.txt", "55C_8_new.txt", "55C_9_new.txt", "60A_10_new.txt", "60A_1_new.txt", "60A_2_new.txt", "60A_3_new.txt", "60A_4_new.txt", "60A_5_new.txt", "60A_6_new.txt", "60A_7_new.txt", "60A_8_new.txt", "60A_9_new.txt", "60B_10_new.txt", "60B_1_new.txt", "60B_2_new.txt", "60B_3_new.txt", "60B_4_new.txt", "60B_5_new.txt", "60B_6_new.txt", "60B_7_new.txt", "60B_8_new.txt", "60B_9_new.txt", "60C_10_new.txt", "60C_1_new.txt", "60C_2_new.txt", "60C_3_new.txt", "60C_4_new.txt", "60C_5_new.txt", "60C_6_new.txt", "60C_7_new.txt", "60C_8_new.txt", "60C_9_new.txt"]
index = []

for file_name in instance_files:
	instance = {}
	with open(F"instances/{file_name}", "r") as ff:
		instance = ff.read()
	print(file_name)

	lines = instance.split("\n")
		
	# Parse instance file.
	n = int(lines[0])
	d = [[float(x) for x in lines[i].strip().split(" ")] for i in range(1, n+1)]
	a = [float(lines[i].strip().split("\t")[0]) for i in range(n+1, n+1+n)]
	b = [float(lines[i].strip().split("\t")[1]) for i in range(n+1, n+1+n)]
	C_START = n+1+n+1
	C = [[int(x)-1 for x in lines[i].strip().split(" ")] for i in range(C_START, C_START+n+1)]
	T_START = C_START+n+1+1
	T = [[float(x) for x in lines[i].strip().split("\t")] for i in range(T_START, T_START + speed_zone_count)]
	S_START = T_START + speed_zone_count + 1
	s = [[float(x) for x in lines[i].strip().split("\t")] for i in range(S_START, S_START + cluster_count)]

	# Split depot.
	a.append(a[0])
	b.append(b[0])
	for i in range(0, n):
		d[i].append(d[i][0])
		d[i][0] = 0.0
	d.append([0.0 for j in range(0, n+1)])
	n = n+1

	for j in range(0, len(jams)):
		print("\t" + jam_files[j])
		jam = jams[j]

		# Compute speeds with jams.
		v = [[s[i][k] * jam[i][k] for k in range(0, speed_zone_count)] for i in range(0, cluster_count)]

		# Create instance.
		I = {
			"digraph": {
				"vertex_count": n,
				"arc_count": (n-1)*(n-2),
				"arcs": [[1 if (i != j and j != 0 and i != n-1 and not(i == 0 and j == n)) else 0 for j in range(0, n)] for i in range(0, n)]
			},
			"cluster_count": cluster_count,
			"cluster_speeds": v,
			"clusters": C,
			"speed_zone_count": speed_zone_count,
			"speed_zones": T,
			"distances": d,
			"time_windows": [[a[i], b[i]] for i in range(0, n)],
			"start_depot": 0,
			"end_depot": n-1,
			"horizon": [0.0, b[0]]
		}
		instance_number = file_name.replace("_new.txt", "").replace("15", "").replace("20", "").replace("25", "").replace("30", "").replace("35", "").replace("40", "").replace("45", "").replace("50", "").replace("55", "").replace("60", "").replace("_", "")
		instance_name = F"{n-2}_{jam_numbers[j]}_{jam_types[j]}_{instance_number}"
		with open(F"output/{instance_name}.json", "w") as ff:
			ff.write(json.dumps(I))

		index.append({"instance_name": instance_name, "file_name": F"{instance_name}.json", "tags": [instance_name, F"N{n-2}", F"J{jam_numbers[j]}", F"P{jam_types[j]}"]})

def sort_key(item):
	return item["instance_name"]

index.sort(key = sort_key)

# Create index.json
with open(F"output/index.json", "w") as ff:
		ff.write(json.dumps(index))

with open(F"output/solutions.json", "w") as ff:
		ff.write(json.dumps([]))