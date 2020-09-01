import os, json, random, math

speed_zone_count = 73
cluster_count = 3

# Load Jams.
sizes = [60, 80, 100]
jam_files = ["new_Jam_70_C1", "new_Jam_70_C3", "new_Jam_80_C1", "new_Jam_80_C3", "new_Jam_90_C1", "new_Jam_90_C3", "new_Jam_98_C1", "new_Jam_98_C3"]
jam_names = ["j70c1", "j70c3", "j80c1", "j80c3", "j90c1", "j90c3", "j98c1", "j98c3"]
jam_numbers = ["70_A", "70_B", "80_A", "80_B", "90_A", "90_B", "98_A", "98_B"]
widths = [40, 60, 80, 100, 120, 150, 180]
numbers = [0,1,2,3,4]

jams = []
for jam_file in jam_files:
	with open(F"instances/{jam_file}", "r") as ff:
		print(jam_file)
		jams.append([[float(x) for x in line.strip().split("\t")] for line in ff.read().split("\n") if line != ""])

# Read instances.
index = []
for S in sizes:
	for J in range(len(jam_names)):
		for W in widths:
			for I in numbers:
				instance = {}
				with open(F"instances/n{S}{jam_names[J]}_w{W}_{I}.txt", "r") as ff:
					instance = ff.read()
				print(F"instances/n{S}{jam_names[J]}_w{W}_{I}.txt")

				lines = instance.split("\n")
					
				# Parse instance file.
				n = int(lines[0])
				d = [[float(x) for x in lines[i].strip().split(" ")] for i in range(1, n+1)]
				a = [float(lines[i].strip().replace(" ", "\t").split("\t")[0]) for i in range(n+1, n+1+n)]
				b = [float(lines[i].strip().replace(" ", "\t").split("\t")[1]) for i in range(n+1, n+1+n)]
				C_START = n+1+n+1
				C = [[int(x)-1 for x in lines[i].strip().split(" ")] for i in range(C_START, C_START+n+1)]
				T_START = C_START+n+1+1
				T = [[float(x) for x in lines[i].strip().replace(" ", "\t").split("\t")] for i in range(T_START, T_START + speed_zone_count)]
				S_START = T_START + speed_zone_count + 1
				s = [[float(x) for x in lines[i].strip().replace(" ", "\t").split("\t")] for i in range(S_START, S_START + cluster_count)]

				# Split depot.
				a.append(a[0])
				b.append(b[0])
				for i in range(0, n):
					d[i].append(d[i][0])
					d[i][0] = 0.0
				d.append([0.0 for j in range(0, n+1)])
				n = n+1

				jam = jams[J]

				# Compute speeds with jams.
				v = [[s[i][k] * jam[i][k] for k in range(0, speed_zone_count)] for i in range(0, cluster_count)]

				# Create instance.
				instance = {
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
				instance_name = F"{S}_{jam_numbers[J]}_{W}_{I}"
				with open(F"output/{instance_name}.json", "w") as ff:
					ff.write(json.dumps(instance))

				index.append({"instance_name": instance_name, "file_name": F"{instance_name}.json", "order": (S, jam_numbers[J], W, I), "tags": [instance_name, F"N{S}", F"J{jam_numbers[J]}", F"W{W}"]})

def sort_key(item):
	return item["order"]

index.sort(key = sort_key)

for item in index:
	del item["order"]

# Create index.json
with open(F"output/index.json", "w") as ff:
		ff.write(json.dumps(index))

with open(F"output/solutions.json", "w") as ff:
		ff.write(json.dumps([]))