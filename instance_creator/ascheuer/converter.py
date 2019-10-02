import os, json, random

speed_profiles = [
	[1.5, 1.0, 1.67, 1.17, 1.33],  # fast
	[1.17, 0.67, 1.33, 0.83, 1.0], # normal
	[1.0, 0.33, 0.67, 0.5, 0.83]   # slow
]
zone_count = 5
cluster_count = 3
# cluster_matrix = [[random.randrange(cluster_count) for j in range(0, max_n)] for i in range(0, max_n)]

# Read AFG instance index.
index_afg = []
with open("afg_1995/index.json", "r") as ff:
	index_afg = json.loads(ff.read())

for entry in index_afg:
	instance_name = entry["instance_name"]
	afg_instance = {}
	with open(F"afg_1995/{entry['file_name']}", "r") as ff:
		afg_instance = json.loads(ff.read())

	n = afg_instance["digraph"]["vertex_count"]
	T = afg_instance["time_windows"][-1][1]
	speed_zones = [[0.0, 0.2*T], [0.2*T,0.3*T], [0.3*T, 0.7*T], [0.7*T,0.8*T], [0.8*T, T]]

	afg_instance["cluster_count"] = cluster_count
	afg_instance["cluster_speeds"] = speed_profiles
	afg_instance["clusters"] = [[random.randrange(cluster_count) for j in range(0, n)] for i in range(0, n)]
	afg_instance["speed_zone_count"] = zone_count
	afg_instance["speed_zones"] = speed_zones
	for i in range(0, n):
		for j in range(0, n):
			if afg_instance["distances"][i][j] == 0.0:
				afg_instance["distances"][i][j] = 1.0

	with open(F"lms_2019/{instance_name}.json", "w") as ff:
		ff.write(json.dumps(afg_instance))
