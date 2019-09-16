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
	afg_instance_name = entry["instance_name"]
	afg_instance = {}
	with open(F"afg_1995/{entry['file_name']}", "r") as ff:
		afg_instance = json.loads(ff.read())
	print(afg_instance)
	exit(0)