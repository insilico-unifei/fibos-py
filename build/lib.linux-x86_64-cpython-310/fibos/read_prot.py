def read_OS(file):
	filter = []
	file_read = open("prot.srf","r")
	for line in file_read:
		if "INF" in line:
			filter.append(line)
	return(filter)

if __name__ == '__main__':
	file = read_OS("prot.srf")
	for i in range(0,len(file)):
		print(file[i])
