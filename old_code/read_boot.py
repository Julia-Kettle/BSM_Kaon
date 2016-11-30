import readfiles as rf


def boot_reader(filename):
	n ,m, data = rf.read_results_Jamie(filename)
	print data
	print data

def logfileread(filename,param):
	fi = open(filename,'r')
	lines = fi.readlines()
	data = []
	error = 0
	save=0
	count=0
	trig2 = "params[" + str(param) + "]"
	trig1 = trig2 + ":"
	print trig1, trig2
	print trig1, trig2
	for line in lines:
		if count == 500:
			print "Hit last bootstrap"
			print "Hit last bootstrap"
			save=0 #have reached last bootstrap - stop saving
		if save == 1:
			print "saving a bootstrap"
			print "saving a bootstrap"
			data.append(((line.split(' '))[1])[:-2])
			count = count + 1 
		if line.startswith(trig1):
			print "Hit first bootstrap"
			print "Hit first bootstrap"
			save=1
		if count >= 500 and line.startswith(trig2):
			print "Hit central Value"
			print "Hit central Value"
			data.append((line.split(' '))[2]) #central value
			error = ((line.split(' '))[3][:-2])

	
	print data
	print data
	print error
	print error
		

def main():
	#boot_reader("../Fits/48cubed/mass/boots/mass-s0.0362-l0.00078.boot")
	logfileread("../Fits/48cubed/mass/logs/mass-s0.0362-l0.00078.log",1)
main()
