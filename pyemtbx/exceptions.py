def exception_type(runtime_err):
	fullerr = runtime_err.args[0]
	first_space = fullerr.find(" ")
	firstword = fullerr[0:first_space]
	return firstword

