
def parse_filter_params(filterparams):
    params = filterparams.split(":")
    filtername = params[0]

    if len(params) == 1:
        return (filtername, None)
    else:
        d = Dict()
        for param in params[1:]:
            key_values = param.split("=")
            d[key_values[0]] = EMObject(key_values[1])
        return (filtername, d)


def get_optionlist(argv):
    optionlist = []
    for arg1 in argv:
        if arg1[0] == "-":
            argname = arg1.split("=")
            optionlist.append(argname[0].lstrip("-"))
    return optionlist
