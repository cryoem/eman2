


































































































































































































































































































































































































































































































































































































































































































































































































































"""0
	# For DEBUG
	if sxcmd.name == "sxwindow":
		print "><><>< DEBUG OUTPUT ><><><"
		print ""
		print "------"
		print "GLOBAL"
		print "------"
		print "name            : %s" % sxcmd.name
		print "subname         : %s" % sxcmd.subname
		print "label           : %s" % sxcmd.label
		print "short_info      : %s" % sxcmd.short_info
		print "mpi_support     : %s" % sxcmd.mpi_support
		print "mpi_add_flag    : %s" % sxcmd.mpi_add_flag
		print "type            : %s" % sxcmd.type
		print "len(token_list) : %d" % len(sxcmd.token_list)
		print "len(token_dict) : %d" % len(sxcmd.token_dict)
		print ""
		print "--------------"
		print "cmd_token_list"
		print "--------------"
		for token in sxcmd.token_list:
			print "%s%s (group=%s, required=%s, default=%s, type=%s, restore=%s) <%s>" % (token.key_prefix, token.key_base, token.group, token.is_required, token.default, token.type, token.restore, token.label, token.help)
		print ""
	"""















































































































































































































































































































































































































































"""1
	# For DEBUG
	if sxcmd.name in ["sxrelion2sphire", "sxpipe"]:
		print("><><>< DEBUG OUTPUT ><><><")
		print("")
		print("------")
		print("GLOBAL")
		print("------")
		print("name            : \'{}\'".format(sxcmd.name))
		print("subname         : \'{}\'".format(sxcmd.subname))
		print("mode            : \'{}\'".format(sxcmd.mode))
		print("subset_config   : \'{}\'".format(sxcmd.subset_config))
		print("label           : \'{}\'".format(sxcmd.label))
		print("short_info      : \'{}\'".format(sxcmd.short_info))
		print("mpi_support     : \'{}\'".format(sxcmd.mpi_support))
		print("mpi_add_flag    : \'{}\'".format(sxcmd.mpi_add_flag))
		print("category        : \'{}\'".format(sxcmd.category))
		print("role            : \'{}\'".format(sxcmd.role))
		print("is_submittable  : \'{}\'".format(sxcmd.is_submittable))
		print("len(token_list) : {}".format(len(sxcmd.token_list)))
		print("len(token_dict) : {}".format(len(sxcmd.token_dict)))
		print("")
		print("--------------")
		print("cmd_token_list")
		print("--------------")
		for token in sxcmd.token_list:
			print("\'{}{}\' (group=\'{}\', is_required=\'{}\', is_locked=\'{}\', is_reversed=\'{}\', default=\'{}\', restore=\'{}\', type=\'{}\', label=\'{}\', help=\'{}\'".format(token.key_prefix, token.key_base, token.group, token.is_required, token.is_locked, token.is_reversed, token.default, token.restore, token.type, token.label, token.help))
		print("")
	"""


















































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































