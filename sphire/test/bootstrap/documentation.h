#ifndef MYMPIDOC_H
#define MYMPIDOC_H
/******************************************************************************
 *
 *	File:			documentation.h
 *
 *	Function:		Documentation for mympi Python MPI module
 *
 *	Author(s):		<put your name here>
 *
 *	Copyright:		Copyright (c) 2006 The Regents of the University of California
 *					All Rights Reserved.
 *
 *	Source:			Original.
 *
 *	Notes:			
 *
 *	Change History:
 *			2006_oct_9	Started source.
 *	
 ******************************************************************************/

/*
$Date$
$Revision$
$Source$
*/

char DATE_DOC[]="$Date$";

#pragma once


char COPYWRITE_STR__[]="This fuction returns the Copywrite for the module.\n\nprint mpi.copywrite()";

char mpi_alltoall__[] = " Python version of MPI_Alltoall\n" \
"	Sends data from all to all processes\n" \
"recvbuf=mpi_alltoall(sendbuf,sendcount,sendtype,recvcount,recvtype,comm)\n" \
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	recvcount\n" \
"		number of elements received from any process (integer)\n" \
"\n" \
"	recvtype\n" \
"		data type of receive array elements (handle)\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_alltoallv__[] = "Python version of MPI_Alltoallv\n" \
"	Sends different amounts of data from all to all processes\n" \
"recvbuf=mpi_alltoall(sendbuf,sendcount,senddisp,sendtype,recvcount,recvdisp,recvtype,comm)\n" \
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	senddisp\n" \
"		integer array (of length group size). Entry j specifies the displacement\n" \
"		, relative to beginning of sendbuf, from which to take the outgoing data\n" \
"		destined for process j \n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	recvcount\n" \
"		number of elements received from any process (integer)\n" \
"\n" \
"	recvdisp\n" \
"		integer array (of length group size). Entry i specifies the displacement,\n" \
"       relative beginning of recvbuf, at which to place the incoming data from process i\n" \
"\n" \
"	recvtype\n" \
"		data type of receive array elements (handle)\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";



char mpi_barrier__[]= " Python version of MPI_Barrier\n" \
"	Blocks until all process in the communicator have reached this routine.\n"\
"mpi_barrier(comm)\n"\
"	comm\n" \
"		communicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_bcast__[]=	" Python version of MPI_Bcast\n" \
"	Broadcasts a message from the process with rank \"root\" to all other processes of the group.\n"\
"recvbuf=mpi_bast(sendbuf,sendcount,sendtype, root,comm )\n"\
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	root\n" \
"		process holding the data to be broadcast\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_send__[]=" Python version of MPI_Send\n" \
"	Send a message from a process to another.\n"\
"mpi_send(sendbuf,sendcount,sendtype, destination,tag,comm)\n"\
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	destination\n" \
"		process to which the data is being sent\n" \
"\n" \
"	tag\n" \
"		message tag\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"Sets error code read by mpi_error()";
char mpi_recv__[]=" Python version of MPI_Recv\n" \
"	Receive a message from nother process.\n"\
"recvbuf=mpi_send(recvcount,recvtype, source,tag,comm)\n"\
"	recvcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	recvtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	source\n" \
"		process from which the data is being received\n" \
"\n" \
"	tag\n" \
"		message tag\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets status as read by mpi_status()\n" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_start__[]=		    		"mpi_start - deprecated, use mpi_init";

char mpi_init__[]=" Python version of MPI_Init\n" \
"	Initialize the MPI library.\n"  \
"Takes as input the number of command line arguments and a list of\n"  \
"command line arguments.\n"  \
"\n"  \
"Returns a potentially modified array of command line arguments\n"  \
"\n"  \
"USAGE:\n"  \
"import mpi\n"  \
"import sys\n"  \
"sys.argv =  mpi.mpi_init(len(sys.argv),sys.argv)\n" \
"\n" \
"Sets error code read by mpi_error()";



char mpi_finalize__[]=" Python version of MPI_Finalize\n" \
"	Terminates MPI execution environment.  All processes must call this routine before exiting.\n"  \
"mpi_finalize()\n"\
"\n" \
"Sets error code read by mpi_error()";



char mpi_reduce__[]=" Python version of MPI_Reduce\n" \
"	Reduces values on all processes to a single valuer.\n"\
"recvbuf=mpi_send(sendbuf,sendcount,sendtype, op, root,comm)\n"\
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	op\n" \
"		reduce operation, MPI_LOR,MPI_LXOR,MPI_SUM,MPI_PROD,MPI_MIN,MPI_MAX\n" \
"\n" \
"	root\n" \
"		process to which the data is being reduced\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_scatter__[]=" Python version of MPI_Scatter\n" \
"	Sends different data from root to all processes\n" \
"recvbuf=mpi_scatter(sendbuf,sendcount,sendtype,recvcount,recvtype,root,comm)\n" \
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	recvcount\n" \
"		number of elements received from any process (integer)\n" \
"\n" \
"	recvtype\n" \
"		data type of receive array elements (handle)\n" \
"\n" \
"	root\n" \
"		process holding the data to be scattered\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_gather__[]=" Python version of MPI_Gather\n" \
"	Sends different data from all processes to the root process\n" \
"recvbuf=mpi_gather(sendbuf,sendcount,sendtype,recvcount,recvtype,root,comm)\n" \
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send from each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	recvcount\n" \
"		number of elements received from any process (integer)\n" \
"\n" \
"	recvtype\n" \
"		data type of receive array elements (handle)\n" \
"\n" \
"	root\n" \
"		process to which the data is to be gathered\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_scatterv__[]=" Python version of MPI_Scatterv\n" \
"	Sends different amounts of data from root to all processes\n" \
"recvbuf=mpi_scatterv(sendbuf,sendcount,senddisp,sendtype,recvcount,recvtype,root,comm)\n" \
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	senddisp\n" \
"		integer array (of length group size). Entry j specifies the displacement\n" \
"		, relative to beginning of sendbuf, from which to take the outgoing data\n" \
"		destined for process j \n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	recvcount\n" \
"		number of elements received from any process (integer)\n" \
"\n" \
"	recvtype\n" \
"		data type of receive array elements (handle)\n" \
"\n" \
"	root\n" \
"		process holding the data to be scattered\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_gatherv__[]=" Python version of MPI_Gatterv\n" \
"	Sends differing amounts of data from all processes to the root\n" \
"recvbuf=mpi_gatterv(sendbuf,sendcount,sendtype,recvcount,recvdisp,recvtype,root,comm)\n" \
"	sendbuf\n" \
"		send array (choice)\n" \
"\n" \
"	sendcount\n" \
"		number of elements to send to each process (integer)\n" \
"\n" \
"	sendtype\n" \
"		data type of send array elements (handle)\n" \
"\n" \
"	recvcount\n" \
"		number of elements received from any process (integer)\n" \
"\n" \
"	recvdisp\n" \
"		integer array (of length group size). Entry i specifies the displacement,\n" \
"       relative beginning of recvbuf, at which to place the incoming data from process i\n" \
"\n" \
"	recvtype\n" \
"		data type of receive array elements (handle)\n" \
"\n" \
"	root\n" \
"		process holding the data to be scattered\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	recvbuf\n" \
"		receive array\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_wtick__[]=" Python version of MPI_Wtick\n" \
"	Returns the resolution of mpi_wtime\n" \
"tres=mpi_wtick()\n" \
"\n" \
"	tres\n" \
"		Time in seconds of resolution of MPI_Wtime\n" \
"\n";


char mpi_wtime__[]=" Python version of MPI_Wtime\n" \
"	Returns an elapsed time on the calling processor\n" \
"thetime=mpi_wtime()\n" \
"\n" \
"	thetime\n" \
"		Time in seconds since an arbitrary time in the past.\n" \
"\n" \
"	This is intended to be a high-resolution, elapsed (or wall) clock.\n" \
"	See MPI_WTICK to determine the resolution of MPI_WTIME. If the attribute.\n" \
"	MPI_WTIME_IS_GLOBAL is defined and true, then the value is synchronized\n" \
"	across all processes in MPI_COMM_WORLD.\n" \
"\n";


char mpi_error__[]=" This routine is not part of the MPI standard\n" \
"	Returns an error code from the previous MPI call\n" \
"errcode=mpi_error()\n" \
"	errcode\n" \
"		The error code from the previous MPI call.\n";


char mpi_comm_rank__[]=" Python version of MPI_Comm_rank\n" \
"	Determines the rank of the calling process in the communicator\n" \
"myid=mpi_comm_rank(comm)\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	myid\n" \
"		rank of the calling process in group of comm\n";
char mpi_comm_size__[]=" Python version of MPI_Comm_size\n" \
"	Determines the size of the group associated with a communictor\n" \
"numprocs=mpi_comm_size(comm)\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	myid\n" \
"		number of processes in the group of comm\n" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_get_processor_name__[]=" Python version of MPI_Get_processor_name\n" \
"	Gets the name of the processor\n" \
"name=mpi_get_processor_name()\n" \
"\n" \
"	name\n" \
"		A unique specifier for the actual (as opposed to virtual) node.\n" \
"\n" \
"Sets error code read by mpi_error()";



char mpi_comm_create__[]=" Python version of MPI_Comm_create\n" \
"	Creates a new communicator\n" \
"comm_out=mpi_comm_create(comm,group)\n" \
"	comm\n" \
"		the old communicator\n" \
"\n" \
"	group\n" \
"		group, which is a subset of the group of comm\n" \
"\n" \
"	comm_out\n" \
"		the new communicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_comm_dup__[]=" Python version of MPI_Comm_dup\n" \
"	Duplicates an existing communicator with all its cached information\n" \
"comm_out=mpi_comm_dup(comm)\n" \
"	comm\n" \
"		the old communicator\n" \
"\n" \
"	comm_out\n" \
"		the new communicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_comm_group__[]=" Python version of MPI_Comm_group\n" \
"	Accesses the group associated with given communicator\n" \
"group_out=mpi_comm_group(comm)\n" \
"	comm\n" \
"		the  communicator\n" \
"\n" \
"	group_out\n" \
"		the group in the communicator\n" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_comm_split__[]=" Python version of MPI_Comm_split\n" \
"	Creates new communicators based on colors and keys\n" \
"comm_out=mpi_comm_split(comm,color,key)\n" \
"	comm\n" \
"		the  communicator\n" \
"\n" \
"	color\n" \
"		control of subset assignment (nonnegative integer). Processes with the same color are in the same new communicator\n" \
"\n" \
"	key\n" \
"		control of rank assigment\n" \
"\n" \
"	comm_out\n" \
"		the new communicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_status__[]=" This routine is not part of the MPI standard\n" \
"	Returns the status associated with the previous mpi_send or mpi_[i]probe call\n" \
"statray=mpi_status()\n" \
"	statray\n" \
"		the  status array\n" \
"		  statray[0]=MPI_SOURCE\n" \
"		  statray[1]=MPI_TAG\n" \
"		  statray[2]=MPI_ERROR";


char mpi_group_incl__[]=" Python version of MPI_Group_incl\n" \
"	Produces a group by reordering an existing group and taking only listed members\n" \
"group_out=mpi_group_incl(group,n,ranks)\n" \
"	group\n" \
"		the old group\n" \
"\n" \
"	n\n" \
"		number of elements in array ranks (and size of newgroup ) \n" \
"\n" \
"	ranks\n" \
"		array of ranks of processes in group to appear in newgroup\n" \
"\n" \
"	group_out\n" \
"		the group in the communicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_group_rank__[]=" Python version of MPI_Group_rank\n" \
"	Determines the rank of the calling process in the group\n" \
"myid=mpi_group_rank(group)\n" \
"	group\n" \
"		group\n" \
"\n" \
"	myid\n" \
"		rank of the calling process in group\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_get_count__[]=" Python version of MPI_Get_count\n" \
"	Gets the number of \"top level\" elements\n" \
"count=mpi_get_count(status,datatype)\n" \
"	status\n" \
"		status associated with the previous mpi_send or mpi_[i]probe call\n" \
"\n" \
"	datatype\n" \
"		datatype of each receive buffer element\n" \
"\n" \
"	count\n" \
"		number of received elements\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_probe__[]=" Python version of MPI_Probe\n" \
"	Blocking test for a message\n" \
"mpi_probe(source,tag,comm)\n"\
"	source\n" \
"		process from which the data is being received\n" \
"\n" \
"	tag\n" \
"		message tag\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"Sets status as read by mpi_status()";

char mpi_iprobe__[]=" Python version of MPI_Iprobe\n" \
"	Nonblocking test for a message\n" \
"flag=mpi_iprobe(source,tag,comm)\n"\
"	source\n" \
"		process from which the data is being received\n" \
"\n" \
"	tag\n" \
"		message tag\n" \
"\n" \
"	comm\n" \
"		communicator\n" \
"\n" \
"	flag\n" \
"		output flag indicates if a message has arrived\n" \
"\n" \
"Sets status as read by mpi_status()";


char mpi_attr_get__[]=" Python version of MPI_Attr_get\n" \
"	Retrieves attribute value by key\n" \
"attr_value=mpi_iprobe(comm,keyvalue)\n"\
"	comm\n" \
"		communicator\n" \
"\n" \
"	keyvalue\n" \
"		key value\n" \
"\n" \
"	attr_value\n" \
"		attribute value if it available, else this routine will cause an exceptionn" \
"\n" \
"Sets status as read by mpi_status()";



char mpi_comm_get_parent__[]=" Python version of MPI_Comm_get_parent\n" \
"	Return the parent communicator for this process\n" \
"parent=mpi.mpi_comm_get_parent()\n"\
"	parent\n" \
"		the parent communicator\n" \
"\n" \
"Remarks:\n" \
"	If a process was started with MPI_Comm_spawn or MPI_Comm_spawn_multiple,\n" \
"	MPI_Comm_get_parent returns the parent intercommunicator of the current\n" \
"	process. This parent intercommunicator is created implicitly inside of MPI_Init\n" \
"	and is the same intercommunicator returned by MPI_Comm_spawn in the parents.\n" \
"\n" \
"	If the process was not spawned, MPI_Comm_get_parent returns MPI_COMM_NULL.\n" \
"	After the parent communicator is freed or disconnected, MPI_Comm_get_parent\n" \
"	returns MPI_COMM_NULL.\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_comm_spawn__[]=" Python version of MPI_Comm_spawn\n" \
"	Spawn up to maxprocs instances of a single MPI application\n" \
"intercomm=mpi.mpi_comm_spawn(command,argv,maxprocs,info,root,comm)\n"\
"	command\n" \
"		arguments  to  command (array of strings, significant only at root)\n" \
"		As an extension, the Python version will replicate argv maxprocs times" \
"		if is a single string\n" \
"\n" \
"	argv\n" \
"		name of program to be spawned  (string,  significant  only  at root)\n" \
"\n" \
"	maxprocs\n" \
"		maximum number of processes to start\n" \
"\n" \
"	info\n" \
"		a set of key-value pairs telling the runtime system where and\n"
"		how to start the processes. normally MPI_INFO_NULL\n" \
"\n" \
"	root\n" \
"		rank of process  in  which  previous  arguments  are  examined\n" \
"\n" \
"	comm\n" \
"		the old communicator\n" \
"\n" \
"	intercomm\n" \
"		intercommunicator between original group and the newly spawned group\n" \
"\n" \
"Sets an array of error codes read by mpi_array_of_errcodes()" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_array_of_errcodes__[]=	" This routine is not part of the MPI standard\n" \
"	Returns an array of error codes from the previous call of mpi_comm_spawn\n" \
"errcodes=mpi_array_of_errcodes()\n" \
"	errcodes\n" \
"		The array of error codes from the previous call of mpi_comm_spawn.\n";

char mpi_comm_free__[]=" Python version of MPI_Comm_free\n" \
"	Marks the communicator object for deallocation\n" \
"mpi.mpi_comm_free(comm)\n"\

"	comm\n" \
"		the old communicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_open_port__[]=" Python version of MPI_Open_port\n" \
"	Establish  an  address that can be used to establish\n" \
"	connections between groups of MPI processes\n" \
"port_name=mpi.mpi_open_port(info)\n"\

"	info\n" \
"		implementation-specific information on how to establish a port\n" \
"		normally MPI_INFO_NULL\n" \
"\n" \
"	port_name\n" \
"		newly established port (string)\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_intercomm_merge__[]=" Python version of MPI_Intercomm_merge\n" \
"	 Creates an intracommuncator from an intercommunicator\n" \
"	 connections between groups of MPI processes\n" \
"comm_out=mpi.mpi_intercomm_merge(intercomm,high)\n"\
"	intercomm\n" \
"		Intercommunicator\n" \
"\n" \
"	high\n" \
"		Used to order the groups within comm (logical)  when  creating\n" \
"		the  new  communicator.  This is a boolean value; the group that\n" \
"		sets high true has its processes ordered after  the  group  that\n" \
"		sets this value to false.  If all processes in the intercommuni-\n" \
"		cator provide the same value,  the  choice  of  which  group  is\n" \
"		ordered first is arbitrary.\n" \
"\n" \
"	comm_out\n" \
"		Created intracommunicator\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_close_port__[]=" Python version of MPI_Close_port\n" \
"	 close port\n" \
"mpi.mpi_close_port(port_name)\n"\
"	port_name\n" \
"		a port name (string)\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_comm_disconnect__[]=" Python version of MPI_Comm_disconnect\n" \
"	 Disconnect from a communicator\n" \
"mpi.mpi_comm_disconnect(comm)\n"\
"	comm\n" \
"		communicator\n" \
"\n" \
"Sets error code read by mpi_error()";



char mpi_comm_accept__[]=" Python version of MPI_Comm_accept\n" \
"	 Accept a request to form a new intercommunicator\n" \
"newcomm=mpi.mpi_comm_accept(port_name,info,root,comm)\n"\
"	port_name\n" \
"		a port name (string used only on root)\n" \
"\n" \
"	info\n" \
"		implementation-dependent  information  (handle,  used only on root)\n" \
"		normally MPI_INFO_NULL\n" \
"\n" \
"	root\n" \
"		rank in comm of root node\n" \
"\n" \
"	comm\n" \
"		intracommunicator over which call is collective\n" \
"\n" \
"	newcomm\n" \
"		intercommunicator with client as remote group\n" \
"\n" \
"Sets error code read by mpi_error()";


char mpi_comm_connect__[]=" Python version of MPI_Comm_connect\n" \
"	 Make a request to form a new intercommunicator\n" \
"mpi.mpi_comm_connect(comm,flag)\n"\
"	port_name\n" \
"		network address (string, used only on root)\n" \
"\n" \
"	info\n" \
"		implementation-dependent  information  (handle,  used only on root)\n" \
"		normally MPI_INFO_NULL\n" \
"\n" \
"	root\n" \
"		rank in comm of root node\n" \
"\n" \
"	comm\n" \
"		intracommunicator over which call is collective\n" \
"\n" \
"	newcomm\n" \
"		intercommunicator with client as remote group\n" \
"\n" \
"Sets error code read by mpi_error()";

char mpi_comm_set_errhandler__[]=" Python version of MPI_Comm_set_errhandler\n" \
"	 Sets the behavior if an MPI error occurs using the given communicator\n" \
"newcomm=mpi.mpi_comm_set_errhandler(comm,flag)\n"\
"	comm\n" \
"		intracommunicator over which call is collective\n" \
"\n" \
"	flag(0-2)\n" \
"		set the type of error handler\n" \
"			0 = MPI_ERRORS_ARE_FATAL, an MPI error will cause the program\n" \
"				to exit calling MPI_Abort\n" \
"			1 = MPI_ERRORS_RETURN, MPI call returns setting the error code\n" \
"			2 = The error code is printed to stdout.MPI call returns\n" \
"				setting the error code.\n" \
"Note:\n" \
"	This is not the way the C and Fortran MPI_Comm_set_errhandler call\n" \
"   is made.  For the original version the second argument is a address\n" \
"	of a routine to be called on error.  This might be implemented in the future.\n"
"\n" \
"Sets error code read by mpi_error()";





char mpi_irecv__[]=		    		"mpi_irecv  not yet implemented";
char mpi_isend__[]=		    		"mpi_isend  not yet implemented";
char mpi_test__[]=		     		"mpi_test  not yet implemented";
char mpi_wait__[]=		     		"mpi_wait  not yet implemented";


#endif /* MYMPIDOC_H */
