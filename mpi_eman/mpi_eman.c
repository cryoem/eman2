// mpi_eman.c - Steve Ludtke 11/23/2010

// Simple MPI bindings for Python for use with EMAN2

#include <python2.6/Python.h>
#include <math.h>
#include <mpi.h>

/*
Calls MPI_Init and returns the current task's rank, and the group size
*/
static PyObject *mpi_init(PyObject *self,PyObject *args) {
MPI_Init(NULL,NULL);

int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

int size;
MPI_Comm_size (MPI_COMM_WORLD, &size);

return Py_BuildValue("ii",myrank,size);
}

static PyObject *mpi_comm_rank(PyObject *self,PyObject *args) {
int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

return Py_BuildValue("i",myrank);
}

static PyObject *mpi_comm_size(PyObject *self,PyObject *args) {
int size;
MPI_Comm_size(MPI_COMM_WORLD, &size);

return Py_BuildValue("i",size);		
}

/*
Calls MPI_Send. Expects (string,dest rank(int),tag (int))
Note that if tag is negative, receiver will not be constrained to the specific tag.
*/
static PyObject *mpi_send(PyObject *self,PyObject *args) {
const char * data;
int datalen,dest,tag;

// Get a string argument
if (!PyArg_ParseTuple(args,"s#ii",&data,&datalen,&dest,&tag)) return NULL;

MPI_Send((void *)data,datalen,MPI_CHAR,dest,tag,MPI_COMM_WORLD);

	Py_RETURN_NONE;

}

/*
Calls MPI_Recv. Expects (source rank(int),tag). 
If negative, ANY source/tag will be acceptable
Returns (data,source rank,tag)
*/
static PyObject *mpi_recv(PyObject *self,PyObject *args) {
MPI_Status status;
const char * data;
int datalen,src,tag;

// Get a string argument
if (!PyArg_ParseTuple(args,"ii",&src,&tag)) return NULL;

if (src<0) src=MPI_ANY_SOURCE;
if (tag<0) tag=MPI_ANY_TAG;

MPI_Probe(src,tag,MPI_COMM_WORLD,&status);

MPI_Get_count(&status,MPI_CHAR,&datalen);
data=(const char *)malloc(datalen);
MPI_Recv((void *)data,datalen,MPI_CHAR,src,tag,MPI_COMM_WORLD,&status);

PyObject *ret = Py_BuildValue("s#ii",(void *)data,datalen,status.MPI_SOURCE,status.MPI_TAG);
free((void *)data);

return ret;
}

/*
Calls MPI_Probe. Expects (source rank(int),tag). 
If negative, ANY source/tag will be acceptable
Returns (length,source rank,tag)

Useful for waiting for input from an unknown source
*/
static PyObject *mpi_probe(PyObject *self,PyObject *args) {
MPI_Status status;
int datalen,src,tag;

// Get a string argument
if (!PyArg_ParseTuple(args,"ii",&src,&tag)) return NULL;

if (src<0) src=MPI_ANY_SOURCE;
if (tag<0) tag=MPI_ANY_TAG;

MPI_Probe(src,tag,MPI_COMM_WORLD,&status);
MPI_Get_count(&status,MPI_CHAR,&datalen);

PyObject *ret = Py_BuildValue("iii",datalen,status.MPI_SOURCE,status.MPI_TAG);

return ret;
}

/*
Calls MPI_Iprobe (nonblocking). Expects (source rank(int),tag). 
If negative, ANY source/tag will be acceptable
Returns (length,source rank,tag)

Returns None immediately if no suitable receives are pending
*/
static PyObject *mpi_iprobe(PyObject *self,PyObject *args) {
MPI_Status status;
int datalen,src,tag,flag;

// Get a string argument
if (!PyArg_ParseTuple(args,"ii",&src,&tag)) return NULL;

if (src<0) src=MPI_ANY_SOURCE;
if (tag<0) tag=MPI_ANY_TAG;

MPI_Iprobe(src,tag,MPI_COMM_WORLD,&flag,&status);
if (!flag) Py_RETURN_NONE;

MPI_Get_count(&status,MPI_CHAR,&datalen);

PyObject *ret = Py_BuildValue("iii",datalen,status.MPI_SOURCE,status.MPI_TAG);

return ret;
}


/* Calls MPI_Bcast
 * Different args on source and destinations. 
 * On source, pass in a string to transmit.
 * On destinations, pass in the data length and source rank.
 * The destinations must know the length before transmission. ALL nodes in the group must call this function
 * before it will complete.
*/
static PyObject *mpi_bcast_send(PyObject *self,PyObject *args) {
	const char * data;
	int datalen;
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	PyArg_ParseTuple(args,"s#",&data,&datalen);
//			printf("Sent: %d\n",datalen);
	MPI_Bcast((void *)data,datalen,MPI_CHAR,myrank,MPI_COMM_WORLD);
	Py_RETURN_NONE;
}

/* Calls MPI_Bcast
 * Different args on source and destinations. 
 * On source, pass in a string to transmit.
 * On destinations, pass in the data length and source rank.
 * The destinations must know the length before transmission. ALL nodes in the group must call this function
 * before it will complete.
*/
static PyObject *mpi_bcast_recv(PyObject *self,PyObject *args) {
	const char * data;
	int datalen;
	int root;

	PyArg_ParseTuple(args,"ii",&datalen,&root);
	data=(const char *)malloc(datalen);
	MPI_Bcast((void *)data,datalen,MPI_CHAR,root,MPI_COMM_WORLD);

	PyObject *ret = Py_BuildValue("s#",(void *)data,datalen);
	free((void *)data);

	return ret;
}


static PyObject * mpi_barrier(PyObject *self,PyObject *args) {

	MPI_Barrier(MPI_COMM_WORLD);
	Py_RETURN_NONE;
}

static PyObject * mpi_finalize(PyObject *self,PyObject *args) {
	MPI_Finalize();
	Py_RETURN_NONE;
}

static PyMethodDef EmanMpiMethods[] = {
	{"mpi_init",mpi_init,METH_VARARGS,"MPI_Init command. No arguments. Returns the rank id and process count to each node."},
	{"mpi_comm_rank",mpi_comm_rank,METH_VARARGS,"This will return the rank id, same as returned by mpi_init."},
	{"mpi_comm_size",mpi_comm_size,METH_VARARGS,"This will return the number of processes, same as returned by mpi_init."},
	{"mpi_send",mpi_send,METH_VARARGS,"MPI_Send(string,destination rank,tag)"},
	{"mpi_recv",mpi_recv,METH_VARARGS,"MPI_Recv(source rank,tag). If either is negative, arbitrary values accepted. Returns (data,src,tag)."},
	{"mpi_probe",mpi_probe,METH_VARARGS,"MPI_Probe(source rank,tag). If either is negative, arbitrary values accepted. Returns (len,src,tag)."},
	{"mpi_iprobe",mpi_iprobe,METH_VARARGS,"MPI_Iprobe(source rank,tag). If either is negative, arbitrary values accepted. Returns (len,src,tag) or None if no recieves are pending."},
	{"mpi_bcast_send",mpi_bcast_send,METH_VARARGS,"MPI_Bcast(string). Provide a string to broadcast to other nodes"},
	{"mpi_bcast_recv",mpi_bcast_recv,METH_VARARGS,"MPI_Bcast(data_length,source). Receives a broadcast of length data_length from source and returns a string."},
	{"mpi_barrier",mpi_barrier,METH_VARARGS,"MPI_Barrier(). No arguments or return. Blocks until all nodes call it."},
	{"mpi_finalize",mpi_finalize,METH_VARARGS,"MPI_Finalize(). No arguments or return. Call before exiting an MPI Python program exactly once."},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initmpi_eman_c(void) {
	(void) Py_InitModule("mpi_eman_c", EmanMpiMethods);
}
	