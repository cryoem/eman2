// mpi_eman.c - Steve Ludtke 11/23/2010

// Simple MPI bindings for Python for use with EMAN2

#include <python2.6/Python.h>
#include <math.h>
#include <mpi.h>

/*
Calls MPI_Init and returns the current task's rank
*/
static PyObject *mpi_init(PyObject *self,PyObject *args) {
MPI_Init(NULL,NULL);

int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

return Py_BuildValue("i",myrank);
}

static PyObject *mpi_comm_rank(PyObject *self,PyObject *args) {
int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

return Py_BuildValue("i",myrank);
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

return NULL;
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


static PyObject *mpi_bcast(PyObject *self,PyObject *args) {

	return NULL;
}

static PyObject * mpi_barrier(PyObject *self,PyObject *args) {

	return NULL;
}

static PyObject * mpi_finalize(PyObject *self,PyObject *args) {
	MPI_Finalize();
	return NULL;
}

static PyMethodDef EmanMpiMethods[] = {
	{"mpi_init",mpi_init,METH_VARARGS,"MPI_Init command. No arguments. Returns the rank id to each node."},
	{"mpi_comm_rank",mpi_comm_rank,METH_VARARGS,"This will return the rank id, same as returned by mpi_init."},
	{"mpi_send",mpi_send,METH_VARARGS,"MPI_Send(string,destination rank,tag)"},
	{"mpi_recv",mpi_recv,METH_VARARGS,"MPI_Recv(source rank,tag). If either is negative, arbitrary values accepted. Returns the received string."},
	{"mpi_bcast",mpi_bcast,METH_VARARGS,"MPI_Bcast(string). Pass data on source node. Pass None on others & return data."},
	{"mpi_barrier",mpi_barrier,METH_VARARGS,"MPI_Barrier(). No arguments or return. Blocks until all nodes call it."},
	{"mpi_finalize",mpi_finalize,METH_VARARGS,"MPI_Finalize(). No arguments or return. Call before exiting an MPI Python program exactly once."},
	{NULL,NULL,0,NULL}
};

PyMODINIT_FUNC initmpi_eman(void) {
	(void) Py_InitModule("mpi_eman", EmanMpiMethods);
}
	