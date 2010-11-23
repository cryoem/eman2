// mpi_eman.c - Steve Ludtke 11/23/2010

// Simple MPI bindings for Python for use with EMAN2

#include <python2.6/Python.h>
#include <math.h>
#include <mpi.h>

/*
Calls MPI_Init and returns the current task's rank
*/
static PyObject *MPI_Init(PyObject *self,PyObject *args) {
MPI_Init(NULL,NULL);

int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

return Py_BuildValue("i",myrank);
}

static PyObject *MPI_Comm_rank(PyObject *self,PyObject *args) {
int myrank;
MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

return Py_BuildValue("i",myrank);
}


/*
Calls MPI_Send. Expects (string,dest rank(int),tag (int))
Note that if tag is negative, receiver will not be constrained to the specific tag.
*/
static PyObject *MPI_Send(PyObject *self,PyObject *args) {
const char * data;
int datalen,dest,tag;

// Get a string argument
if (!PyArg_ParseTuple(args,"s#ii",&data,&datalen,&dest,&tag)) return NULL;

MPI_Send(data,datalen,MPI_CHAR,dest,tag,MPI_COMM_WORLD);

return NULL;
}

/*
Calls MPI_Recv. Expects (source rank(int),tag). 
If negative, ANY source/tag will be acceptable
Returns (data,source rank,tag)
*/
static PyObject *MPI_Recv(PyObject *self,PyObject *args) {
MPI_Status status;
const char * data;
int datalen,src,tag;

// Get a string argument
if (!PyArg_ParseTuple(args,"ii",&src,&tag)) return NULL;

if (src<0) src=MPI_ANY_SOURCE;
if (tag<0) tag=MPI_ANY_TAG;

MPI_PROBE(src,tag,MPI_COMM_WORLD,&status);

MPI_GET_COUNT(&status,MPI_CHAR,&datalen);
data=(const char *)malloc(datalen);
MPI_Recv(data,datalen,MPI_CHAR,src,tag,MPI_COMM_WORLD,&status);

PyObject *ret = Py_BuildValue("s#ii",data,datalen,status.MPI_SOURCE,status.MPI_TAG);
free(data);

return ret;
}


static PyObject *MPI_Broadcast(PyObject *self,PyObject *args) {

}
