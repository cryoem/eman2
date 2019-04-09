/*
Copyright (c) 2005 The Regents of the University of California
All Rights Reserved
Permission to use, copy, modify and distribute any part of this
	Continuity 6 for educational, research and non-profit purposes,
	without fee, and without a written agreement is hereby granted,
	provided that the above copyright notice, this paragraph and the
	following three paragraphs appear in all copies.
Those desiring to incorporate this software into commercial products or
	use for commercial purposes should contact the Technology
	Transfer & Intellectual Property Services, University of
	California, San Diego, 9500 Gilman Drive, Mail Code 0910, La
	Jolla, CA 92093-0910, Ph: (858) 534-5815, FAX: (858) 534-7345,
	E-MAIL:invent@ucsd.edu.
IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
	PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
	CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
	THE USE OF THIS SOFTWARE  EVEN IF THE UNIVERSITY OF
	CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
THE SOFTWARE PROVIDED HEREIN IS ON AN AS IS BASIS, AND THE
	UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE
	MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
	THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND
	EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS,
	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
	MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE
	USE OF THE Continuity WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
	OTHER RIGHTS.";
*/


/*
$Date$
$Revision$
$Source$
*/

/************** numpy or Numeric? **************/


/* if NUMPY is defined then we look for arrayobject.h
   in <arrayobject.h>.  If not then in <Numeric/arrayobject.h>

   The second case would be the norm if you are linking against
   Numeric instead of numpy.

   We default to Numeric but this might change in the future.

   You can see the which was used at run time by printing
   mpi.ARRAY_LIBRARY.

*/

/************** numpy or Numeric **************/
char DATE_SRC[]="$Date$";

#define VERSION "1.13.0"
#define COPYWRITE "Copyright (c) 2005 The Regents of the University of California All Rights Reserved.  print mpi.copywrite() for details."
#undef DO_UNSIGED
#define DATA_TYPE long
#define COM_TYPE long
#define ARG_ARRAY
/* #define ARG_STR */
/* #define SIZE_RANK */
#include <Python.h>
#include "documentation.h"

#ifdef NUMPY
#include <arrayobject.h>
#define LIBRARY "NUMPY"
#else
#include <Numeric/arrayobject.h>
#define LIBRARY "Numeric"
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

#include <string.h>
#include <mpi.h>

static int calls = 0;
static int errs = 0;
MPI_Errhandler newerr;
static MPI_Comm mycomm;



#ifdef MPI_VERSION
#if (MPI_VERSION == 2)
#define MPI2
#endif
#else
#define MPI_VERSION 1
#endif

#ifndef MPI_SUBVERSION
#define MPI_SUBVERSION 0
#endif

char cw[2160];

#define com_ray_size 20
MPI_Comm com_ray[com_ray_size];
char errstr[256];
char version[8];

#ifdef MPI2
int *array_of_errcodes=NULL;
int array_of_errcodes_size=0;
#endif

MPI_Status status;
int ierr;
void char_func (char *ret, int retl, char *str, int slen, char *str2, int slen2, int *offset);
void the_date(double *since,char* date_str);
void myerror(char *s);
int getptype(long mpitype);
int erroron;

static PyObject *mpiError;

void eh( MPI_Comm *comm, int *err, ... )
{
char string[256];
int len;

/*
 if (*err != MPI_ERR_OTHER) {
 errs++;
 printf( "Unexpected error code\n" );fflush(stdout);
 }
 if (*comm != mycomm) {
 errs++;
 printf( "Unexpected communicator\n" );fflush(stdout);
 }
 calls++;
 */
 ierr=*err;
 /*
 MPI_Error_string(ierr, string,  &len);
 printf( "mpi generated the error %d %s\n",*err,string );fflush(stdout);
 */
 printf( "mpi generated the error %d\n",*err );fflush(stdout);
 return;
}


static PyObject *mpi_get_processor_name(PyObject *self, PyObject *args)
{
/* int MPI_Get_processor_name( char *name, int *resultlen) */
char c_name[MPI_MAX_PROCESSOR_NAME];
int c_len;
        ierr=MPI_Get_processor_name((char *)c_name,&c_len);
      return PyString_FromStringAndSize(c_name,c_len);
}

int getptype(long mpitype) {
	if(mpitype == (long)MPI_INT)    return(PyArray_INT);
	if(mpitype == (long)MPI_FLOAT)  return(PyArray_FLOAT);
	if(mpitype == (long)MPI_DOUBLE) return(PyArray_DOUBLE);
	printf("could not find type input: %ld  available: MPI_FLOAT %ld MPI_INT %ld MPI_DOUBLE %ld\n",mpitype,(long)MPI_FLOAT,(long)MPI_INT,(long)MPI_DOUBLE);
	return(PyArray_INT);
}

static PyObject *mpi_test(PyObject *self, PyObject *args)
{
/* int MPI_Test (MPI_Request  *request,int *flag, MPI_Status *status) */
	printf("this routine does not work yet\n");
    return NULL;
}

static PyObject *mpi_wait(PyObject *self, PyObject *args)
{
/* int MPI_Wait (MPI_Request  *request, MPI_Status *status) */
	printf("this routine does not work yet\n");
    return NULL;
}
static PyObject *mpi_isend(PyObject *self, PyObject *args)
{
/* int MPI_Isend( void *buf, int count, MPI_Datatype datatype, int dest, int tag,MPI_Comm comm, MPI_Request *request ) */
	printf("this routine does not work yet\n");
    return NULL;
}
static PyObject *mpi_irecv(PyObject *self, PyObject *args)
{
/* int MPI_Irecv( void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request ) */
	printf("this routine does not work yet\n");
    return NULL;
}

static PyObject *mpi_group_rank(PyObject *self, PyObject *args)
{
/* int MPI_Group_rank ( MPI_Group group, int *rank ) */
long in_group;
int rank;
MPI_Group group;

	if (!PyArg_ParseTuple(args, "l", &in_group))
		return NULL;
	group=(MPI_Group)in_group;
	ierr=MPI_Group_rank (group, &rank );
	return PyInt_FromLong((long)rank);
}

static PyObject *mpi_group_incl(PyObject *self, PyObject *args)
{
/* int MPI_Group_incl ( MPI_Group group, int n, int *ranks, MPI_Group *group_out ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
long in_group;
int *ranks,n;
MPI_Group group,out_group;
PyObject *ranks_obj;
PyArrayObject *array;

	if (!PyArg_ParseTuple(args, "liO", &in_group,&n,&ranks_obj))
        return NULL;
    group=(MPI_Group)in_group;
	array = (PyArrayObject *) PyArray_ContiguousFromObject(ranks_obj, PyArray_INT, 1, 1);
	if (array == NULL)
		return NULL;
	if(array->dimensions[0] < n)
		return NULL;
	ranks=(int*)malloc((size_t) (n*sizeof(int)));
	memcpy((void*)ranks,(void *)(array->data),  (size_t) (n*sizeof(int)));
	ierr=MPI_Group_incl ( (MPI_Group )group, n, ranks, &out_group);
	free(ranks);
	Py_DECREF(array);
	return VERT_FUNC((CAST)out_group);
}

static PyObject *mpi_comm_group(PyObject *self, PyObject *args)
{
/* int MPI_Comm_group ( MPI_Comm comm, MPI_Group *group ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
MPI_Group group;
COM_TYPE comm;
int ierr;

	if (!PyArg_ParseTuple(args, "l", &comm))
        return NULL;
    if((sizeof(MPI_Group) != sizeof(long))  &&  (sizeof(MPI_Group) != sizeof(int)))
    	printf("can not return MPI_Group as long or int sizes %ld %ld %ld",
    	                    sizeof(MPI_Group),sizeof(long),sizeof(int));
    ierr=MPI_Comm_group ( (MPI_Comm) comm, &group );
	return VERT_FUNC((CAST)group);
}

static PyObject *mpi_comm_dup(PyObject *self, PyObject *args)
{
/* int MPI_Comm_dup ( MPI_Comm comm, MPI_Comm *comm_out ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
long tmpcomm;
COM_TYPE  incomm,outcomm;

	if (!PyArg_ParseTuple(args, "l", &tmpcomm))
        return NULL;
    incomm=(COM_TYPE)tmpcomm;
	ierr=MPI_Comm_dup ((MPI_Comm)incomm,(MPI_Comm*)&outcomm );
	return VERT_FUNC((CAST)outcomm);
}


static PyObject *mpi_comm_set_errhandler(PyObject *self, PyObject *args)
{
/* int MPI_Comm_dup ( MPI_Comm comm, MPI_Comm *comm_out ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
long tmpcomm;
COM_TYPE  incomm;
int choice;

	if (!PyArg_ParseTuple(args, "li", &tmpcomm,&choice))
        return NULL;
    incomm=(COM_TYPE)tmpcomm;
    if(choice == 0)
		ierr= MPI_Comm_set_errhandler( (MPI_Comm)incomm, MPI_ERRORS_ARE_FATAL );
    if(choice == 1)
		ierr= MPI_Comm_set_errhandler( (MPI_Comm)incomm, MPI_ERRORS_RETURN );
    if(choice == 2)
		ierr= MPI_Comm_set_errhandler( (MPI_Comm)incomm, newerr );

	return PyInt_FromLong((long)ierr);
}




static PyObject *mpi_comm_create(PyObject *self, PyObject *args)
{
/* int MPI_Comm_create ( MPI_Comm comm, MPI_Group group, MPI_Comm *comm_out ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  incomm,outcomm;
long group;

	if (!PyArg_ParseTuple(args, "ll", &incomm,&group))
        return NULL;
	ierr=MPI_Comm_create ((MPI_Comm)incomm,(MPI_Group)group,(MPI_Comm*)&outcomm );
	return VERT_FUNC((CAST)outcomm);
}


static PyObject *mpi_barrier(PyObject *self, PyObject *args)
{
/* int MPI_Barrier ( MPI_Comm comm ) */
COM_TYPE comm;
int ierr;

	if (!PyArg_ParseTuple(args, "l", &comm))
        return NULL;
    ierr=MPI_Barrier((MPI_Comm)comm );
	return PyInt_FromLong((long)ierr);
}
static PyObject *mpi_send(PyObject *self, PyObject *args)
{
/* int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm ) */
int count,dest,tag;
DATA_TYPE datatype;
COM_TYPE comm;
int i,n;
PyObject *input;
PyArrayObject *array;
char *aptr;

	if (!PyArg_ParseTuple(args, "Oiliil", &input, &count,&datatype,&dest,&tag,&comm))
        return NULL;
	array = (PyArrayObject *) PyArray_ContiguousFromObject(input, getptype(datatype), 0, 3);
	if (array == NULL)
		return NULL;
/*
	n=1;
	for(i=0;i<array->nd;i++)
		n = n* array->dimensions[i];
	if(array->nd == 0)n=1;
	if (n < count)
		return NULL;
*/
	aptr=(char*)(array->data);
    ierr=MPI_Send(aptr,  count,  (MPI_Datatype)datatype,  dest,  tag,  (MPI_Comm)comm );
    Py_DECREF(array);
	return PyInt_FromLong((long)ierr);
}


static PyObject *mpi_recv(PyObject *self, PyObject *args)
{
/* int MPI_Recv( void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status ) */
int count,source,tag;
DATA_TYPE datatype;
COM_TYPE comm;
PyArrayObject *result;
int dimensions[1];
char *aptr;

	if (!PyArg_ParseTuple(args, "iliil", &count,&datatype,&source,&tag,&comm))
        return NULL;
    dimensions[0]=count;
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(datatype));
	aptr=(char*)(result->data);
	ierr=MPI_Recv( aptr,  count, (MPI_Datatype)datatype,source,tag, (MPI_Comm)comm, &status );
  	return PyArray_Return(result);
}

static PyObject *mpi_status(PyObject *self, PyObject *args)
{
/* int MPI_Recv( void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status ) */
PyArrayObject *result;
int dimensions[1],statray[3];

    dimensions[0]=3;
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, PyArray_INT);
	statray[0]=status.MPI_SOURCE;
	statray[1]=status.MPI_TAG;
	statray[2]=status.MPI_ERROR;
	memcpy((void *)(result->data), (void*)statray, (size_t) (12));
  	return PyArray_Return(result);
}

static PyObject *mpi_error(PyObject *self, PyObject *args)
{
	return PyInt_FromLong((long)ierr);
}


static PyObject * mpi_wtime(PyObject *self, PyObject *args) {

        return( PyFloat_FromDouble(MPI_Wtime()));
}

static PyObject * mpi_wtick(PyObject *self, PyObject *args) {

        return( PyFloat_FromDouble(MPI_Wtick()));
}

static PyObject * mpi_attr_get(PyObject *self, PyObject *args) {
/*
  int MPI_Attr_get ( MPI_Comm comm,int keyval,void *attr_value,int *flag)
Input Parameters
        comm    communicator to which attribute is attached (handle)
        keyval  key value (integer)

Output Parameters
        attr_value      attribute value, unless flag = false
        flag    true if an attribute value was extracted; false if no attribute is associated with the key
*/

        int keyval, *attr_value;
        int flag;
        COM_TYPE comm;
        if (!PyArg_ParseTuple(args, "li", &comm, &keyval))
        return NULL;

/*        printf("mpi_attr_get:  keyval:%d\n",keyval); */

        /* get the keyval for the specified attribute */
        ierr = MPI_Attr_get((MPI_Comm)comm, keyval, &attr_value,&flag);
        if ( !flag ) {
                return NULL;
        }

/*        printf("mpi_attr_get:  attr_val: %d  %d\n",*attr_value,flag); */
        return( PyInt_FromLong((long)*attr_value));
}

#ifdef MPI2
static PyObject *mpi_array_of_errcodes(PyObject *self, PyObject *args)
{
	int dimensions[1];
	if(array_of_errcodes){
	dimensions[0]=array_of_errcodes_size;
	return(PyArray_FromDimsAndData(1,
                        dimensions,
                        PyArray_INT,
                        (char *)array_of_errcodes));
	}
	else {
	dimensions[0]=0;
	return(PyArray_FromDimsAndData(1,
                        dimensions,
                        PyArray_INT,
                        (char *)dimensions));
    }

}

static PyObject *mpi_intercomm_merge(PyObject *self, PyObject *args)
{
/* int MPI_Intercomm_merge ( MPI_Comm comm, int high, MPI_Comm *comm_out ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  incomm,outcomm;
int high;

	if (!PyArg_ParseTuple(args, "li", &incomm,&high))
        return NULL;
	ierr=MPI_Intercomm_merge ((MPI_Comm)incomm,high,(MPI_Comm*)&outcomm );
	return VERT_FUNC((CAST)outcomm);
}
static PyObject *mpi_comm_free(PyObject *self, PyObject *args)
{
/* int MPI_Comm_free ( MPI_Comm *comm) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  incomm;


	if (!PyArg_ParseTuple(args, "l", &incomm))
        return NULL;
	ierr=MPI_Comm_free ((MPI_Comm*)&incomm);
	return VERT_FUNC((CAST)ierr);
}

static PyObject *mpi_comm_spawn(PyObject *self, PyObject *args)
{
/* int MPI_Comm_spawn(char *command, char *argv[], int maxprocs, MPI_Info info,
                  int root, MPI_Comm comm, MPI_Comm *intercomm,
                  int array_of_errcodes[])                    */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
int maxprocs,info,root,comm;
MPI_Comm  outcomm;
char *command;
PyObject *input;
int n,len,i;
char **argv;


if (!PyArg_ParseTuple(args, "sOiiii", &command,&input,&maxprocs,&info,&root,&comm))
	return NULL;
	argv=MPI_ARGV_NULL;
	if((input == NULL )|| (input == Py_None )){
		/* printf("is none\n"); */
	}
	if(strncmp("int",input->ob_type->tp_name,3)==0){
		/* printf("is int\n"); */
	}
	if(strncmp("str",input->ob_type->tp_name,3)==0){
		/* printf("is str\n"); */
		argv=(char**)malloc((maxprocs+2)*sizeof(char*));
		for(i=0;i<maxprocs+2;i++) {
			argv[i]=(char*)0;
		}
		n=maxprocs;
		if(n > 0){
			for(i=0;i<n;i++) {
				len=strlen(PyString_AsString(input));
				argv[i]=(char*)malloc(len+1);
				argv[i][len]=(char)0;
				strncpy(argv[i],PyString_AsString(input),(size_t)len);
				/* printf("%s\n",argv[i]); */
			}
		}
	}

	if(strncmp("list",input->ob_type->tp_name,4)==0){
		printf("is list\n");
		argv=(char**)malloc((maxprocs+2)*sizeof(char*));
		for(i=0;i<maxprocs+2;i++) {
			argv[i]=(char*)0;
		}
		n=PyList_Size(input);
		if(n > 0){
			for(i=0;i<n;i++) {
				len=strlen(PyString_AsString(PyList_GetItem(input,i)));
				argv[i]=(char*)malloc(len+1);
				argv[i][len]=(char)0;
				strncpy(argv[i],PyString_AsString(PyList_GetItem(input,i)),(size_t)len);
				printf("%s\n",argv[i]);
			}
		}
	}

	if(array_of_errcodes)free(array_of_errcodes);
	array_of_errcodes=(int*)malloc(maxprocs*sizeof(int));
	array_of_errcodes_size=maxprocs;

/* int MPI_Comm_spawn(char *command, char *argv[], int maxprocs, MPI_Info info,
                  int root, MPI_Comm comm, MPI_Comm *intercomm,
                  int array_of_errcodes[])                    */
/*	printf("launching %s from %d\n",command,root); */
	ierr=MPI_Comm_spawn(command,
	                    argv,
	                    maxprocs,
	                    (MPI_Info)info,
	                    root,
	                    (MPI_Comm)comm,
	                    &outcomm,
	                    array_of_errcodes);


if(argv != MPI_ARGV_NULL){
		n=maxprocs+2;
		for(i=0;i<n;i++) {
			if(argv[i])free(argv[i]);
		}
		free(argv);
		}

return VERT_FUNC((CAST)outcomm);
}

static PyObject *mpi_comm_get_parent(PyObject *self, PyObject *args)
{
/* int MPI_Comm_get_parent(MPI_Comm *parent)*/
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  outcomm;
ierr=MPI_Comm_get_parent((MPI_Comm*)&outcomm);
return VERT_FUNC((CAST)outcomm);
}

static PyObject *mpi_open_port(PyObject *self, PyObject *args)
{
/* int MPI_Open_port ( MPI_Info info, char *port_name) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
char  port_name[MPI_MAX_PORT_NAME];
int info;
	if (!PyArg_ParseTuple(args, "i", &info))
        return NULL;
	ierr=MPI_Open_port(info,port_name);
	return PyString_FromString(port_name);
}

static PyObject *mpi_close_port(PyObject *self, PyObject *args)
{
/* int MPI_Close_port ( char *port_name) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
char  *port_name;
	if (!PyArg_ParseTuple(args, "s", &port_name))
        return NULL;
	ierr=MPI_Close_port(port_name);
	return VERT_FUNC((CAST)ierr);
}

static PyObject *mpi_comm_accept(PyObject *self, PyObject *args)
{
/* int MPI_Comm_accept(char* port_name, MPI_info info,int root, MPI_Comm comm, MPI_Comm *newcomm)*/
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  comm,newcomm;
char* port_name;
int info,root;
	if (!PyArg_ParseTuple(args, "siii", &port_name,&info,&root,&comm))
        return NULL;
	ierr=MPI_Comm_accept(port_name,  info, root,  comm,  (MPI_Comm*)&newcomm);
	return VERT_FUNC((CAST)newcomm);
}
static PyObject *mpi_comm_connect(PyObject *self, PyObject *args)
{
/* int MPI_Comm_connect(char* port_name, MPI_info info,int root, MPI_Comm comm, MPI_Comm *newcomm)*/
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  comm,newcomm;
char* port_name;
int info,root;
	if (!PyArg_ParseTuple(args, "siii", &port_name,&info,&root,&comm))
        return NULL;
	ierr=MPI_Comm_connect(port_name,  info, root,  comm,  (MPI_Comm*)&newcomm);
	return VERT_FUNC((CAST)newcomm);
}
static PyObject *mpi_comm_disconnect(PyObject *self, PyObject *args)
{
/* int MPI_Comm_disconnect(MPI_Comm *comm)*/
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
COM_TYPE  comm;
	if (!PyArg_ParseTuple(args, "i", &comm))
        return NULL;
	ierr=MPI_Comm_disconnect((MPI_Comm*)&comm);
	return VERT_FUNC((CAST)ierr);
}

#endif

static PyObject *mpi_comm_split(PyObject *self, PyObject *args)
{
/* int MPI_Comm_split ( MPI_Comm comm, int color, int key, MPI_Comm *comm_out ) */
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
int color,key;
COM_TYPE  incomm,outcomm;

	if (!PyArg_ParseTuple(args, "lii", &incomm,&color,&key))
        return NULL;
	ierr=MPI_Comm_split ((MPI_Comm)incomm,color,key,(MPI_Comm*)&outcomm );
	return VERT_FUNC((CAST)outcomm);
}

static PyObject *mpi_probe(PyObject *self, PyObject *args)
{
/* int MPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status ) */
int source,tag;
COM_TYPE comm;

	if (!PyArg_ParseTuple(args, "iil", &source,&tag,&comm))
        return NULL;
    ierr=MPI_Probe(source,tag, (MPI_Comm )comm, &status );
	return PyInt_FromLong((long)0);
}

static PyObject *mpi_get_count(PyObject *self, PyObject *args)
{
/* int MPI_Get_count( MPI_Status *status, MPI_Datatype datatype, int *count ) */
DATA_TYPE datatype;
int count;

	if (!PyArg_ParseTuple(args, "l",&datatype))
        return NULL;
    ierr=MPI_Get_count(&status,(MPI_Datatype)datatype,&count);
	return PyInt_FromLong((long)count);
}


static PyObject *mpi_comm_size(PyObject *self, PyObject *args)
{
/* int MPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status ) */
COM_TYPE comm;
int numprocs;

	if (!PyArg_ParseTuple(args, "l",&comm))
        return NULL;
	ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
	return PyInt_FromLong((long)numprocs);
}

static PyObject *mpi_comm_rank(PyObject *self, PyObject *args)
{
/* int MPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status ) */
COM_TYPE comm;
int rank;

	if (!PyArg_ParseTuple(args, "l",&comm))
        return NULL;
	ierr=MPI_Comm_rank((MPI_Comm)comm,&rank);
	return PyInt_FromLong((long)rank);
}


static PyObject *mpi_iprobe(PyObject *self, PyObject *args)
{
/* int MPI_Iprobe( int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status ) */
int source,tag,flag;
COM_TYPE comm;

	if (!PyArg_ParseTuple(args, "iil", &source,&tag,&comm))
        return NULL;
    ierr=MPI_Iprobe(source,tag, (MPI_Comm )comm, &flag,&status );
	return PyInt_FromLong((long)flag);
}

static PyObject * mpi_init(PyObject *self, PyObject *args) {
	PyObject *input;
	int argc,i,n;
	int did_it;
	char **argv;
#ifdef ARG_STR
	char *argstr;
	int *strides;
	int arglen;
#endif
#ifdef SIZE_RANK
	PyArrayObject *result;
	int dimensions[1],data[2];
	char *aptr;
#endif
#ifdef ARG_ARRAY
	PyObject *result;
#endif
    int myid;
	int len;
	argv=NULL;
        erroron=0;
	ierr=MPI_Initialized(&did_it);
	if(!did_it){
		if (!PyArg_ParseTuple(args, "iO", &argc, &input))
			return NULL;
		argv=(char**)malloc((argc+2)*sizeof(char*));
		n=PyList_Size(input);
		for(i=0;i<n;i++) {
			len=strlen(PyString_AsString(PyList_GetItem(input,i)));
			argv[i]=(char*)malloc(len+1);
			argv[i][len]=(char)0;
			strncpy(argv[i],PyString_AsString(PyList_GetItem(input,i)),(size_t)len);
			/* printf("%s ",argv[i]); */
		}

		/* printf("\n"); */
		Py_DECREF(input);
#ifdef LAM_MPI
		ierr=MPI_Init(NULL,NULL);
#else
		ierr=MPI_Init(&argc,&argv);
		ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid);

/*
for(i=0;i<argc;i++) {
		   printf("from %d init %d %s\n",myid,i,argv[i]);

		}
*/
#endif

#ifdef MPI2
		MPI_Comm_create_errhandler( eh, &newerr );
#endif

/*		free(argv); */
	}

/*  this returns the command line as a string */
#ifdef ARG_STR
		arglen=0;
		strides=(int*)malloc(argc*sizeof(int));
		strides[0]=0;
		for(i=0;i<argc;i++) {
			arglen=arglen+strlen(argv[i])+1;
			strides[i+1]=strides[i]+strlen(argv[i])+1;
		}
		argstr=(char*)malloc(arglen*sizeof(char));
		for(i=0;i<argc;i++) {
		    for(n=0;n<strlen(argv[i]);n++) {
		        argstr[strides[i]+n]=argv[i][n];
		    }
		    argstr[strides[i]+strlen(argv[i])]=(char)32;
/*
			free(argv[i]);
*/
		}
		return PyString_FromString(argstr);
#endif
#ifdef ARG_ARRAY
		result = PyTuple_New(argc);
		for(i=0;i<argc;i++) {
			PyTuple_SetItem(result,i,PyString_FromString(argv[i]));
		}
/*
for(i=0;i<argc;i++) {
			free(argv[i]);
		}
*/
		return result;
#endif
/*  this returns size and rank */
#ifdef SIZE_RANK
    ierr=MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	dimensions[0]=2;
	result = (PyArrayObject *)PyArray_FromDims(1, dimensions, PyArray_INT);
	if (result == NULL)
		return NULL;
	data[0]=myid;
	data[1]=numprocs;
	aptr=(char*)&(data);
	for(i=0;i<8;i++)
		result->data[i]=aptr[i];
	if(erroron){ erroron=0; return NULL;}
	return PyArray_Return(result);
#endif
}

static PyObject * mpi_start(PyObject *self, PyObject *args) {
	PyArrayObject *result;
	int argc,did_it,i;
	int dimensions[1],data[2];
	int numprocs,myid;
	char *command,*aptr;
	char **argv;
        erroron=0;
	if (!PyArg_ParseTuple(args, "is", &argc, &command))
        return NULL;
	ierr=MPI_Initialized(&did_it);
    if(!did_it){
		/* MPI_Init(0,0); */ /* lam mpi will start with this line
		                        mpich requires us to build a real
		                        command line */
		/* MPI_Init(0,0); */
		argv=(char**)malloc((argc+2)*sizeof(char*));
		argv[0]=(char*)malloc(128);
		sprintf(argv[0],"dummy");
		strtok(command, " ");
		for(i=0;i<argc-1;i++) {
				argv[i+1]=strtok(NULL, " ");
		}
		printf("calling mpi init from mpi_start\n");
		for (i=0;i<argc;i++)
			printf("%d %d %s\n",i,(int)strlen(argv[i]),argv[i]);
#ifdef LAM_MPI
		ierr=MPI_Init(NULL,NULL);
#else
		ierr=MPI_Init(&argc,&argv);
#endif
	}
	ierr=MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    ierr=MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#ifdef MPI2
		MPI_Comm_create_errhandler( eh, &newerr );
#endif


	dimensions[0]=2;
	result = (PyArrayObject *)PyArray_FromDims(1, dimensions, PyArray_INT);
	if (result == NULL)
		return NULL;
	data[0]=myid;
	data[1]=numprocs;
	aptr=(char*)&(data);
	for(i=0;i<8;i++)
		result->data[i]=aptr[i];
	if(erroron){ erroron=0; return NULL;}
	return PyArray_Return(result);
}

static PyObject * mpi_bcast(PyObject *self, PyObject *args) {
/* int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ) */
int count,root;
DATA_TYPE datatype;
COM_TYPE comm;
int myid;
int mysize;
PyArrayObject *result;
PyArrayObject *array;
PyObject *input;
int dimensions[1];
char *aptr;

	if (!PyArg_ParseTuple(args, "Oilil", &input, &count,&datatype,&root,&comm))
        return NULL;
    dimensions[0]=count;
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(datatype));
	aptr=(char*)(result->data);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
		array = (PyArrayObject *) PyArray_ContiguousFromObject(input, getptype(datatype), 0, 3);
		if (array == NULL)
			return NULL;
		ierr=MPI_Type_size((MPI_Datatype)datatype,&mysize);
		memcpy((void *)(result->data), (void*)array->data, (size_t) (mysize*count));
		Py_DECREF(array);
	}
	ierr=MPI_Bcast(aptr,count,(MPI_Datatype)datatype,root,(MPI_Comm)comm);
  	return PyArray_Return(result);
}


static PyObject * mpi_scatterv(PyObject *self, PyObject *args) {
/* int MPI_Scatterv(void *sendbuf, int *sendcnts, int *displs, MPI_Datatype sendtype,
                    void *recvbuf, int recvcnt,                MPI_Datatype recvtype,
                    int root, MPI_Comm comm ) */
int root;
COM_TYPE comm;
DATA_TYPE sendtype,recvtype;
PyObject *sendbuf_obj, *sendcnts_obj,*displs_obj;
PyArrayObject *array,*result;
int *sendcnts,*displs,recvcnt;
char *sptr,*rptr;
int *sray;
int numprocs,myid;
int dimensions[1];

	sendcnts=0;
	displs=0;

	array=NULL;
	sptr=NULL;

	if (!PyArg_ParseTuple(args, "OOOlilil", &sendbuf_obj, &sendcnts_obj,&displs_obj,&sendtype,&recvcnt,&recvtype,&root,&comm))
        return NULL;
    /* ! get the number of processors in this comm */
    ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);


#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendcnts_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		sendcnts=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)sendcnts, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		Py_DECREF(array);
		array = (PyArrayObject *) PyArray_ContiguousFromObject(displs_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		displs=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)displs, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		Py_DECREF(array);
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendbuf_obj, getptype(sendtype), 1, 3);
		if (array == NULL)
			return NULL;
		sptr=(char*)(array->data);
	}

    dimensions[0]=recvcnt;
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(recvtype));
    rptr=(char*)(result->data);

	ierr=MPI_Scatterv(sptr, sendcnts, displs, (MPI_Datatype)sendtype,rptr,recvcnt,(MPI_Datatype)recvtype, root, (MPI_Comm )comm );
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
		Py_DECREF(array);
		free(sendcnts);
		free(displs);
	}
  	return PyArray_Return(result);
}

static PyObject * mpi_gatherv(PyObject *self, PyObject *args) {
/*
int MPI_Gatherv ( void *sendbuf, int sendcnt,                MPI_Datatype sendtype,
                  void *recvbuf, int *recvcnts, int *displs, MPI_Datatype recvtype,
                  int root, MPI_Comm comm )
 */
int root;
COM_TYPE comm;
DATA_TYPE sendtype,recvtype;
PyObject *sendbuf_obj, *recvcnts_obj,*displs_obj;
PyArrayObject *array,*result;
int sendcnt,*displs,*recvcnts,rtot,i;
char *sptr,*rptr;
int numprocs,myid;
int dimensions[1];

	displs=0;

	array=NULL;
	sptr=NULL;

	if (!PyArg_ParseTuple(args, "OilOOlil", &sendbuf_obj, &sendcnt,&sendtype,&recvcnts_obj,&displs_obj,&recvtype,&root,&comm))
        return NULL;
    /* ! get the number of processors in this comm */
    ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
    rtot=0;
    recvcnts=0;
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
    /* printf("  get the recv_counts array \n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(recvcnts_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		recvcnts=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)recvcnts, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		rtot=0;
		for(i=0;i<numprocs;i++)
			rtot=rtot+recvcnts[i];
		Py_DECREF(array);
    /* printf("  get the offset array \n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(displs_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		displs=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)displs, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		Py_DECREF(array);
	}
	/* printf("  allocate the recvbuf \n"); */
		dimensions[0]=rtot;
		result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(recvtype));
		rptr=(char*)(result->data);
   /* printf("  get sendbuf\n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendbuf_obj, getptype(sendtype), 1, 3);
		if (array == NULL)
			return NULL;
		sptr=array->data;


   /* printf("   do the call %d \n",recvcnt); */
	ierr=MPI_Gatherv(sptr, sendcnt, (MPI_Datatype)sendtype,rptr,recvcnts,displs,(MPI_Datatype)recvtype, root, (MPI_Comm )comm );
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
		Py_DECREF(array);
		free(recvcnts);
		free(displs);
	}
   /* printf("   did the call  %d \n",myid); */
  	return PyArray_Return(result);
}

static PyObject * mpi_gather(PyObject *self, PyObject *args) {
/*
int MPI_Gather ( void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                  void *recvbuf, int recvcnts,
                  MPI_Datatype recvtype,
                  int root, MPI_Comm comm )
 */
int root;
COM_TYPE comm;
DATA_TYPE sendtype,recvtype;
PyObject *sendbuf_obj;
PyArrayObject *array,*result;
int sendcnt,recvcnt,rtot;
char *sptr,*rptr;
int numprocs,myid;
int dimensions[1];

	array=NULL;
	sptr=NULL;

	if (!PyArg_ParseTuple(args, "Oililil", &sendbuf_obj, &sendcnt,&sendtype,&recvcnt,&recvtype,&root,&comm))
        return NULL;
    /* ! get the number of processors in this comm */
    ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
    rtot=0;
   /* printf("  get sendbuf\n"); */
	array = (PyArrayObject *) PyArray_ContiguousFromObject(sendbuf_obj, getptype(sendtype), 0, 3);
	if (array == NULL)
		return NULL;
	sptr=array->data;
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
		rtot=recvcnt*numprocs;
    }
	/* printf("  allocate the recvbuf \n"); */
	dimensions[0]=rtot;
	result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(recvtype));
	rptr=(char*)(result->data);


   /* printf("   do the call %d \n",recvcnt); */
	ierr=MPI_Gather(sptr, sendcnt, (MPI_Datatype)sendtype,rptr,recvcnt,(MPI_Datatype)recvtype, root, (MPI_Comm )comm );
	Py_DECREF(array);
   /* printf("   did the call  %d \n",myid); */
  	return PyArray_Return(result);
}

static PyObject * mpi_scatter(PyObject *self, PyObject *args) {
/*
  int MPI_Scatter ( void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                    void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                    int root, MPI_Comm comm )
*/
int root;
COM_TYPE comm;
DATA_TYPE sendtype,recvtype;
PyObject *sendbuf_obj;
PyArrayObject *array,*result;
int sendcnts,recvcnt;
int numprocs,myid;
int dimensions[1];
char *sptr,*rptr;

	sendcnts=0;

	array=NULL;
	sptr=NULL;

	if (!PyArg_ParseTuple(args, "Oililil", &sendbuf_obj, &sendcnts,&sendtype,&recvcnt,&recvtype,&root,&comm))
        return NULL;
    /* ! get the number of processors in this comm */
    ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);

#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
    /* get sendbuf */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendbuf_obj, getptype(sendtype), 1, 3);
		if (array == NULL)
			return NULL;
		    sptr=(char*)(array->data);

    }

    /* allocate the recvbuf */
    dimensions[0]=recvcnt;
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(recvtype));
    rptr=(char*)(result->data);

   /*  do the call */
	ierr=MPI_Scatter(sptr, sendcnts, (MPI_Datatype)sendtype,rptr,recvcnt,(MPI_Datatype)recvtype, root, (MPI_Comm )comm );
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
		Py_DECREF(array);
	}
  	return PyArray_Return(result);
}

static PyObject * mpi_reduce(PyObject *self, PyObject *args) {
/* int MPI_Reduce ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm ) */
int count,root;
DATA_TYPE datatype;
COM_TYPE comm;
long op;
int myid;
PyArrayObject *result;
PyArrayObject *array;
PyObject *input;
int dimensions[1];
char *sptr,*rptr;

	if (!PyArg_ParseTuple(args, "Oillil", &input, &count,&datatype,&op,&root,&comm))
        return NULL;
    MPI_Comm_rank((MPI_Comm)comm,&myid);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);
#ifdef MPI2
    if(myid == root || root == MPI_ROOT) {
#else
    if(myid == root) {
#endif
    	dimensions[0]=count;
	}
	else {
		dimensions[0]=0;
	}
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(datatype));
	rptr=(char*)(result->data);

	array = (PyArrayObject *) PyArray_ContiguousFromObject(input, getptype(datatype), 0, 3);
	if (array == NULL)
		return NULL;
	sptr=(char*)(array->data);
	ierr=MPI_Reduce(sptr,rptr,count,(MPI_Datatype)datatype,(MPI_Op)op,root,(MPI_Comm)comm);
	Py_DECREF(array);
  	return PyArray_Return(result);
}


static PyObject * mpi_finalize(PyObject *self, PyObject *args) {
/* int MPI_Finalize() */
	if(erroron){ erroron=0; return NULL;}
    return PyInt_FromLong((long)MPI_Finalize());
}
static PyObject * mpi_alltoall(PyObject *self, PyObject *args) {
/*
   int MPI_Alltoall( void *sendbuf, int sendcount, MPI_Datatype sendtype,
                     void *recvbuf, int recvcnt,   MPI_Datatype recvtype,
                     MPI_Comm comm )
 */
COM_TYPE comm;
DATA_TYPE sendtype,recvtype;
PyObject *sendbuf_obj;
PyArrayObject *array,*result;
int sendcnts,recvcnt;
int numprocs,myid;
int dimensions[1];
char *sptr,*rptr;
sendcnts=0;

	array=NULL;
	sptr=NULL;

	if (!PyArg_ParseTuple(args, "Oilill", &sendbuf_obj, &sendcnts,&sendtype,&recvcnt,&recvtype,&comm))
        return NULL;
    /* ! get the number of processors in this comm */
    ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
    ierr=MPI_Comm_rank((MPI_Comm)comm,&myid);

    /* get sendbuf */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendbuf_obj, getptype(sendtype), 1, 3);
		if (array == NULL)
			return NULL;
		    sptr=(char*)(array->data);


    /* allocate the recvbuf */
    dimensions[0]=recvcnt*numprocs;
    result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(recvtype));
    rptr=(char*)(result->data);

   /*  do the call */
	ierr=MPI_Alltoall(sptr, sendcnts, (MPI_Datatype)sendtype,rptr,recvcnt,(MPI_Datatype)recvtype, (MPI_Comm )comm );
	Py_DECREF(array);
  	return PyArray_Return(result);
}
static PyObject * mpi_alltoallv(PyObject *self, PyObject *args) {
/*
  int MPI_Alltoallv ( void *sendbuf, int *sendcnts, int *sdispls, MPI_Datatype sendtype,
                      void *recvbuf, int *recvcnts, int *rdispls, MPI_Datatype recvtype,
                      MPI_Comm comm )
*/
COM_TYPE comm;
DATA_TYPE sendtype,recvtype;
PyObject *sendbuf_obj, *recvcnts_obj,*rdispls_obj,*sdispls_obj,*sendcnts_obj;
PyArrayObject *array,*result;
int *sendcnts,*sdispls,*rdispls,*recvcnts,rtot,i;
char *sptr,*rptr;
int numprocs;
int dimensions[1];

	rdispls=0;

	array=NULL;
	sptr=NULL;

	if (!PyArg_ParseTuple(args, "OOOlOOll", &sendbuf_obj, &sendcnts_obj,&sdispls_obj,&sendtype,&recvcnts_obj,&rdispls_obj,&recvtype,&comm))
        return NULL;
    /* ! get the number of processors in this comm */
    ierr=MPI_Comm_size((MPI_Comm)comm,&numprocs);
    rtot=0;
    recvcnts=0;

    /* printf("  get the recvcnts array \n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(recvcnts_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		recvcnts=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)recvcnts, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		rtot=0;
		for(i=0;i<numprocs;i++)
			rtot=rtot+recvcnts[i];
		Py_DECREF(array);

    /* printf("  get the recv offset array \n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(rdispls_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		rdispls=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)rdispls, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		Py_DECREF(array);

	/* printf("  allocate the recvbuf \n"); */
		dimensions[0]=rtot;
		result = (PyArrayObject *)PyArray_FromDims(1, dimensions, getptype(recvtype));
		rptr=(char*)(result->data);



    /* printf("  get the sendcnts array \n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendcnts_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		sendcnts=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)sendcnts, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		Py_DECREF(array);

    /* printf("  get the send offset array \n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sdispls_obj, PyArray_INT, 1, 1);
		if (array == NULL)
			return NULL;
		sdispls=(int*)malloc((size_t) (sizeof(int)*numprocs));
		memcpy((void *)sdispls, (void*)array->data, (size_t) (sizeof(int)*numprocs));
		Py_DECREF(array);

   /* printf("  get sendbuf\n"); */
		array = (PyArrayObject *) PyArray_ContiguousFromObject(sendbuf_obj, getptype(sendtype), 1, 3);
		if (array == NULL)
			return NULL;
		sptr=(char*)array->data;

   /* printf("   do the call %d \n"); */
   /*
        MPI_Comm_rank((MPI_Comm)comm,&myid);
   		printf("myid =%d ",myid);
   		for(i=0;i<numprocs;i++)
   			printf("%d ",sendcnts[i]);
   		printf(" | ");
   		for(i=0;i<numprocs;i++)
   			printf("%d ",sdispls[i]);
   		printf(" | ");
   		for(i=0;i<numprocs;i++)
   			printf("%d ",recvcnts[i]);
   		printf(" | ");
   		for(i=0;i<numprocs;i++)
   			printf("%d ",rdispls[i]);
   		printf("\n");
   */
       ierr=MPI_Alltoallv(sptr, sendcnts, sdispls, (MPI_Datatype)sendtype,
                          rptr, recvcnts, rdispls, (MPI_Datatype)recvtype,
                          (MPI_Comm)comm);

		Py_DECREF(array);
		free(recvcnts);
		free(rdispls);
		free(sendcnts);
		free(sdispls);

  	return PyArray_Return(result);
}
#include <string.h>
static PyObject * copywrite(PyObject *self, PyObject *args) {
int i;
for(i=0;i<2160;i++) {
	cw[i]=(char)0;
}
strncat(cw,"Copyright (c) 2005 The Regents of the University of California\n",80);
strncat(cw,"All Rights Reserved\n",80);
strncat(cw,"Permission to use, copy, modify and distribute any part of this\n",80);
strncat(cw,"	software for educational, research and non-profit purposes,\n",80);
strncat(cw,"	without fee, and without a written agreement is hereby granted,\n",80);
strncat(cw,"	provided that the above copyright notice, this paragraph and the\n",80);
strncat(cw,"	following three paragraphs appear in all copies.\n",80);
strncat(cw,"Those desiring to incorporate this software into commercial products or\n",80);
strncat(cw,"	use for commercial purposes should contact the Technology\n",80);
strncat(cw,"	Transfer & Intellectual Property Services, University of\n",80);
strncat(cw,"	California, San Diego, 9500 Gilman Drive, Mail Code 0910, La\n",80);
strncat(cw,"	Jolla, CA 92093-0910, Ph: (858) 534-5815, FAX: (858) 534-7345,\n",80);
strncat(cw,"	E-MAIL:invent@ucsd.edu.\n",80);
strncat(cw,"IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY\n",80);
strncat(cw,"	PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR\n",80);
strncat(cw,"	CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF\n",80);
strncat(cw,"	THE USE OF THIS SOFTWARE  EVEN IF THE UNIVERSITY OF\n",80);
strncat(cw,"	CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n",80);
strncat(cw,"THE SOFTWARE PROVIDED HEREIN IS ON AN AS IS BASIS, AND THE\n",80);
strncat(cw,"	UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE\n",80);
strncat(cw,"	MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.\n",80);
strncat(cw,"	THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND\n",80);
strncat(cw,"	EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR EXPRESS,\n",80);
strncat(cw,"	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF\n",80);
strncat(cw,"	MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE\n",80);
strncat(cw,"	USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR\n",80);
strncat(cw,"	OTHER RIGHTS.\n",80);
return PyString_FromString(cw);
}
/*
**** to add ****
mpi_isend
mpi_irecv
mpi_test
mpi_wait
*/


static PyMethodDef mpiMethods[] = {
    {"mpi_alltoall",	mpi_alltoall,		METH_VARARGS,     	 mpi_alltoall__},
    {"mpi_alltoallv",	mpi_alltoallv,      METH_VARARGS,     	 mpi_alltoallv__},
    {"mpi_barrier",		mpi_barrier,		METH_VARARGS,     	 mpi_barrier__},
    {"mpi_bcast",		mpi_bcast,			METH_VARARGS,     	 mpi_bcast__},
    {"mpi_comm_create",	mpi_comm_create,	METH_VARARGS,     	 mpi_comm_create__},
    {"mpi_comm_dup",	mpi_comm_dup,		METH_VARARGS,     	 mpi_comm_dup__},
    {"mpi_comm_group",	mpi_comm_group,		METH_VARARGS,     	 mpi_comm_group__},
    {"mpi_comm_rank",	mpi_comm_rank,		METH_VARARGS,     	 mpi_comm_rank__},
    {"mpi_comm_size",	mpi_comm_size,		METH_VARARGS,     	 mpi_comm_size__},
    {"mpi_comm_split",	mpi_comm_split,		METH_VARARGS,     	 mpi_comm_split__},
    {"mpi_error",		mpi_error,			METH_VARARGS,     	 mpi_error__},
    {"mpi_finalize",	mpi_finalize,		METH_VARARGS,     	 mpi_finalize__},
    {"mpi_gather",		mpi_gather,			METH_VARARGS,     	 mpi_gather__},
    {"mpi_gatherv",		mpi_gatherv,		METH_VARARGS,     	 mpi_gatherv__},
    {"mpi_get_count",	mpi_get_count,		METH_VARARGS,     	 mpi_get_count__},
    {"mpi_group_incl",	mpi_group_incl,		METH_VARARGS,     	 mpi_group_incl__},
    {"mpi_group_rank",	mpi_group_rank,		METH_VARARGS,     	 mpi_group_rank__},
    {"mpi_init",		mpi_init,			METH_VARARGS,     	 mpi_init__},
    {"mpi_iprobe",		mpi_iprobe,			METH_VARARGS,     	 mpi_iprobe__},
    {"mpi_irecv",		mpi_irecv,			METH_VARARGS,     	 mpi_irecv__},
    {"mpi_isend",		mpi_isend,			METH_VARARGS,     	 mpi_isend__},
    {"mpi_probe",		mpi_probe,			METH_VARARGS,     	 mpi_probe__},
    {"mpi_recv",		mpi_recv,			METH_VARARGS,     	 mpi_recv__},
    {"mpi_reduce",		mpi_reduce,			METH_VARARGS,     	 mpi_reduce__},
    {"mpi_scatter",		mpi_scatter,		METH_VARARGS,     	 mpi_scatter__},
    {"mpi_scatterv",	mpi_scatterv,		METH_VARARGS,     	 mpi_scatterv__},
    {"mpi_send",		mpi_send,			METH_VARARGS,     	 mpi_send__},
    {"mpi_start",		mpi_start,			METH_VARARGS,     	 mpi_start__},
    {"mpi_status",		mpi_status,			METH_VARARGS,     	 mpi_status__},
    {"mpi_get_processor_name",		mpi_get_processor_name,			METH_VARARGS,     	 mpi_get_processor_name__},
    {"mpi_test",		mpi_test,			METH_VARARGS,     	 mpi_test__},
    {"mpi_wait",		mpi_wait,			METH_VARARGS,     	 mpi_wait__},
    {"mpi_wtime",		mpi_wtime,			METH_VARARGS,     	 mpi_wtime__},
    {"mpi_wtick",		mpi_wtick,			METH_VARARGS,     	 mpi_wtick__},
    {"mpi_attr_get",	mpi_attr_get,		METH_VARARGS,     	 mpi_attr_get__},


#ifdef MPI2
    {"mpi_comm_spawn",		    mpi_comm_spawn,			METH_VARARGS,     	 mpi_comm_spawn__},
    {"mpi_array_of_errcodes",	mpi_array_of_errcodes,	METH_VARARGS,     	 mpi_array_of_errcodes__},
    {"mpi_comm_get_parent",		mpi_comm_get_parent,	METH_VARARGS,     	 mpi_comm_get_parent__},
    {"mpi_comm_free",			mpi_comm_free,			METH_VARARGS,     	 mpi_comm_free__},
    {"mpi_intercomm_merge",		mpi_intercomm_merge,	METH_VARARGS,     	 mpi_intercomm_merge__},
    {"mpi_open_port",			mpi_open_port,			METH_VARARGS,     	 mpi_open_port__},
    {"mpi_close_port",			mpi_close_port,			METH_VARARGS,     	 mpi_close_port__},
    {"mpi_comm_accept",			mpi_comm_accept,		METH_VARARGS,     	 mpi_comm_accept__},
    {"mpi_comm_connect",		mpi_comm_connect,		METH_VARARGS,     	 mpi_comm_connect__},
    {"mpi_comm_disconnect",		mpi_comm_disconnect,	METH_VARARGS,     	 mpi_comm_disconnect__},
    {"mpi_comm_set_errhandler",	mpi_comm_set_errhandler,METH_VARARGS,     	 mpi_comm_set_errhandler__},

#endif
    {"copywrite",		copywrite,			METH_VARARGS,     	COPYWRITE_STR__},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

/*
*/

PyMODINIT_FUNC initmpi(void)
{
#ifdef DO_UNSIGED
#define CAST unsigned long
#define VERT_FUNC PyLong_FromUnsignedLong
#else
#define CAST long
#define VERT_FUNC PyInt_FromLong
#endif
	PyObject *m, *d;
    PyObject *tmp;
	import_array();
	m=Py_InitModule("mpi", mpiMethods);
	mpiError = PyErr_NewException("mpi.error", NULL, NULL);
    Py_INCREF(mpiError);
    PyModule_AddObject(m, "error", mpiError);
    d = PyModule_GetDict(m);
    tmp = PyString_FromString(VERSION);
    PyDict_SetItemString(d,   "VERSION", tmp);  Py_DECREF(tmp);
    tmp = PyInt_FromLong((long)MPI_VERSION);
    PyDict_SetItemString(d,   "MPI_VERSION", tmp);  Py_DECREF(tmp);
    tmp = PyInt_FromLong((long)MPI_SUBVERSION);
    PyDict_SetItemString(d,   "MPI_SUBVERSION", tmp);  Py_DECREF(tmp);
    tmp = PyString_FromString(COPYWRITE);
    PyDict_SetItemString(d,   "COPYWRITE", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_CHAR);
    PyDict_SetItemString(d,   "MPI_CHAR", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_BYTE);
    PyDict_SetItemString(d,   "MPI_BYTE", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_SHORT);
    PyDict_SetItemString(d,   "MPI_SHORT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_INT);
    PyDict_SetItemString(d,   "MPI_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LONG);
    PyDict_SetItemString(d,   "MPI_LONG", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_FLOAT);
    PyDict_SetItemString(d,   "MPI_FLOAT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_DOUBLE);
    PyDict_SetItemString(d,   "MPI_DOUBLE", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_UNSIGNED_CHAR);
    PyDict_SetItemString(d,   "MPI_UNSIGNED_CHAR", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_UNSIGNED_SHORT);
    PyDict_SetItemString(d,   "MPI_UNSIGNED_SHORT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_UNSIGNED);
    PyDict_SetItemString(d,   "MPI_UNSIGNED", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_UNSIGNED_LONG);
    PyDict_SetItemString(d,   "MPI_UNSIGNED_LONG", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LONG_DOUBLE);
    PyDict_SetItemString(d,   "MPI_LONG_DOUBLE", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_FLOAT_INT);
    PyDict_SetItemString(d,   "MPI_FLOAT_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LONG_INT);
    PyDict_SetItemString(d,   "MPI_LONG_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_DOUBLE_INT);
    PyDict_SetItemString(d,   "MPI_DOUBLE_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_SHORT_INT);
    PyDict_SetItemString(d,   "MPI_SHORT_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_2INT);
    PyDict_SetItemString(d,   "MPI_2INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LONG_DOUBLE_INT);
    PyDict_SetItemString(d,   "MPI_LONG_DOUBLE_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LONG_LONG_INT);
    PyDict_SetItemString(d,   "MPI_LONG_LONG_INT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_PACKED);
    PyDict_SetItemString(d,   "MPI_PACKED", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_Pack);
    PyDict_SetItemString(d,   "MPI_Pack", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_Unpack);
    PyDict_SetItemString(d,   "MPI_Unpack", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_UB);
    PyDict_SetItemString(d,   "MPI_UB", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LB);
    PyDict_SetItemString(d,   "MPI_LB", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MAX);
    PyDict_SetItemString(d,   "MPI_MAX", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MIN);
    PyDict_SetItemString(d,   "MPI_MIN", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_SUM);
    PyDict_SetItemString(d,   "MPI_SUM", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_PROD);
    PyDict_SetItemString(d,   "MPI_PROD", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LAND);
    PyDict_SetItemString(d,   "MPI_LAND", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_BAND);
    PyDict_SetItemString(d,   "MPI_BAND", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LOR);
    PyDict_SetItemString(d,   "MPI_LOR", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_BOR);
    PyDict_SetItemString(d,   "MPI_BOR", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_LXOR);
    PyDict_SetItemString(d,   "MPI_LXOR", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_BXOR);
    PyDict_SetItemString(d,   "MPI_BXOR", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MINLOC);
    PyDict_SetItemString(d,   "MPI_MINLOC", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MAXLOC);
    PyDict_SetItemString(d,   "MPI_MAXLOC", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MAXLOC);
    PyDict_SetItemString(d,   "MPI_MAXLOC", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_COMM_NULL);
    PyDict_SetItemString(d,   "MPI_COMM_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_OP_NULL);
    PyDict_SetItemString(d,   "MPI_OP_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_GROUP_NULL);
    PyDict_SetItemString(d,   "MPI_GROUP_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_DATATYPE_NULL);
    PyDict_SetItemString(d,   "MPI_DATATYPE_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_REQUEST_NULL);
    PyDict_SetItemString(d,   "MPI_REQUEST_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_ERRHANDLER_NULL);
    PyDict_SetItemString(d,   "MPI_ERRHANDLER_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MAX_PROCESSOR_NAME);
    PyDict_SetItemString(d,   "MPI_MAX_PROCESSOR_NAME", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_MAX_ERROR_STRING);
    PyDict_SetItemString(d,   "MPI_MAX_ERROR_STRING", tmp);  Py_DECREF(tmp);
	tmp = PyInt_FromLong((long)MPI_UNDEFINED);
    PyDict_SetItemString(d,   "MPI_UNDEFINED", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_KEYVAL_INVALID);
    PyDict_SetItemString(d,   "MPI_KEYVAL_INVALID", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_BSEND_OVERHEAD);
    PyDict_SetItemString(d,   "MPI_BSEND_OVERHEAD", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_PROC_NULL);
    PyDict_SetItemString(d,   "MPI_PROC_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_ANY_SOURCE);
    PyDict_SetItemString(d,   "MPI_ANY_SOURCE", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_ANY_TAG);
    PyDict_SetItemString(d,   "MPI_ANY_TAG", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_BOTTOM);
    PyDict_SetItemString(d,   "MPI_BOTTOM", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_COMM_WORLD);
    PyDict_SetItemString(d,   "MPI_COMM_WORLD", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_TAG_UB);
    PyDict_SetItemString(d,   "MPI_TAG_UB", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_HOST);
    PyDict_SetItemString(d,   "MPI_HOST", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_IO);
    PyDict_SetItemString(d,   "MPI_IO", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_WTIME_IS_GLOBAL);
    PyDict_SetItemString(d,   "MPI_WTIME_IS_GLOBAL", tmp);  Py_DECREF(tmp);

    tmp = PyString_FromString(LIBRARY);
    PyDict_SetItemString(d,   "ARRAY_LIBRARY", tmp);  Py_DECREF(tmp);

    tmp = PyString_FromString(DATE_DOC);
    PyDict_SetItemString(d,   "DATE_DOC", tmp);  Py_DECREF(tmp);
    tmp = PyString_FromString(DATE_SRC);
    PyDict_SetItemString(d,   "DATE_SRC", tmp);  Py_DECREF(tmp);

#ifdef MPI2
    tmp = VERT_FUNC((CAST)MPI_UNIVERSE_SIZE);
    PyDict_SetItemString(d,   "MPI_UNIVERSE_SIZE", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_ROOT);
    PyDict_SetItemString(d,   "MPI_ROOT", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_ARGV_NULL);
    PyDict_SetItemString(d,   "MPI_ARGV_NULL", tmp);  Py_DECREF(tmp);
    tmp = VERT_FUNC((CAST)MPI_INFO_NULL);
    PyDict_SetItemString(d,   "MPI_INFO_NULL", tmp);  Py_DECREF(tmp);
#endif
}
void myerror(char *s) {
	erroron=1;
	PyErr_SetString(mpiError,s);
}

/*int MPIR_Err_return_comm( MPI_Comm  *comm_ptr, const char fcname[], int errcode ) {} */
