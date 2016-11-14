#include <iostream>
#include "emdata.h"
#include "transform.h"

using namespace std;
using namespace EMAN;

int main() {
	cout << "Hello" << std::endl;

	Transform * t = new Transform();

	t->printme();

	delete t;

	vector<Transform> vt;

	Transform t1, t2, t3;
	t1.to_identity();
	t2.to_identity();
	t3.to_identity();
	t3.set_trans(7.0,8.0,9.0);

	vt.push_back(t1);
	vt.push_back(t2);
	vt.push_back(t3);

	vt[2].printme();

	EMObject emobj = EMObject(vt);
	vector<Transform> vt2 = emobj;
	vt2[2].printme();

	return 0;
}

