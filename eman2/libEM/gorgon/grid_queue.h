// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju, Refactored by Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Grid queue

#ifndef SKELETON_MAKER_GRID_QUEUE_H
#define SKELETON_MAKER_GRID_QUEUE_H

#include <cstdio>
#include <cstdlib>
using namespace std;

namespace wustl_mm {
	namespace SkeletonMaker {
		struct gridQueueEle
		{
			int x, y, z;
			int score ;
			gridQueueEle* next ;
		};

		class GridQueue
		{
		public:
			GridQueue();
			gridQueueEle* getHead();
			int getNumElements();
			void sort(int eles);
			void pushQueue(int xx, int yy, int zz);
			int popQueue(int& xx, int& yy, int& zz);


		private:
			void swapEle(gridQueueEle* pre, gridQueueEle* e1, gridQueueEle* e2);
		private:
			gridQueueEle* head ;
			gridQueueEle* tail ;
			int numEles ;
		};

	}
}
#endif
