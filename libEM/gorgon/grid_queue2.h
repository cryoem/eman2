// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju, Refactored by Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Grid queue

#ifndef SKELETON_MAKER_GRID_QUEUE2_H
#define SKELETON_MAKER_GRID_QUEUE2_H

#include "grid_queue.h"

namespace wustl_mm {
	namespace SkeletonMaker {

		class GridQueue2
		{
		public:
			GridQueue2();
			~GridQueue2();
			gridQueueEle* getNext();
			void reset();
			int getNumElements();
			void prepend(int xx, int yy, int zz);
			gridQueueEle * remove();
			gridQueueEle* swap();
		private:
			gridQueueEle* head ;
			gridQueueEle* pre ;
			gridQueueEle* prepre ;
			gridQueueEle* cur ;
			int numEles ;
		};

	}
}
#endif
