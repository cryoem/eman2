// Copyright (C) 2005-2008 Washington University in St Louis, Baylor College of Medicine.  All rights reserved
// Author:        Tao Ju, Refactored by Sasakthi S. Abeysinghe (sasakthi@gmail.com)
// Description:   Grid queue

#include "grid_queue.h"
using namespace wustl_mm::SkeletonMaker;

		GridQueue::GridQueue( ) {
			head = NULL ;
			tail = NULL ;
			numEles = 0 ;
		}

		gridQueueEle* GridQueue::getHead( ) {
			return head ;
		}

		int GridQueue::getNumElements( ) {
			return numEles ;
		}

		/* Naive bubble sort */
		void GridQueue::sort( int eles ) {
			//printf("Sorting elements with descending scores...\n") ;
			gridQueueEle* pre ;
			gridQueueEle* e1 ;
			gridQueueEle* e2 ;

			for ( int i = eles - 1 ; i > 0 ; i -- )
			{
				pre = NULL ;
				e1 = head ;
				e2 = e1->next ;

				for ( int j = 0 ; j < i ; j ++ )
				{
					if ( e1->score < e2->score )
					{
						swapEle( pre, e1, e2 ) ;
					}
					else if ( e1->score == e2->score && rand() < RAND_MAX / 2)
					{
						swapEle( pre, e1, e2 ) ;
					}

					if ( pre == NULL )
					{
						pre = head ;
					}
					else
					{
						pre = pre->next;
					}
					e1 = pre->next ;
					e2 = e1->next ;
				}
			}

			/* Debugging
			pre = head ;
			while ( pre != NULL )
			{
				printf("%d ", pre->score ) ;
				pre = pre->next ;
			}*/
		}


		void GridQueue::pushQueue( int xx, int yy, int zz ) {
			gridQueueEle* ele = new gridQueueEle ;
			ele->x = xx ;
			ele->y = yy ;
			ele->z = zz ;
			ele->score = 0 ;
			ele->next = NULL ;
			if ( head == NULL )
			{
				head = ele ;
			}
			else
			{
				tail->next = ele ;
			}
			tail = ele ;
			numEles ++ ;
		}

		int GridQueue::popQueue( int& xx, int& yy, int& zz ) {
			if ( head == NULL )
			{
				return 0 ;
			}

			xx = head->x ;
			yy = head->y ;
			zz = head->z ;

			gridQueueEle* temp = head ;
			head = head->next ;
			delete temp ;

			if ( head == NULL )
			{
				tail = NULL ;
			}
			numEles -- ;

			return 1 ;
		}

		/* Switching two elements */
		void GridQueue::swapEle( gridQueueEle* pre, gridQueueEle* e1, gridQueueEle* e2 )
		{
			if ( pre != NULL )
			{
				e1->next = e2->next ;
				e2->next = pre->next ;
				pre->next = e2 ;

				if ( tail == e2 )
				{
					tail = e1 ;
				}
			}
			else
			{
				e1->next = e2->next ;
				e2->next = e1 ;
				head = e2 ;

				if ( tail == e2 )
				{
					tail = e1 ;
				}
			}
		}
