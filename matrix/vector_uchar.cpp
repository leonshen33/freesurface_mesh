/*=============================================================================
        File: vector.cpp
     Purpose:
    Revision: $Id: vector_uchar.cpp,v 1.2 2002/05/13 21:07:45 philosophil Exp $
  Created by:    Philippe Lavoie          (3 Oct, 1996)
 Modified by: 

 Copyright notice:
          Copyright (C) 1996-1998 Philippe Lavoie
 
	  This library is free software; you can redistribute it and/or
	  modify it under the terms of the GNU Library General Public
	  License as published by the Free Software Foundation; either
	  version 2 of the License, or (at your option) any later version.
 
	  This library is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	  Library General Public License for more details.
 
	  You should have received a copy of the GNU Library General Public
	  License along with this library; if not, write to the Free
	  Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
=============================================================================*/
#ifndef vector_uchar_h__
#define vector_uchar_h__



#include "vector.h"

namespace PLib {

#ifdef NO_IMPLICIT_TEMPLATES

  template class Vector<unsigned char> ;
  
  template Vector<unsigned char> operator+(const Vector<unsigned char>&, const Vector<unsigned char>&);
  template Vector<unsigned char> operator-(const Vector<unsigned char>&, const Vector<unsigned char>&);
  template unsigned char operator*(const Vector<unsigned char>&,const Vector<unsigned char>&);
  template Vector<unsigned char> operator*(const Vector<unsigned char>& v, const double d);
  template Vector<unsigned char> operator*(const Vector<unsigned char>& v, const Complex d);
  template int operator==(const Vector<unsigned char>&,const Vector<unsigned char>&);
  template int operator!=(const Vector<unsigned char>&,const Vector<unsigned char>&);
  
#endif

}

#endif // vector_uchar_h__