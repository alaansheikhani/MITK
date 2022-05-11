/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef ABLATIONZONE_H
#define ABLATIONZONE_H

namespace mitk
{
struct AblationZone
  {
    itk::Index<3> indexCenter;
    double radius;
    bool IsEqual(AblationZone a) {
      if ((radius == a.radius) && indexCenter[0] == a.indexCenter[0] && indexCenter[1] == a.indexCenter[1] &&
          indexCenter[2] == a.indexCenter[2])
        {return true;}
      else
        {return false;}
    }
  };
}

#endif // ABLATIONZONE_H
