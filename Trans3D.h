// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_BASE_TRANS3D_H_
#define ALICEO2_BASE_TRANS3D_H_

#include <Math/GenVector/DisplacementVector3D.h>
#include <Math/GenVector/PositionVector3D.h>
#include <Math/GenVector/Transform3D.h>

template <typename T>
using Point3D = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;
template <typename T>
using Vector3D = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;


namespace o2 {
namespace Base {

class Trans3D : ROOT::Math::Transform3D/*<double>*/ {

 public: 

  using ROOT::Math::Transform3D::Transform3D;
  using ROOT::Math::Transform3D::SetComponents;
  using ROOT::Math::Transform3D::GetComponents;  
  
  template <typename T>
    Point3D<T> operator() (const Point3D<T> & p) const {
    /*
     * Local to global transformation. Once the matrix coefs will be protected, we will use them directly
     */
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Point3D<T> ( fm[kXX]*p.X() + fm[kXY]*p.Y() + fm[kXZ]*p.Z() + fm[kDX],
			fm[kYX]*p.X() + fm[kYY]*p.Y() + fm[kYZ]*p.Z() + fm[kDY],
			fm[kZX]*p.X() + fm[kZY]*p.Y() + fm[kZZ]*p.Z() + fm[kDZ] );
  }

  template <typename T>
    Point3D<T> operator^ (const Point3D<T> & p) const {
    /*
     * global to local (inverse) transformation. Once the matrix coefs will be protected, we will use them directly
     */
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    Point3D<double> tmp(p.X() - fm[kDX], p.Y() - fm[kDY], p.Z() - fm[kDZ]);
    return Point3D<T> (fm[kXX] * tmp.X() + fm[kYX] * tmp.Y() + fm[kZX] * tmp.Z(),
		       fm[kXY] * tmp.X() + fm[kYY] * tmp.Y() + fm[kZY] * tmp.Z(),
		       fm[kXZ] * tmp.X() + fm[kYZ] * tmp.Y() + fm[kZZ] * tmp.Z());
  }
  
  

 
  ClassDef(Trans3D,1)
};

}
}

#endif
