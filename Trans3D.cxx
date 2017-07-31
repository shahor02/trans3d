// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Trans3D.h"

using namespace o2::base;

ClassImp(o2::base::Trans3D)

//_________________________________________________
Trans3D::Trans3D(const TGeoMatrix &m)
{
  /*
   * Construct from TGeoMatrix or its derived classes
   */
  const double* t = m.GetTranslation();
  const double *r  = m.GetRotationMatrix();
  SetComponents(r[0],r[1],r[2],t[0],
		r[3],r[4],r[5],t[1],
		r[6],r[7],r[8],t[2]);
}

//_________________________________________________
void Trans3D::print() const
{
  /*
   * print itself
   */
  std::cout << *this << std::endl;

}
/*
  template <typename T>
  Point3D<T> operator() (const Point3D<T> & p) const {
    /// Local to global transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    for (int i=0;i<12;i++) printf("%+e ",fm[i]); printf("\n");
    return Point3D<T> ( fm[kXX]*p.X() + fm[kXY]*p.Y() + fm[kXZ]*p.Z() + fm[kDX],
			fm[kYX]*p.X() + fm[kYY]*p.Y() + fm[kYZ]*p.Z() + fm[kDY],
			fm[kZX]*p.X() + fm[kZY]*p.Y() + fm[kZZ]*p.Z() + fm[kDZ] );
  }

  template <typename T>
    Point3D<T> operator^ (const Point3D<T> & p) const {
    /// global to local (inverse) transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    for (int i=0;i<12;i++) printf("%+e ",fm[i]); printf("\n");
    Point3D<double> tmp(p.X() - fm[kDX], p.Y() - fm[kDY], p.Z() - fm[kDZ]);
    return Point3D<T> (fm[kXX] * tmp.X() + fm[kYX] * tmp.Y() + fm[kZX] * tmp.Z(),
		       fm[kXY] * tmp.X() + fm[kYY] * tmp.Y() + fm[kZY] * tmp.Z(),
		       fm[kXZ] * tmp.X() + fm[kYZ] * tmp.Y() + fm[kZZ] * tmp.Z());
  }
*/
