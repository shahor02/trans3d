// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_BASE__TRANS3D_H_
#define ALICEO2_BASE__TRANS3D_H_

#include <Math/GenVector/DisplacementVector3D.h>
#include <Math/GenVector/PositionVector3D.h>
#include <Math/GenVector/Transform3D.h>
#include <TGeoMatrix.h>

template <typename T>
using Point3D = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;
template <typename T>
using Vector3D = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;


namespace o2 {
namespace base {

class Trans3D : public ROOT::Math::Transform3D {

 public: 

  Trans3D() = default;
  Trans3D(const TGeoMatrix &m);
  ~Trans3D() = default;

  template <typename T>
    Point3D<T> operator() (const Point3D<T>& p) const; // local->master

  template <typename T>
    Point3D<T> operator^ (const Point3D<T>& p) const; // master->local

  template <typename T>
    Vector3D<T> operator() (const Vector3D<T> & v) const; // master->local

  template <typename T>
    Vector3D<T> operator^ (const Vector3D<T> & v) const; // local->master

  // TGeoHMatrix-like aliases
  template <typename T>
    void LocalToMaster(const Point3D<T>& loc, Point3D<T>& mst) const {mst = operator()(loc);}

  template <typename T>
    void MasterToLocal(const Point3D<T>& mst, Point3D<T>& loc) const {loc = operator^(loc);}
  
  template <typename T>
    void LocalToMasterVect(const Point3D<T>& loc, Point3D<T>& mst) const {mst = operator()(loc);}

  template <typename T>
    void MasterToLocalVect(const Point3D<T>& mst, Point3D<T>& loc) const {loc = operator^(loc);}

  void print() const;
 
  ClassDef(Trans3D,1)
};

 
 template <typename T>
   Point3D<T> Trans3D::operator() (const Point3D<T> & p) const {
    // Local to master transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    for (int i=0;i<12;i++) printf("%.3e ",fm[i]); printf("\n");
    return Point3D<T> ( fm[kXX]*p.X() + fm[kXY]*p.Y() + fm[kXZ]*p.Z() + fm[kDX],
			fm[kYX]*p.X() + fm[kYY]*p.Y() + fm[kYZ]*p.Z() + fm[kDY],
			fm[kZX]*p.X() + fm[kZY]*p.Y() + fm[kZZ]*p.Z() + fm[kDZ] );
  }

  template <typename T>
    Point3D<T> Trans3D::operator^ (const Point3D<T> & p) const {
      // master to local (inverse) transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    for (int i=0;i<12;i++) printf("%.3e ",fm[i]); printf("\n");
    Point3D<double> tmp(p.X() - fm[kDX], p.Y() - fm[kDY], p.Z() - fm[kDZ]);
    return Point3D<T> (fm[kXX] * tmp.X() + fm[kYX] * tmp.Y() + fm[kZX] * tmp.Z(),
		       fm[kXY] * tmp.X() + fm[kYY] * tmp.Y() + fm[kZY] * tmp.Z(),
		       fm[kXZ] * tmp.X() + fm[kYZ] * tmp.Y() + fm[kZZ] * tmp.Z());
  }

  template <typename T>
    Vector3D<T> Trans3D::operator() (const Vector3D<T> & v) const {
    // Local to master transformation.
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Vector3D<T>( fm[kXX]*v.X() + fm[kXY]*v.Y() + fm[kXZ]*v.Z() ,
		   fm[kYX]*v.X() + fm[kYY]*v.Y() + fm[kYZ]*v.Z() ,
		   fm[kZX]*v.X() + fm[kZY]*v.Y() + fm[kZZ]*v.Z()  );
  }

  template <typename T>
    Vector3D<T> Trans3D::operator^ (const Vector3D<T> & v) const {
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Vector3D<T>(fm[kXX] * v.X() + fm[kYX] * v.Y() + fm[kZX] * v.Z(),
		       fm[kXY] * v.X() + fm[kYY] * v.Y() + fm[kZY] * v.Z(),
		       fm[kXZ] * v.X() + fm[kYZ] * v.Y() + fm[kZZ] * v.Z());
  }
    

/*
  RS: templated version w/o explicit specialization misteriously does not work: 
  GetComponents method called from the templated method fills fm array with shift of 1, 
  while Trans3D::GetComponents(...) returns correct array, i.e.
  using namespace o2::base;
  Trans3D tr; // unit matrix
  Point3D<double> p;
  tr(p);
  //-> 6.907e-310 1.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.000e+00
  double ad[12];
  tr.GetComponents(ad[0],ad[1],ad[2],ad[3],ad[4],ad[5],ad[6],ad[7],ad[8],ad[9],ad[10],ad[11]);
  for (int i=0;i<12;i++) printf("%.3e ",ad[i]); printf("\n");
  //-> 1.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.000e+00 0.000e+00 0.000e+00 0.000e+00 0.000e+00 1.000e+00 0.000e+00

  Therefore we define explicit specializations
*/
/*
  template <>
  Point3D<float> Trans3D::operator() (const Point3D<float> & p) const {
    // Local to master transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Point3D<float> ( fm[kXX]*p.X() + fm[kXY]*p.Y() + fm[kXZ]*p.Z() + fm[kDX],
			    fm[kYX]*p.X() + fm[kYY]*p.Y() + fm[kYZ]*p.Z() + fm[kDY],
			    fm[kZX]*p.X() + fm[kZY]*p.Y() + fm[kZZ]*p.Z() + fm[kDZ] );
  }

  template <>
  Point3D<double> Trans3D::operator() (const Point3D<double> & p) const {
    // Local to master transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Point3D<double> ( fm[kXX]*p.X() + fm[kXY]*p.Y() + fm[kXZ]*p.Z() + fm[kDX],
			     fm[kYX]*p.X() + fm[kYY]*p.Y() + fm[kYZ]*p.Z() + fm[kDY],
			     fm[kZX]*p.X() + fm[kZY]*p.Y() + fm[kZZ]*p.Z() + fm[kDZ] );
  }

  template <>
  Point3D<float> Trans3D::operator^ (const Point3D<float> & p) const {
      // master to local (inverse) transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    Point3D<float> tmp(p.X() - fm[kDX], p.Y() - fm[kDY], p.Z() - fm[kDZ]);
    return Point3D<float> (fm[kXX] * tmp.X() + fm[kYX] * tmp.Y() + fm[kZX] * tmp.Z(),
			   fm[kXY] * tmp.X() + fm[kYY] * tmp.Y() + fm[kZY] * tmp.Z(),
			   fm[kXZ] * tmp.X() + fm[kYZ] * tmp.Y() + fm[kZZ] * tmp.Z());
  }

  template <>
  Point3D<double> Trans3D::operator^ (const Point3D<double> & p) const {
      // master to local (inverse) transformation. Once the matrix coefs will be protected, we will use them directly
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    Point3D<double> tmp(p.X() - fm[kDX], p.Y() - fm[kDY], p.Z() - fm[kDZ]);
    return Point3D<double> (fm[kXX] * tmp.X() + fm[kYX] * tmp.Y() + fm[kZX] * tmp.Z(),
			    fm[kXY] * tmp.X() + fm[kYY] * tmp.Y() + fm[kZY] * tmp.Z(),
			    fm[kXZ] * tmp.X() + fm[kYZ] * tmp.Y() + fm[kZZ] * tmp.Z());
  }

  template <>
    Vector3D<float> Trans3D::operator() (const Vector3D<float> & v) const {
    // Local to master transformation.
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Vector3D<float>( fm[kXX]*v.X() + fm[kXY]*v.Y() + fm[kXZ]*v.Z() ,
			    fm[kYX]*v.X() + fm[kYY]*v.Y() + fm[kYZ]*v.Z() ,
			    fm[kZX]*v.X() + fm[kZY]*v.Y() + fm[kZZ]*v.Z()  );
  }

  template <>
    Vector3D<float> Trans3D::operator^ (const Vector3D<float> & v) const {
    // Master to local transformation.
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Vector3D<float>(fm[kXX] * v.X() + fm[kYX] * v.Y() + fm[kZX] * v.Z(),
			   fm[kXY] * v.X() + fm[kYY] * v.Y() + fm[kZY] * v.Z(),
			   fm[kXZ] * v.X() + fm[kYZ] * v.Y() + fm[kZZ] * v.Z());
  }

  template <>
    Vector3D<double> Trans3D::operator() (const Vector3D<double> & v) const {
    // Local to master transformation.
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Vector3D<double>( fm[kXX]*v.X() + fm[kXY]*v.Y() + fm[kXZ]*v.Z() ,
			    fm[kYX]*v.X() + fm[kYY]*v.Y() + fm[kYZ]*v.Z() ,
			    fm[kZX]*v.X() + fm[kZY]*v.Y() + fm[kZZ]*v.Z()  );
  }

  template <>
    Vector3D<double> Trans3D::operator^ (const Vector3D<double> & v) const {
    // Master to local transformation.
    double fm[12];
    GetComponents(fm[kXX],fm[kXY],fm[kXZ],fm[kDX],fm[kYX],fm[kYY],fm[kYZ],fm[kDY],fm[kZX],fm[kZY],fm[kZZ],fm[kDZ]);
    return Vector3D<double>(fm[kXX] * v.X() + fm[kYX] * v.Y() + fm[kZX] * v.Z(),
			    fm[kXY] * v.X() + fm[kYY] * v.Y() + fm[kZY] * v.Z(),
			    fm[kXZ] * v.X() + fm[kYZ] * v.Y() + fm[kZZ] * v.Z());
  }

*/
  
}
}

#endif
