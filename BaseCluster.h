// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_BASE_BASECLUSTER_H
#define ALICEO2_BASE_BASECLUSTER_H
#include <Math/GenVector/DisplacementVector3D.h>
#include <Math/GenVector/PositionVector3D.h>
#include <iostream>
#include <iomanip>
#include <bitset>
#include <ios>

template <typename T>
using Point3D = ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;
template <typename T>
using Vector3D = ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<T>, ROOT::Math::DefaultCoordinateSystemTag>;

namespace o2
{
namespace Base
{ 
// Basic cluster class with X,Y,Z position detector ID information + user fields
// DetectorID should correspond to continuous (no jumps between detector layers
// planes etc.) internal sensor ID within detector 
// Detector specific clusters should be composed by including it as data member
template <typename T>
class BaseCluster 
{
 private:
  
  Point3D<T>    mPos;          // cartesian position
  std::int16_t  mSensorID;     // the sensor id
  std::int8_t   mCount=0;      // user field reserved for counting
  std::uint8_t  mBits=0;       // user field reserved for bit flags
  enum masks_t : std::int32_t {kMisalignMask=(0x1<<7), kUserBitsMask=~kMisalignMask};
  
 public:
  BaseCluster() = default;

  // constructor
  BaseCluster(T x, T y, T z, std::int16_t sensid)
    :mPos(x, y, z), mSensorID(sensid)
  {
  }

  // getting the cartesian coordinates
  T getX()                            const { return mPos.X(); }
  T getY()                            const { return mPos.Y(); }
  T getZ()                            const { return mPos.Z(); }
  Point3D<T>  getPos()                const { return mPos; }
  Point3D<T>& getPos()                      { return mPos; }
  
  // get sensor id
  std::int16_t getSensorID()          const { return mSensorID; }
  // get count field
  std::int8_t getCount()              const { return mCount; }
  // get bit field
  std::uint8_t getBits()              const { return mBits; }
  bool isBitSet(int bit)              const { return mBits & (0xff & (0x1<<bit)); }
  // check special reserved bit to flag cluster misalignment
  bool IsMisaligned()                 const  { return mBits & kMisalignMask;}
  
  // cast to Point3D
  operator Point3D<T>()               const { return mPos; } 
  operator Point3D<T>&()                    { return mPos; } 
  
  // modifiers

  // set sensor id
  void setSensorID(std::int16_t sid)        { mSensorID = sid; }
  // set count field
  void setCount(std::int8_t c)              { mCount = c; }
  // set bit field
  void setBits(std::uint8_t b)              { mBits = b; }
  void setBit(int bit)                      { mBits |= kUserBitsMask & (0x1<<bit); }
  void resetBit(int bit)                    { mBits &= ~(kUserBitsMask & (0x1<<bit)); }

  // set special reserved bit to flag cluster misalignment
  void setMisaligned()                      { mBits |= kMisalignMask;}
  
  // set position
  void setX(T x)                            { mPos.SetX(x); }
  void setY(T y)                            { mPos.SetY(y); }
  void setZ(T z)                            { mPos.SetZ(z); }
  void setXYZ(T x, T y, T z)
  {
    setX(x);
    setY(y);
    setZ(z);
  }
  void setPos(Point3D<T> const &p)          { mPos = p; }
 
  ClassDefNV(BaseCluster, 1);
};

template <class T>
std::ostream &operator<<(std::ostream &os, const BaseCluster<T> &c)
{
   // stream itself
  os << "SId" << std::setw(5) << c.getSensorID() << " ("
     << std::showpos 
     << std::scientific << c.getX() << ","
     << std::scientific << c.getY() << ","
     << std::scientific << c.getZ() << ") cnt:"
     << std::setw(4) << +c.getCount() << " bits:"
     << std::bitset<8>(c.getBits());
  return os;
}

 
}
} // end namespace AliceO2
#endif
